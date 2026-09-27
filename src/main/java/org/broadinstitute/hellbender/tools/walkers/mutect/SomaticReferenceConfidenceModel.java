package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceResult;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public class SomaticReferenceConfidenceModel extends ReferenceConfidenceModel {

    private final SampleList samples;
    private final Optional<BetaDistributionShape> afPrior;



    /**
     * Create a new ReferenceConfidenceModel
     *  @param samples the list of all samples we'll be considering with this model
     * @param header the SAMFileHeader describing the read information (used for debugging)
     * @param indelInformativeDepthIndelSize the max size of indels to consider when calculating indel informative depths
     * @param minAF soft threshold for allele fractions -- above this value prior is nearly flat, below, prior is nearly zero
     * @param refModelDelQual reference model deletion quality (if constant)
     * @param useSoftClippedBases should soft clipped bases be counted against the reference
     * @param isFlowBasedModel is the error model flow based
     */
    SomaticReferenceConfidenceModel(final SampleList samples, final SAMFileHeader header, final int indelInformativeDepthIndelSize,
                                    final double minAF, final byte refModelDelQual, final boolean useSoftClippedBases,
                                    final boolean isFlowBasedModel){
        super(samples, header, indelInformativeDepthIndelSize, 0, refModelDelQual, useSoftClippedBases, isFlowBasedModel);
        Utils.validateArg(minAF >= 0.0 && minAF < 1, "minAF must be < 1 and >= 0");

        // To softly cut off allele fractions below minAF, we use a Beta prior of the form Beta(1+epsilon, 1); that is
        // the prior on allele fraction f is proportional to f^epsilon.  If epsilon is small this prior vanishes as f -> 0
        // and very rapidly becomes flat.  We choose epsilon such that minAF^epsilon = 0.5.
        afPrior = minAF == 0.0 ? Optional.empty() : Optional.of(new BetaDistributionShape(1 - Math.log(2)/Math.log(minAF), 1));
        this.samples = samples;
    }

        /**
         * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
         *
         * @param ploidy target sample ploidy.
         * @param pileup the read backed pileup containing the data we want to evaluate
         * @param refBase the reference base at this pileup position
         * @param minBaseQual the min base quality for a read in the pileup at the pileup position to be included in the calculation
         * @param hqSoftClips running average data structure (can be null) to collect information about the number of high quality soft clips
         * @return a RefVsAnyResult genotype call.
         */
    @Override
    public ReferenceConfidenceResult calcGenotypeLikelihoodsOfRefVsAny(final int ploidy,
                                                                       final ReadPileup pileup,
                                                                       final byte refBase,
                                                                       final byte minBaseQual,
                                                                       final MathUtils.RunningAverage hqSoftClips,
                                                                       final boolean readsWereRealigned) {

        final SomaticRefVsAnyResult result = new SomaticRefVsAnyResult();
        final Map<String, List<GATKRead>> perSampleReadMap = new HashMap<>();
        perSampleReadMap.put(samples.getSample(0), pileup.getReads());

        final List<Byte> altQuals = new ArrayList<>(pileup.size() / 20);

        for (final PileupElement element : pileup) {
            if (!element.isDeletion() && element.getQual() <= minBaseQual) {
                continue;
            }

            final boolean isAlt = readsWereRealigned ? isAltAfterAssembly(element, refBase) : isAltBeforeAssembly(element, refBase);
            if (isAlt) {
                altQuals.add(element.getQual());
                result.nonRefDepth++;
            } else {
                result.refDepth++;
            }
        }

        final double logOdds = Mutect2Engine.logLikelihoodRatio(result.refDepth, altQuals, 1, afPrior);
        result.lods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        result.lods.set(Allele.NON_REF_ALLELE, logOdds);

        return result;
    }

    @Override
    public void addGenotypeData(final ReferenceConfidenceResult result, final GenotypeBuilder gb) {
        gb.attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, MathUtils.logToLog10(((SomaticRefVsAnyResult)result).lods.get(Allele.NON_REF_ALLELE)));
    }

    /**
     * Evaluates each position from its pileup. The somatic model rates the SNP evidence alone: the germline indel
     * model counts reads that rule out indels of up to 10bp and folds a constant-quality diploid indel likelihood into
     * the site, which does not represent low allele fraction variants or PCR errors well, and the first use of this
     * model is mitochondrial calling, where the high-complexity reference lets the SNP model dominate.
     */
    @Override
    protected List<ReferenceConfidenceResult> calculateSiteResults(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                   final AssemblyRegion activeRegion,
                                                                   final byte[] ref,
                                                                   final int ploidy) {
        final SimpleInterval span = activeRegion.getSpan();
        final int globalRefOffset = span.getStart() - activeRegion.getPaddedSpan().getStart();
        final List<ReadPileup> pileups = AssemblyBasedCallerUtils.getPileupsOverReference(activeRegion.getHeader(), span, readLikelihoods, samples);
        final List<ReferenceConfidenceResult> results = new ArrayList<>(pileups.size());
        for (int i = 0; i < pileups.size(); i++) {
            results.add(calcGenotypeLikelihoodsOfRefVsAny(ploidy, pileups.get(i), ref[i + globalRefOffset], BASE_QUAL_THRESHOLD, null, true));
        }
        return results;
    }

}
