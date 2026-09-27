package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingDebugger;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Helper to calculate genotype likelihoods for DRAGEN advanced genotyping models (BQD - Base Quality Dropout, and FRD - Foreign Reads Detection).
 *
 * This object is simply a thin wrapper on top of a regular GenotypeLikelihoods object with some extra logic for handling new inputs to the genotyper:
 *  - both BQD and FRD rely on per-read per-genotype scores as would be computed for the standard genotyper, rather than pay the cost of recomputing these
 *    for each of the 3 independent models this GenotypeLikelihoodCalculator simply makes the computation once and relies on the fact that the underlying
 *    readLikelihoodsByGenotypeIndex is still populated from the previous call. To this end strict object equality tests have been implemented to ensure
 *    that the cache is populated with the correct likelihoods before running either of the advanced models.
 */
public final class GenotypeLikelihoodCalculatorDRAGEN extends GenotypeLikelihoodCalculator {
    // measure of the likelihood of an error base occurring if we have triggered a base quality dropout
    static final double BQD_FIXED_ERROR_RATE = 0.5;

    // PhredScaled adjustment applied to the BQD score (this controls the weight of the base quality prior term in the BQD calculation)
    static final double PHRED_SCALED_ADJUSTMENT_FOR_BQ_SCORE = 2.5;

    private static final double CACHED_LOG_10_ERROR_RATE = Math.log10(BQD_FIXED_ERROR_RATE);
    private static final double CACHED_LOG_10_NON_ERROR_RATE = Math.log10(1 - BQD_FIXED_ERROR_RATE);

    // FRD's strand models: forward reads only, reverse reads only, and all reads
    private static final int FORWARD_STRAND_MODEL = 0;
    private static final int REVERSE_STRAND_MODEL = 1;
    private static final int BOTH_STRANDS_MODEL = 2;
    private static final int STRAND_MODEL_COUNT = 3;
    private static final String[] STRAND_MODEL_DEBUG_HEADERS = {"\nForwards Strands: ", "\nReverse Strands: ", "\nBoth Strands: "};

    private GenotypeLikelihoodCalculatorDRAGEN() {
        super();
    }

    /**
     * Calculate the FRD model outputs to the likelihoods array.
     * This method handles splitting the model by strand and selecting the best scoring parameters across the two for return in the likelihoods array.
     *
     * BQD needs to see reads that have been disqualified in {@link org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods#filterPoorlyModeledEvidence(ToDoubleFunction)} as
     * well as reads that only overlap the variant in question in their low quality ends. Reads in the former category do not have their hmm scores
     * accounted for in the genotyping model, whereas reads in the later category do. All reads, (disqualified, low quality ends, and all others)
     * are sorted by the cycle-count of the SNP being genotyped, and the average base qualities for the SNP base are computed across partitions
     * of the reads in aggregate.
     *
     * NOTES:
     * - The model will not handle indel alleles
     * - The model currently does not support mixed-allele mode (modifying 0/1 GTs in addition to 0/0 GTs)
     *
     * @param sampleLikelihoods allele likelihoods containing data for reads
     * @param strandForward list of reads in the forwards orientation overlapping the site
     * @param strandReverse list of reads in the reverse orientation overlapping the site
     * @param paddedReference reference bases (with padding) used for calculating homopolymer adjustemnt
     * @param offsetForRefIntoEvent offset of the variant into the reference event
     * @return An array corresponding to the likelihoods array score for BQD, with Double.NEGATIVE_INFINITY filling all mixed allele/indel allelse
     */
    public static <A extends Allele> double[] calculateBQDLikelihoods(final int ploidy, final LikelihoodMatrix<GATKRead, A> sampleLikelihoods,
                                                               final List<DRAGENGenotypesModel.DragenReadContainer> strandForward,
                                                               final List<DRAGENGenotypesModel.DragenReadContainer> strandReverse,
                                                               final byte[] paddedReference,
                                                               final int offsetForRefIntoEvent) {
        final int alleleCount = sampleLikelihoods.numberOfAlleles();
        final double[] outputArray = new double[GenotypeIndexCalculator.genotypeCount(ploidy, alleleCount)];
        Arrays.fill(outputArray, Double.NEGATIVE_INFINITY);

        final Allele refAllele = sampleLikelihoods.getAllele(0);

        for (int gtAlleleIndex = 0; gtAlleleIndex < sampleLikelihoods.numberOfAlleles(); gtAlleleIndex++) {
            // find the index of the homozygous gtAllele genotype
            final int indexForGT = GenotypeIndexCalculator.alleleCountsToIndex(gtAlleleIndex, ploidy);

            for (int errorAlleleIndex = 0; errorAlleleIndex < sampleLikelihoods.numberOfAlleles(); errorAlleleIndex++) {
                // We only want to make calls on SNPs for now
                if (sampleLikelihoods.getAllele(gtAlleleIndex) == sampleLikelihoods.getAllele(errorAlleleIndex) ||
                        sampleLikelihoods.getAllele(gtAlleleIndex).length() != refAllele.length() ||
                        sampleLikelihoods.getAllele(errorAlleleIndex).length() != refAllele.length()) {
                    continue;
                }
                // TODO super validate this
                final byte baseOfErrorAllele = sampleLikelihoods.getAllele(errorAlleleIndex).getBases()[0];

                final double forwardHomopolymerAdjustment = FRDBQDUtils.computeForwardHomopolymerAdjustment(paddedReference, offsetForRefIntoEvent, baseOfErrorAllele);
                final double reverseHomopolymerAdjustment = FRDBQDUtils.computeReverseHomopolymerAdjustment(paddedReference, offsetForRefIntoEvent, baseOfErrorAllele);

                // BQD scores by strand
                final double minScoreFoundForwardsStrand = computeBQDModelForStrandData(sampleLikelihoods, strandForward, forwardHomopolymerAdjustment, true, gtAlleleIndex, errorAlleleIndex);
                final double minScoreFoundReverseStrand = computeBQDModelForStrandData(sampleLikelihoods, strandReverse, reverseHomopolymerAdjustment, false, gtAlleleIndex, errorAlleleIndex);

                final double modelScoreInLog10 = (minScoreFoundForwardsStrand + minScoreFoundReverseStrand) * -0.1;
                //////
                // NOTE we have not applied the prior here, this is because that gets applied downstream to each genotype in the array.
                // since the prior is applied evenly to the defualt and error models this should not change this selection here.
                outputArray[indexForGT] = Math.max(outputArray[indexForGT], modelScoreInLog10);
            }
        }
        return outputArray;
    }

    /**
     * Helper function that actually manages the math for BQD;
     *
     * This method works by combining the computed genotype scores for reads with the raw allele likelihoods scores for each evidence
     *
     * @param sampleLikelihoods allele likelihoods containing data for reads
     * @param positionSortedReads  Reads pairs objects (Pair<Pair<read,readBaseOffset>, sampleReadIndex>) objects sorted in the correct order for partitioning.
     *                             This means that the "error" reads in the partition are sorted by read cycle first in the provided list
     * @param homopolymerAdjustment  Penalty to be applied to reads based on the homopolymer run (this should be precomputed for the ref site in quesiton)
     * @return phred scale likelihood for a BQD error mode for reads in the given direction according to the offsets requested
     */
    private static <A extends Allele> double computeBQDModelForStrandData(final LikelihoodMatrix<GATKRead, A> sampleLikelihoods,
                                                                   final List<DRAGENGenotypesModel.DragenReadContainer> positionSortedReads,
                                                                   final double homopolymerAdjustment,
                                                                   final boolean forwards, final int homozygousAlleleIndex, final int errorAlleleIndex) {
        // If no reads are found for a particular strand direction return no adjusted likelihoods score for those (non-existent) reads
        if (positionSortedReads.isEmpty()) {
            return 0.0;
        }
        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println("errorAllele index: " + errorAlleleIndex + " theta: " + (forwards ? "1" : "2") + " homopolymerAdjustment: " + homopolymerAdjustment);
        }

        // Forwards strand tables (all in phred space for the sake of conveient debugging with provided scripts)
        final int evidenceSize = positionSortedReads.size();
        final double[] cumulativeProbReadForErrorAllele = new double[evidenceSize + 1];
        final double[] cumulativeMeanBaseQualityPhredAdjusted = new double[evidenceSize + 1];
        final double[] cumulativeProbGenotype = new double[evidenceSize + 1];

        double totalBaseQuality = 0;
        int baseQualityDenominator = 0; // We track this separately because not every read overlaps the SNP in question due to padding.
        // Iterate over the reads and populate the cumulative arrays
        for (int i = 1; i < cumulativeProbReadForErrorAllele.length; i++) {
            final DRAGENGenotypesModel.DragenReadContainer container = positionSortedReads.get(i - 1);
            final int readIndex = container.getIndexInLikelihoodsObject();

            // Retrieve the homozygous genotype score and the error allele scores for the read in question
            final double homozygousGenotypeContribution;
            final double errorAlleleContribution;
            if (readIndex != -1) {
                homozygousGenotypeContribution = sampleLikelihoods.get(homozygousAlleleIndex, readIndex);
                errorAlleleContribution = sampleLikelihoods.get(errorAlleleIndex, readIndex);
            } else {
                // If read index == -1 then we are evaluating a read that was rejected by the HMM and therefore doesn't have genotype scores
                homozygousGenotypeContribution = 0;
                errorAlleleContribution = 0;
            }

            // Populate the error probability array in phred space
            // Calculation: Alpha * P(r|E_allele) + (1 - Alpha) * P(r | G_homozygousGT))
            double phredContributionForRead = (homozygousGenotypeContribution==0 && errorAlleleContribution==0) ? 0 : -10 *
                    MathUtils.approximateLog10SumLog10(errorAlleleContribution + CACHED_LOG_10_ERROR_RATE,
                                                       homozygousGenotypeContribution + CACHED_LOG_10_NON_ERROR_RATE);
            cumulativeProbReadForErrorAllele[i] = cumulativeProbReadForErrorAllele[i-1] + phredContributionForRead;

            // Populate the cumulative genotype contribution array with the score for this read
            // Calculation: (P(r | G_A1) + P(r | G_A2)) / 2
            cumulativeProbGenotype[i] = cumulativeProbGenotype[i - 1] + -10 * homozygousGenotypeContribution;

            // Populate the mean base quality array
            if (container.hasValidBaseQuality()) {
                totalBaseQuality += container.getBaseQuality();
                baseQualityDenominator++;
            }
            cumulativeMeanBaseQualityPhredAdjusted[i] = Math.max(0,
                    ((totalBaseQuality / (baseQualityDenominator == 0 ? 1 : baseQualityDenominator)) * PHRED_SCALED_ADJUSTMENT_FOR_BQ_SCORE) - homopolymerAdjustment);
        }

        // Now we find the best partitioning N for the forwards evaluation of the data
        double minScoreFound = Double.POSITIVE_INFINITY;
        int nIndexUsed = 0;
        for (int n = 0; n < cumulativeMeanBaseQualityPhredAdjusted.length; n++) {
            final double bqdScore = cumulativeMeanBaseQualityPhredAdjusted[n] + cumulativeProbReadForErrorAllele[n] + (cumulativeProbGenotype[cumulativeProbGenotype.length-1] - cumulativeProbGenotype[n]);
            if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                HaplotypeCallerGenotypingDebugger.println(String.format("n=%d: %.2f, cum_phred_bq=%.2f, cum_prob_r_Error=%.2f, prob_G_remaining=%.2f",
                        n, bqdScore, cumulativeMeanBaseQualityPhredAdjusted[n], cumulativeProbReadForErrorAllele[n],
                        (cumulativeProbGenotype[cumulativeProbGenotype.length - 1] - cumulativeProbGenotype[n])));
            }
            if (minScoreFound > bqdScore) {
                minScoreFound = bqdScore;
                nIndexUsed = n;
            }
        }

        // Debug output for the genotyper to see into the calculation itself
        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println(String.format("theta=%d n%d=%2d, best_phred_score =%5.2f q_mean=%5.2f, alpha=%4.2f, Ph(E)=%4.2f;  ", forwards ? 1 : 0,
                    (forwards ? 1 : 2), nIndexUsed, minScoreFound, cumulativeMeanBaseQualityPhredAdjusted[nIndexUsed],
                    BQD_FIXED_ERROR_RATE, cumulativeProbReadForErrorAllele[nIndexUsed]));
        }

        return minScoreFound;
    }

    /**
     * Calculate the BQD model outputs to the likelihoods array.
     *
     * This method is responsible for computing critical phred-mapping quality adjustments for the entire pool of reads (Disqualified reads,
     * reads only overlapping in low quality ends, and otherwise) and selecting true-allele/error-allele combinations as well as strand model
     * combinations (all forward reads/ all reverse reads/ all reads), calling {@link #computeFRDModelsForStrands} for each
     * true-allele/error-allele combination and selecting the best scoring columns in the final likelihoods array output.
     *
     * Like BQD this model genotypes with all reads that overlap the site in either their accepted bases or low quality ends, but it does
     * not include disqualified reads for genotyping. All reads are used for computing the critical values for the mapping quality cutoffs.
     *
     * NOTES:
     * - The model currently does not support mixed-allele mode (modifying 0/1 GTs in addition to 0/0 GTs)
     * - The model will not treat symbolic alleles specially, always treating them as indels. This might or might not be the best way to handle them.
     *
     * @param sampleLikelihoods the likelihoods object with allele likelihoods for the reads to be genotyped
     * @param ploidyModelLikelihoods standard genotyping model allele likelihoods (to be used for maxEffectiveDepthAdjustment)
     * @param readContainers reads (both forwards and reverse orientation as well as disqualified reads) overlapping the site in question
     * @param snipAprioriHet prior for heterozygus SNP allele
     * @param indelAprioriHet prior for heterozygus indel alleles based on the STRE tables if present
     * @param maxEffectiveDepthForHetAdjustment maxEffectiveDepthAdjustment used to reduce the effect of FRD at high depth sites (0 means no adjustment)
     * @return a likelihoods array corrsponding to the log10 likelihoods scores for the best combination of model parameters for each Genotype (Double.NEGATIVE_INFINITY for Genotypes not considered)
     */
    public static <A extends Allele> double[] calculateFRDLikelihoods(final int ploidy, final LikelihoodMatrix<GATKRead, A> sampleLikelihoods, final double[] ploidyModelLikelihoods,
                                                               final List<DRAGENGenotypesModel.DragenReadContainer> readContainers,
                                                               final double snipAprioriHet, final double indelAprioriHet, final int maxEffectiveDepthForHetAdjustment) {
        final int alleleCount = sampleLikelihoods.numberOfAlleles();
        final double[] outputArray = new double[GenotypeIndexCalculator.genotypeCount(ploidy, alleleCount)];
        Arrays.fill(outputArray, Double.NEGATIVE_INFINITY);

        final Allele refAllele = sampleLikelihoods.getAllele(0);
        final FRDReads reads = new FRDReads(sampleLikelihoods, readContainers);

        for (int fAlleleIndex = 0; fAlleleIndex < alleleCount; fAlleleIndex++) {
            // ignore symbolic alleles
            final boolean isIndel = sampleLikelihoods.getAllele(fAlleleIndex).length() != refAllele.length();

            // Here we generate a set of the critical log10(P(F)) values that we will iterate over
            final double[] criticalThresholds = reads.updateCriticalValues(fAlleleIndex == 0 ? 0 : (isIndel? indelAprioriHet : snipAprioriHet) * -0.1); // simplified in line with DRAGEN, uses 1 alleledist for both snp and indels

            if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                HaplotypeCallerGenotypingDebugger.println("fIndex: " + fAlleleIndex + " criticalValues: \n" + Arrays.stream(criticalThresholds).mapToObj(Double::toString).collect(Collectors.joining("\n")));
            }
            // iterate over all of the homozygous genotypes for the given allele
            for (int gtAlleleIndex = 0; gtAlleleIndex < alleleCount; gtAlleleIndex++) {
                // Skip over the allele corresponding to the "foreign" allele
                if (gtAlleleIndex == fAlleleIndex) {
                    continue;
                }

                //This is crufty, it just so happens that the index of the homozygous genotype corresponds to the maximum genotype count per field.
                //This should be pulled off as a calculator in some genotyping class.
                final int indexForGT = GenotypeIndexCalculator.alleleCountsToIndex(gtAlleleIndex, ploidy);

                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                    HaplotypeCallerGenotypingDebugger.println("indexForGT "+indexForGT);
                }
                final double[][] strandModels = computeFRDModelsForStrands(reads, gtAlleleIndex, fAlleleIndex, criticalThresholds);
                final double[] maxLog10FForwardsStrand = strandModels[FORWARD_STRAND_MODEL];
                final double[] maxLog10FReverseStrand = strandModels[REVERSE_STRAND_MODEL];
                final double[] maxLog10FBothStrands = strandModels[BOTH_STRANDS_MODEL];

                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                    HaplotypeCallerGenotypingDebugger.println("gtAlleleIndex : "+gtAlleleIndex+ " fAlleleIndex: "+fAlleleIndex +" forwards: "+Arrays.toString(maxLog10FForwardsStrand)+" reverse: "+Arrays.toString(maxLog10FReverseStrand)+" both: "+Arrays.toString(maxLog10FBothStrands));
                }
                double[] localBestModel = maxLog10FForwardsStrand;
                if (localBestModel[0] < maxLog10FReverseStrand[0]) {
                    localBestModel = maxLog10FReverseStrand;
                }
                if (localBestModel[0] < maxLog10FBothStrands[0]) {
                    localBestModel = maxLog10FBothStrands;
                }

                // Handle max effective depth adjustment if specified
                if (maxEffectiveDepthForHetAdjustment > 0) {
                    // Use the index corresponding the mixture of F and
                    final double localBestModelScore = localBestModel[0] - localBestModel[1];
                    final int closestGTAlleleIndex = GenotypeIndexCalculator.allelesToIndex(gtAlleleIndex, fAlleleIndex);
                    final double log10LikelihoodsForPloyidyModel = ploidyModelLikelihoods[closestGTAlleleIndex] - -MathUtils.LOG10_ONE_HALF;
                    final int depthForGenotyping = sampleLikelihoods.evidenceCount();
                    final double adjustedBestModel = log10LikelihoodsForPloyidyModel + ((localBestModelScore - log10LikelihoodsForPloyidyModel)
                            * ((Math.min(depthForGenotyping, maxEffectiveDepthForHetAdjustment) * 1.0) / depthForGenotyping));
                    outputArray[indexForGT] = Math.max(outputArray[indexForGT], adjustedBestModel + localBestModel[1]);

                    if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                        HaplotypeCallerGenotypingDebugger.println("best FRD likelihoods: "+localBestModelScore+" P(F) score used: "+localBestModel[1]+"  use MaxEffectiveDepth: "+maxEffectiveDepthForHetAdjustment);
                        HaplotypeCallerGenotypingDebugger.println("Using array index "+closestGTAlleleIndex+" for mixture gt with likelihood of "+log10LikelihoodsForPloyidyModel+" adjusted based on depth: "+depthForGenotyping);
                        HaplotypeCallerGenotypingDebugger.println("p_rG_adj : "+adjustedBestModel);
                    }
                } else {
                    outputArray[indexForGT] = Math.max(outputArray[indexForGT], localBestModel[0]);
                }

            }

        }


        return outputArray;
    }

    /**
     * Computes the FRD model for one homozygous genotype and foreign allele under each strand combination: forward
     * reads only, reverse reads only, and all reads. Reads outside a model's strand still contribute their genotype
     * likelihood to it. Every critical threshold is evaluated for every model, as DRAGEN does, even where a threshold
     * only arises from reads on the other strand.
     *
     * The three models share the two passes over the reads made for each threshold; each model's sums still add its
     * reads in read-container order, so the results equal those of separate passes per model.
     *
     * @param reads the reads to genotype, with critical values set for {@code fAlleleIndex}
     * @param homozygousAlleleIndex index of allele in homzygous genotype whose likelihood is to be adjusted
     * @param fAlleleIndex index of foreign allele within likelihoods matrix
     * @param criticalThresholdsSorted distinct critical thresholds in ascending order
     * @return for each strand model (indexed by {@link #FORWARD_STRAND_MODEL}, {@link #REVERSE_STRAND_MODEL} and
     *         {@link #BOTH_STRANDS_MODEL}), two doubles: index 0 is the frd score and the second is log p(F()) score
     *         used to adjust the score
     */
    private static double[][] computeFRDModelsForStrands(final FRDReads reads, final int homozygousAlleleIndex, final int fAlleleIndex,
                                                         final double[] criticalThresholdsSorted) {
        if (!reads.anyReads) {
            if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                Arrays.stream(STRAND_MODEL_DEBUG_HEADERS).forEach(HaplotypeCallerGenotypingDebugger::println);
            }
            return new double[][]{{Double.NEGATIVE_INFINITY, 0}, {Double.NEGATIVE_INFINITY, 0}, {Double.NEGATIVE_INFINITY, 0}};
        }

        final int readCount = reads.genotypedCount;
        final int forwardCount = reads.forwardCount;
        final boolean[] isReverseStrand = reads.isReverseStrand;
        final double[] criticalValues = reads.criticalValuesWithTolerance;
        final double[] log10LikelihoodsForF = reads.log10LikelihoodsByAllele[fAlleleIndex];
        final double[] log10LikelihoodsForGT = reads.log10LikelihoodsByAllele[homozygousAlleleIndex];

        // A read's support for the foreign allele depends on the threshold only through whether the threshold
        // excludes the read's foreign-allele likelihood, so both possible values are computed once per read. The
        // excluded value is 0.0 unless the genotype likelihood is -Infinity; it is computed rather than written as a
        // constant so that the NaN of that corner case is exactly the value the formula gives.
        final double[] supportIfIncluded = new double[readCount];
        final double[] supportIfExcluded = new double[readCount];
        for (int i = 0; i < readCount; i++) {
            supportIfIncluded[i] = Math.pow(10, log10LikelihoodsForF[i] - MathUtils.approximateLog10SumLog10(log10LikelihoodsForF[i], log10LikelihoodsForGT[i]));
            supportIfExcluded[i] = Math.pow(10, Double.NEGATIVE_INFINITY - MathUtils.approximateLog10SumLog10(Double.NEGATIVE_INFINITY, log10LikelihoodsForGT[i]));
        }

        final double[] maxLpspi = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
        final double[] lpfApplied = new double[STRAND_MODEL_COUNT];
        final double[] foreignAlleleLikelihoods = new double[STRAND_MODEL_COUNT];
        final double[] cumulativeLikelihoods = new double[STRAND_MODEL_COUNT];
        final List<List<String>> debugLines = HaplotypeCallerGenotypingDebugger.isEnabled() ?
                List.of(new ArrayList<>(), new ArrayList<>(), new ArrayList<>()) : null;

        for (int thresholdIndex = 0; thresholdIndex < criticalThresholdsSorted.length; thresholdIndex++) {
            final double logProbFAllele = criticalThresholdsSorted[thresholdIndex];

            // the foreign allele's alpha for each model: the mean support for the foreign allele over the model's
            // reads, where a read whose critical value is at or below the threshold gives no support
            double fAlleleProbRatioForward = 0.0;
            double fAlleleProbRatioReverse = 0.0;
            double fAlleleProbRatioBoth = 0.0;
            for (int i = 0; i < readCount; i++) {
                final double support = criticalValues[i] <= logProbFAllele ? supportIfExcluded[i] : supportIfIncluded[i];
                fAlleleProbRatioBoth += support;
                if (isReverseStrand[i]) {
                    fAlleleProbRatioReverse += support;
                } else {
                    fAlleleProbRatioForward += support;
                }
            }

            // Don't learn the beta but approximate it based on the read support for the alt
            foreignAlleleLikelihoods[FORWARD_STRAND_MODEL] = Math.min(fAlleleProbRatioForward / forwardCount, 0.5);
            foreignAlleleLikelihoods[REVERSE_STRAND_MODEL] = Math.min(fAlleleProbRatioReverse / (readCount - forwardCount), 0.5);
            foreignAlleleLikelihoods[BOTH_STRANDS_MODEL] = Math.min(fAlleleProbRatioBoth / readCount, 0.5);
            final double log10ForeignForward = Math.log10(foreignAlleleLikelihoods[FORWARD_STRAND_MODEL]);
            final double log10NotForeignForward = Math.log10(1.0 - foreignAlleleLikelihoods[FORWARD_STRAND_MODEL]);
            final double log10ForeignReverse = Math.log10(foreignAlleleLikelihoods[REVERSE_STRAND_MODEL]);
            final double log10NotForeignReverse = Math.log10(1.0 - foreignAlleleLikelihoods[REVERSE_STRAND_MODEL]);
            final double log10ForeignBoth = Math.log10(foreignAlleleLikelihoods[BOTH_STRANDS_MODEL]);
            final double log10NotForeignBoth = Math.log10(1.0 - foreignAlleleLikelihoods[BOTH_STRANDS_MODEL]);

            // iterate over the reads again using the approximated beta constraint; LP_R_GF for each model
            double cumulativeForward = 0.0;
            double cumulativeReverse = 0.0;
            double cumulativeBoth = 0.0;
            for (int i = 0; i < readCount; i++) {
                final double log10LikelihoodReadForGenotype = log10LikelihoodsForGT[i];
                final double log10LikelihoodOfForeignAlleleGivenLPFCutoff = criticalValues[i] <= logProbFAllele ?
                        Double.NEGATIVE_INFINITY : log10LikelihoodsForF[i];
                cumulativeBoth += MathUtils.approximateLog10SumLog10(log10ForeignBoth + log10LikelihoodOfForeignAlleleGivenLPFCutoff, log10NotForeignBoth + log10LikelihoodReadForGenotype);
                if (isReverseStrand[i]) {
                    cumulativeReverse += MathUtils.approximateLog10SumLog10(log10ForeignReverse + log10LikelihoodOfForeignAlleleGivenLPFCutoff, log10NotForeignReverse + log10LikelihoodReadForGenotype);
                    cumulativeForward += log10LikelihoodReadForGenotype;
                } else {
                    cumulativeForward += MathUtils.approximateLog10SumLog10(log10ForeignForward + log10LikelihoodOfForeignAlleleGivenLPFCutoff, log10NotForeignForward + log10LikelihoodReadForGenotype);
                    cumulativeReverse += log10LikelihoodReadForGenotype;
                }
            }
            cumulativeLikelihoods[FORWARD_STRAND_MODEL] = cumulativeForward;
            cumulativeLikelihoods[REVERSE_STRAND_MODEL] = cumulativeReverse;
            cumulativeLikelihoods[BOTH_STRANDS_MODEL] = cumulativeBoth;

            for (int model = 0; model < STRAND_MODEL_COUNT; model++) {
                // Allele prior for error allele, plus posterior for foreign event, plus model posterior
                final double lpsi = logProbFAllele + cumulativeLikelihoods[model]; // NOTE unlike DRAGEN we apply the prior to the combined likelihoods array after the fact so gtAllelePrior is not included at this stage
                if (debugLines != null) {
                    debugLines.get(model).add("beta: " + foreignAlleleLikelihoods[model] + " localMaxLpspi: " + lpsi + " for lpf: " + logProbFAllele + " with LP_R_GF: " + cumulativeLikelihoods[model] + " index: " + thresholdIndex);
                }
                if (lpsi > maxLpspi[model]) {
                    maxLpspi[model] = lpsi;
                    lpfApplied[model] = logProbFAllele;
                }
            }
        }

        if (debugLines != null) {
            for (int model = 0; model < STRAND_MODEL_COUNT; model++) {
                HaplotypeCallerGenotypingDebugger.println(STRAND_MODEL_DEBUG_HEADERS[model]);
                debugLines.get(model).forEach(HaplotypeCallerGenotypingDebugger::println);
            }
        }

        // TODO soon should not need to use the LPF applied here...
        return new double[][]{
                {maxLpspi[FORWARD_STRAND_MODEL], lpfApplied[FORWARD_STRAND_MODEL]},
                {maxLpspi[REVERSE_STRAND_MODEL], lpfApplied[REVERSE_STRAND_MODEL]},
                {maxLpspi[BOTH_STRANDS_MODEL], lpfApplied[BOTH_STRANDS_MODEL]}};
    }

    /**
     * The reads FRD genotypes a sample with, held in primitive arrays in read-container order.
     *
     * Reads disqualified by the HMM have no likelihoods and are not genotyped, but their mapping qualities still
     * contribute critical thresholds.
     */
    private static final class FRDReads {
        /** Whether there are any read containers at all, disqualified or not. */
        final boolean anyReads;
        /** Phred-scaled mapping quality of every read container, in order. */
        final double[] allPhredScaledMappingQualities;
        /** Number of genotyped (not disqualified) reads; the per-read arrays below are indexed by genotyped read. */
        final int genotypedCount;
        /** Number of genotyped reads on the forward strand. */
        final int forwardCount;
        final boolean[] isReverseStrand;
        final double[] phredScaledMappingQualities;
        /** log10 likelihood of each genotyped read, by allele index and then read. */
        final double[][] log10LikelihoodsByAllele;
        /**
         * Each genotyped read's critical log10(P(F)) value for the current foreign allele, plus the tolerance below
         * which a threshold excludes the read's foreign-allele likelihood; set by {@link #updateCriticalValues}.
         * Reads disqualified by the HMM have no entry here, although their values still count as thresholds.
         */
        final double[] criticalValuesWithTolerance;

        <A extends Allele> FRDReads(final LikelihoodMatrix<GATKRead, A> sampleLikelihoods, final List<DRAGENGenotypesModel.DragenReadContainer> readContainers) {
            anyReads = !readContainers.isEmpty();
            allPhredScaledMappingQualities = new double[readContainers.size()];
            int genotyped = 0;
            for (int i = 0; i < readContainers.size(); i++) {
                final DRAGENGenotypesModel.DragenReadContainer container = readContainers.get(i);
                allPhredScaledMappingQualities[i] = container.getPhredScaledMappingQuality();
                if (!container.wasFilteredByHMM()) {
                    genotyped++;
                }
            }

            genotypedCount = genotyped;
            isReverseStrand = new boolean[genotypedCount];
            phredScaledMappingQualities = new double[genotypedCount];
            criticalValuesWithTolerance = new double[genotypedCount];
            log10LikelihoodsByAllele = new double[sampleLikelihoods.numberOfAlleles()][genotypedCount];
            int readIndex = 0;
            int forward = 0;
            for (int i = 0; i < readContainers.size(); i++) {
                final DRAGENGenotypesModel.DragenReadContainer container = readContainers.get(i);
                if (container.wasFilteredByHMM()) {
                    continue;
                }
                isReverseStrand[readIndex] = container.isReverseStrand();
                if (!isReverseStrand[readIndex]) {
                    forward++;
                }
                phredScaledMappingQualities[readIndex] = allPhredScaledMappingQualities[i];
                for (int allele = 0; allele < log10LikelihoodsByAllele.length; allele++) {
                    log10LikelihoodsByAllele[allele][readIndex] = sampleLikelihoods.get(allele, container.getIndexInLikelihoodsObject());
                }
                readIndex++;
            }
            forwardCount = forward;
        }

        /**
         * Sets each genotyped read's critical value for a foreign allele and returns the distinct critical thresholds
         * of all reads, disqualified ones included, in ascending order.
         *
         * @param log10MapqPriorAdjustment log10 prior adjustment for the foreign allele
         */
        double[] updateCriticalValues(final double log10MapqPriorAdjustment) {
            for (int i = 0; i < genotypedCount; i++) {
                criticalValuesWithTolerance[i] = (phredScaledMappingQualities[i] * -0.1 + log10MapqPriorAdjustment) + 0.0000001;
            }
            final double[] thresholds = new double[allPhredScaledMappingQualities.length];
            for (int i = 0; i < thresholds.length; i++) {
                thresholds[i] = allPhredScaledMappingQualities[i] * -0.1 + log10MapqPriorAdjustment;
            }
            // Arrays.sort orders doubles as Double.compareTo does, and dropping each value that Double.compare finds
            // equal to its predecessor keeps one of each value Double.equals distinguishes, -0.0 and NaN included.
            Arrays.sort(thresholds);
            int distinct = 0;
            for (final double threshold : thresholds) {
                if (distinct == 0 || Double.compare(thresholds[distinct - 1], threshold) != 0) {
                    thresholds[distinct++] = threshold;
                }
            }
            return Arrays.copyOf(thresholds, distinct);
        }
    }
}
