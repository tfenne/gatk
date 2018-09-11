package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of alt reads with an OA tag that doesn't match the current alignment contig.")
public class OriginalAlignment extends GenotypeAnnotation implements Annotation {

    public static final String OA_NOT_CURRENT_CONTIG = "OA_NOT_CURRENT_CONTIG";
    private static final int CONTIG_INDEX_OF_OA_TAG = 0;
    private static final String OA_TAG_SEPARATOR = ",";

    @Override
    public void annotate(ReferenceContext ref, VariantContext vc, Genotype g, GenotypeBuilder gb, ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        Utils.nonNull(likelihoods);

        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY, () -> null, -1);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);
        final Allele altAlelle = vc.getAlternateAllele(indexOfMaxTumorLod);
        final Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAllelesBreakingTies(g.getSampleName());
        final String currentContig = ref.getInterval().getContig();

        int nonChrMAlt = (int) bestAlleles.stream().filter(ba -> ba.read.hasAttribute("OA") && ba.isInformative() &&
                ba.allele.equals(altAlelle) && !ba.read.getAttributeAsString("OA").split(OA_TAG_SEPARATOR)[CONTIG_INDEX_OF_OA_TAG].equals(currentContig)).count();
        gb.attribute(OA_NOT_CURRENT_CONTIG, nonChrMAlt);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(OA_NOT_CURRENT_CONTIG, 1, VCFHeaderLineType.Integer, "Number of alt reads whose original alignment doesn't match the current contig."));
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(OA_NOT_CURRENT_CONTIG);
    }
}
