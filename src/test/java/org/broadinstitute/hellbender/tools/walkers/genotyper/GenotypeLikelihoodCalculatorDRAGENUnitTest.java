package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public final class GenotypeLikelihoodCalculatorDRAGENUnitTest {
    private static final SAMFileHeader HEADER = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1_000_000);
    private static final List<Allele> SNP_ALLELES = List.of(Allele.create("A", true), Allele.create("C"));
    private static final List<Allele> SNP_AND_INDEL_ALLELES = List.of(Allele.create("A", true), Allele.create("C"), Allele.create("AT"));
    private static final int DIPLOID = 2;
    private static final int HOM_REF = 0;
    private static final int HET = 1;
    private static final int HOM_ALT = 2;
    private static final double SNP_HET_PRIOR = 34.77;
    private static final double INDEL_HET_PRIOR = 41.0;

    /** A read overlapping the site: its mapping quality, strand, and log10 likelihood of each allele in order. */
    private record SiteRead(int mappingQuality, boolean reverseStrand, double... log10Likelihoods) { }

    /** The inputs FRD genotypes one sample with. */
    private record Site(int ploidy, LikelihoodMatrix<GATKRead, Allele> likelihoods, List<DRAGENGenotypesModel.DragenReadContainer> containers,
                        double[] standardGenotypeLikelihoods) {
        double[] frd() {
            return frd(0);
        }

        double[] frd(final int maxEffectiveDepthForHetAdjustment) {
            return GenotypeLikelihoodCalculatorDRAGEN.calculateFRDLikelihoods(ploidy, likelihoods, standardGenotypeLikelihoods,
                    containers, SNP_HET_PRIOR, INDEL_HET_PRIOR, maxEffectiveDepthForHetAdjustment);
        }
    }

    private static Site site(final List<SiteRead> genotypedReads, final List<SiteRead> hmmDisqualifiedReads) {
        return site(SNP_ALLELES, DIPLOID, genotypedReads, hmmDisqualifiedReads);
    }

    private static Site site(final List<Allele> alleles, final int ploidy, final List<SiteRead> genotypedReads, final List<SiteRead> hmmDisqualifiedReads) {
        final List<GATKRead> reads = new ArrayList<>();
        final List<DRAGENGenotypesModel.DragenReadContainer> containers = new ArrayList<>();
        for (final SiteRead siteRead : genotypedReads) {
            containers.add(new DRAGENGenotypesModel.DragenReadContainer(read(siteRead, reads.size()), 5, 1000, reads.size()));
            reads.add(containers.get(containers.size() - 1).underlyingRead);
        }
        for (final SiteRead siteRead : hmmDisqualifiedReads) {
            containers.add(new DRAGENGenotypesModel.DragenReadContainer(read(siteRead, containers.size()), 5, 1000, -1));
        }
        final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods = new AlleleLikelihoods<>(new IndexedSampleList("sample"),
                new IndexedAlleleList<>(alleles), Map.of("sample", reads));
        final LikelihoodMatrix<GATKRead, Allele> matrix = alleleLikelihoods.sampleMatrix(0);
        for (int r = 0; r < genotypedReads.size(); r++) {
            for (int a = 0; a < alleles.size(); a++) {
                matrix.set(a, r, genotypedReads.get(r).log10Likelihoods()[a]);
            }
        }
        return new Site(ploidy, matrix, containers, GenotypeLikelihoodCalculator.computeLog10GenotypeLikelihoods(ploidy, matrix));
    }

    private static GATKRead read(final SiteRead siteRead, final int index) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(HEADER, "read" + index, 0, 1000, 20);
        read.setMappingQuality(siteRead.mappingQuality());
        read.setIsReverseStrand(siteRead.reverseStrand());
        return read;
    }

    /** Twenty well-mapped reference reads split across strands, plus the given alt reads. */
    private static List<SiteRead> refReadsPlus(final List<SiteRead> altReads) {
        final List<SiteRead> reads = new ArrayList<>();
        for (int i = 0; i < 20; i++) {
            reads.add(new SiteRead(60, i % 2 == 1, 0.0, -5.0));
        }
        reads.addAll(altReads);
        return reads;
    }

    private static List<SiteRead> altReads(final int count, final int mappingQuality) {
        final List<SiteRead> reads = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            reads.add(new SiteRead(mappingQuality, i % 2 == 1, -5.0, 0.0));
        }
        return reads;
    }

    /**
     * Twelve reference reads and four alt reads on both strands, with seven distinct mapping qualities and likelihoods
     * that are not round numbers, for comparing exact values.
     */
    private static List<SiteRead> mixedReads() {
        final List<SiteRead> reads = new ArrayList<>();
        final int[] refMappingQualities = {60, 45, 20};
        for (int i = 0; i < 12; i++) {
            reads.add(new SiteRead(refMappingQualities[i % 3], i % 2 == 1, -0.02 - 0.001 * i, -4.3 + 0.07 * i));
        }
        reads.add(new SiteRead(3, false, -3.7, -0.05));
        reads.add(new SiteRead(12, true, -2.9, -0.11));
        reads.add(new SiteRead(25, false, -4.1, -0.02));
        reads.add(new SiteRead(40, false, -3.3, -0.3));
        return reads;
    }

    /** Asserts that every value is bit-for-bit the expected one, so any change of summation order is caught. */
    private static void assertBitIdentical(final double[] actual, final double[] expected) {
        Assert.assertEquals(actual.length, expected.length);
        for (int i = 0; i < actual.length; i++) {
            Assert.assertEquals(Double.doubleToLongBits(actual[i]), Double.doubleToLongBits(expected[i]),
                    "index " + i + ": " + Arrays.toString(actual) + " vs expected " + Arrays.toString(expected));
        }
    }

    @Test
    public void noReadsLeavesEveryGenotypeUnscored() {
        final double[] frd = site(List.of(), List.of()).frd();
        Assert.assertEquals(frd, new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY});
    }

    @Test
    public void onlyHomozygousGenotypesAreScored() {
        final double[] frd = site(refReadsPlus(altReads(3, 5)), List.of()).frd();
        Assert.assertEquals(frd[HET], Double.NEGATIVE_INFINITY);
        Assert.assertTrue(Double.isFinite(frd[HOM_REF]), Arrays.toString(frd));
        Assert.assertTrue(Double.isFinite(frd[HOM_ALT]), Arrays.toString(frd));
    }

    @Test
    public void poorlyMappedAltReadsAreExplainedAsForeignReadsUnderHomRef() {
        final Site site = site(refReadsPlus(altReads(3, 5)), List.of());
        Assert.assertTrue(site.frd()[HOM_REF] > site.standardGenotypeLikelihoods()[HOM_REF],
                "FRD hom-ref " + site.frd()[HOM_REF] + " should exceed the standard hom-ref likelihood " + site.standardGenotypeLikelihoods()[HOM_REF]);
    }

    @Test
    public void poorlyMappedAltReadsSupportHomRefMoreThanWellMappedOnes() {
        final double poorlyMapped = site(refReadsPlus(altReads(3, 5)), List.of()).frd()[HOM_REF];
        final double wellMapped = site(refReadsPlus(altReads(3, 60)), List.of()).frd()[HOM_REF];
        Assert.assertTrue(poorlyMapped > wellMapped, "poorly mapped " + poorlyMapped + " vs well mapped " + wellMapped);
    }

    @Test
    public void altReadsOnOneStrandAreExplainedAtLeastAsWellAsAltReadsOnBothStrands() {
        final List<SiteRead> forwardOnly = new ArrayList<>();
        for (int i = 0; i < 3; i++) {
            forwardOnly.add(new SiteRead(5, false, -5.0, 0.0));
        }
        final double oneStrand = site(refReadsPlus(forwardOnly), List.of()).frd()[HOM_REF];
        final double bothStrands = site(refReadsPlus(altReads(3, 5)), List.of()).frd()[HOM_REF];
        Assert.assertTrue(oneStrand >= bothStrands, "one strand " + oneStrand + " vs both strands " + bothStrands);
    }

    @Test
    public void hmmDisqualifiedReadsAreNotGenotyped() {
        final List<SiteRead> reads = refReadsPlus(altReads(3, 5));
        final double[] without = site(reads, List.of()).frd();
        // The disqualified reads' mapping qualities match genotyped reads, so they add no new critical threshold.
        final double[] with = site(reads, List.of(new SiteRead(60, false, -50.0, 0.0), new SiteRead(5, true, 0.0, -50.0))).frd();
        Assert.assertEquals(with, without);
    }

    @Test
    public void hmmDisqualifiedReadsContributeCriticalThresholds() {
        // With no alt support the highest critical threshold scores best, so a disqualified read mapped worse than
        // every genotyped read raises the hom-ref score.
        final List<SiteRead> reads = refReadsPlus(List.of());
        final double without = site(reads, List.of()).frd()[HOM_REF];
        final double with = site(reads, List.of(new SiteRead(1, false, 0.0, 0.0))).frd()[HOM_REF];
        Assert.assertTrue(with > without, "with disqualified read " + with + " vs without " + without);
    }

    // The expected values below are the outputs of the FRD implementation this one replaced, on the same inputs.

    @Test
    public void frdMatchesReferenceValuesForReadsOnBothStrandsWithManyMappingQualities() {
        assertBitIdentical(site(mixedReads(), List.of()).frd(), EXPECTED_MIXED);
    }

    @Test
    public void frdMatchesReferenceValuesWhenEveryReadIsOnOneStrand() {
        final List<SiteRead> forwardReads = new ArrayList<>();
        for (final SiteRead read : mixedReads()) {
            forwardReads.add(new SiteRead(read.mappingQuality(), false, read.log10Likelihoods()));
        }
        assertBitIdentical(site(forwardReads, List.of()).frd(), EXPECTED_ONE_STRAND);
    }

    @Test
    public void frdMatchesReferenceValuesWhenEveryReadIsDisqualified() {
        assertBitIdentical(site(List.of(), mixedReads()).frd(), EXPECTED_ALL_DISQUALIFIED);
    }

    @Test
    public void frdMatchesReferenceValuesWithNegativeInfiniteLikelihoods() {
        final List<SiteRead> reads = mixedReads();
        reads.add(new SiteRead(8, true, 0.0, Double.NEGATIVE_INFINITY));
        reads.add(new SiteRead(15, false, Double.NEGATIVE_INFINITY, -0.4));
        assertBitIdentical(site(reads, List.of(new SiteRead(33, true, 0.0, 0.0))).frd(), EXPECTED_NEGATIVE_INFINITY);
    }

    @Test
    public void frdMatchesReferenceValuesForASnpAndAnIndelAllele() {
        final List<SiteRead> reads = new ArrayList<>();
        for (final SiteRead read : mixedReads()) {
            reads.add(new SiteRead(read.mappingQuality(), read.reverseStrand(), read.log10Likelihoods()[0], read.log10Likelihoods()[1], -3.9));
        }
        reads.add(new SiteRead(7, true, -4.4, -3.8, -0.03));
        reads.add(new SiteRead(18, false, -3.6, -4.2, -0.2));
        assertBitIdentical(site(SNP_AND_INDEL_ALLELES, DIPLOID, reads, List.of()).frd(), EXPECTED_SNP_AND_INDEL);
    }

    @Test
    public void frdMatchesReferenceValuesForAHaploidSample() {
        assertBitIdentical(site(SNP_ALLELES, 1, mixedReads(), List.of()).frd(), EXPECTED_HAPLOID);
    }

    @Test
    public void frdMatchesReferenceValuesForATriploidSample() {
        assertBitIdentical(site(SNP_ALLELES, 3, mixedReads(), List.of()).frd(), EXPECTED_TRIPLOID);
    }

    @Test
    public void frdMatchesReferenceValuesWithTheMaxEffectiveDepthAdjustment() {
        assertBitIdentical(site(mixedReads(), List.of()).frd(5), EXPECTED_MAX_EFFECTIVE_DEPTH);
    }

    private static final double[] EXPECTED_MIXED = {-11.499925625747858, Double.NEGATIVE_INFINITY, -9.60042445240666};
    private static final double[] EXPECTED_ONE_STRAND = {-11.499925625747858, Double.NEGATIVE_INFINITY, -9.60042445240666};
    private static final double[] EXPECTED_ALL_DISQUALIFIED = {-3.777, Double.NEGATIVE_INFINITY, -0.30000000000000004};
    private static final double[] EXPECTED_NEGATIVE_INFINITY = {-12.611704326888452, Double.NEGATIVE_INFINITY, -10.602484443734621};
    private static final double[] EXPECTED_SNP_AND_INDEL = {
            -19.353285504705806, Double.NEGATIVE_INFINITY, -17.40783885631723,
            Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, -23.190052105130555};
    private static final double[] EXPECTED_HAPLOID = {-11.499925625747858, -9.60042445240666};
    private static final double[] EXPECTED_TRIPLOID = {-11.499925625747858, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, -9.60042445240666};
    private static final double[] EXPECTED_MAX_EFFECTIVE_DEPTH = {-12.333080857761438, Double.NEGATIVE_INFINITY, -9.807382574425647};
}
