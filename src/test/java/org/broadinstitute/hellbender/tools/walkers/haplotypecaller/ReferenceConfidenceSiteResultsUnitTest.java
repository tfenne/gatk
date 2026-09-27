package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Tests of {@link ReferenceConfidenceModel#calculateSiteResults}, which accumulates reference confidence evidence read
 * by read. The expected values come from per-position pileups fed through the model's per-pileup likelihood method.
 */
public final class ReferenceConfidenceSiteResultsUnitTest extends GATKBaseTest {
    private static final String CONTIG = "1";
    private static final int CONTIG_LENGTH = 1000;
    private static final SimpleInterval SPAN = new SimpleInterval(CONTIG, 201, 400);
    private static final int PADDING = 100;
    private static final int PLOIDY = 2;
    private static final int MAX_INDEL_SIZE = 10;
    private static final byte[] QUALS_TO_DRAW = {2, 5, 6, 7, 15, 30, 40};

    private final String sample = "NA12878";
    private final SampleList samples = SampleList.singletonSampleList(sample);
    private SAMFileHeader header;
    private SAMReadGroupRecord readGroup;
    private byte[] contigBases;
    private ReferenceConfidenceModel model;
    private ReferenceConfidenceModel flowModel;
    private ReferenceConfidenceModel noSoftClipModel;

    @BeforeClass
    public void setUp() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, CONTIG_LENGTH);
        readGroup = new SAMReadGroupRecord("RG1");
        readGroup.setSample(sample);
        header.addReadGroup(readGroup);
        final Random rng = new Random(17);
        contigBases = new byte[CONTIG_LENGTH];
        for (int i = 0; i < CONTIG_LENGTH; i++) {
            contigBases[i] = randomBase(rng);
        }
        model = new ReferenceConfidenceModel(samples, header, MAX_INDEL_SIZE, -1, (byte) 30, true, false);
        flowModel = new ReferenceConfidenceModel(samples, header, MAX_INDEL_SIZE, -1, (byte) 30, true, true);
        noSoftClipModel = new ReferenceConfidenceModel(samples, header, MAX_INDEL_SIZE, -1, (byte) 30, false, false);
    }

    private List<GATKRead> randomReads(final int seed, final int count) {
        final Random rng = new Random(seed);
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            reads.add(randomRead(rng, "read" + i));
        }
        return reads;
    }

    @Test
    public void siteResultsMatchPileupsForRandomReads() {
        for (int seed = 1; seed <= 20; seed++) {
            assertSiteResultsMatchPileups(model, randomReads(seed, 80), PLOIDY, "seed " + seed);
        }
    }

    @Test
    public void siteResultsMatchPileupsWhenPositionsHaveMoreIndelInformativeReadsThanAreCounted() {
        for (int seed = 1; seed <= 5; seed++) {
            final int maxInformative = assertSiteResultsMatchPileups(model, randomReads(seed, 600), PLOIDY, "deep seed " + seed);
            Assert.assertTrue(maxInformative > 2 * ReferenceConfidenceModel.MAX_N_INDEL_INFORMATIVE_READS,
                    "the pileups should be deep enough to fill the count: " + maxInformative);
        }
    }

    @Test
    public void siteResultsMatchPileupsForRandomReadsInFlowMode() {
        for (int seed = 1; seed <= 10; seed++) {
            assertSiteResultsMatchPileups(flowModel, randomReads(seed, 80), PLOIDY, "flow seed " + seed);
        }
    }

    @Test
    public void siteResultsMatchPileupsWhenSoftClippedBasesAreNotUsed() {
        for (int seed = 1; seed <= 10; seed++) {
            final Random rng = new Random(seed);
            final List<GATKRead> reads = randomReads(seed, 80);
            for (final GATKRead read : reads) {
                read.setAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG, read.getStart() + rng.nextInt(8));
                read.setAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_END_TAG, read.getEnd() - rng.nextInt(8));
            }
            assertSiteResultsMatchPileups(noSoftClipModel, reads, PLOIDY, "no soft clips seed " + seed);
        }
    }

    @Test
    public void siteResultsMatchPileupsAtPloidyOneAndThree() {
        for (final int ploidy : Arrays.asList(1, 3)) {
            for (int seed = 1; seed <= 5; seed++) {
                assertSiteResultsMatchPileups(model, randomReads(seed, 80), ploidy, "ploidy " + ploidy + " seed " + seed);
            }
        }
    }

    @Test
    public void readsAreAccumulatedInCoordinateOrderRegardlessOfInputOrder() {
        final List<GATKRead> reads = randomReads(99, 60);
        final List<ReferenceConfidenceResult> fromSorted = siteResults(model, reads);
        final List<GATKRead> shuffled = new ArrayList<>(reads);
        Collections.shuffle(shuffled, new Random(99));
        assertSameSites(siteResults(model, shuffled), fromSorted, "shuffled input");
    }

    @Test
    public void positionsWithoutReadsHaveNoDepthAndFlatLikelihoods() {
        for (final ReferenceConfidenceResult site : siteResults(model, Collections.emptyList())) {
            final RefVsAnyResult result = (RefVsAnyResult) site;
            Assert.assertEquals(result.getDP(), 0);
            Assert.assertEquals(result.getAD(), new int[]{0, 0});
            for (final double likelihood : result.genotypeLikelihoods) {
                Assert.assertEquals(likelihood, 0.0);
            }
        }
    }

    @Test
    public void readsOutsideTheSpanContributeNothing() {
        final List<GATKRead> reads = Arrays.asList(
                matchingRead("left", SPAN.getStart() - 50, 50),
                matchingRead("right", SPAN.getEnd() + 1, 50));
        assertSameSites(siteResults(model, reads), siteResults(model, Collections.emptyList()), "reads outside the span");
    }

    @Test
    public void everyCoveredPositionCountsTheRead() {
        final int start = SPAN.getStart() + 10;
        final int length = 50;
        final List<ReferenceConfidenceResult> sites = siteResults(model, Collections.singletonList(matchingRead("r", start, length)));
        for (int position = SPAN.getStart(); position <= SPAN.getEnd(); position++) {
            final boolean covered = position >= start && position < start + length;
            Assert.assertEquals(sites.get(position - SPAN.getStart()).getDP(), covered ? 1 : 0, "position " + position);
        }
    }

    @Test
    public void referenceSkipsAreNotEvidence() {
        final int start = SPAN.getStart() + 10;
        final byte[] bases = Arrays.copyOfRange(contigBases, start - 1, start - 1 + 40);
        System.arraycopy(contigBases, start - 1 + 60, bases, 20, 20);
        final GATKRead read = read("skip", start, bases, quals(bases.length, (byte) 30), "20M40N20M");
        final List<ReferenceConfidenceResult> sites = siteResults(model, Collections.singletonList(read));
        for (int position = SPAN.getStart(); position <= SPAN.getEnd(); position++) {
            final boolean aligned = (position >= start && position < start + 20) || (position >= start + 60 && position < start + 80);
            Assert.assertEquals(sites.get(position - SPAN.getStart()).getDP(), aligned ? 1 : 0, "position " + position);
        }
    }

    @Test
    public void deletedPositionsCountAsNonReferenceEvidence() {
        final int start = SPAN.getStart() + 10;
        final byte[] bases = new byte[20];
        System.arraycopy(contigBases, start - 1, bases, 0, 10);
        System.arraycopy(contigBases, start - 1 + 13, bases, 10, 10);
        final GATKRead read = read("del", start, bases, quals(20, (byte) 30), "10M3D10M");
        final List<ReferenceConfidenceResult> sites = siteResults(model, Collections.singletonList(read));
        for (int position = start; position < start + 23; position++) {
            final ReferenceConfidenceResult site = sites.get(position - SPAN.getStart());
            final boolean deleted = position >= start + 10 && position < start + 13;
            Assert.assertEquals(site.nonRefDepth, deleted ? 1 : 0, "position " + position);
            Assert.assertEquals(site.refDepth, deleted ? 0 : 1, "position " + position);
        }
    }

    @Test
    public void basesInsideTheAdaptorAreNotEvidence() {
        final int start = SPAN.getStart() + 10;
        final int fragmentLength = 60;
        final GATKRead read = matchingRead("short-fragment", start, 100);
        read.setIsPaired(true);
        read.setMatePosition(CONTIG, start + 20);
        read.setMateIsReverseStrand(true);
        read.setFragmentLength(fragmentLength);
        final List<ReferenceConfidenceResult> sites = siteResults(model, Collections.singletonList(read));
        for (int position = start; position < start + 100; position++) {
            final boolean insideAdaptor = position >= start + fragmentLength;
            Assert.assertEquals(sites.get(position - SPAN.getStart()).getDP(), insideAdaptor ? 0 : 1, "position " + position);
        }
    }

    @Test
    public void originallySoftClippedBasesAreSkippedWhenSoftClipsAreNotUsed() {
        final ReferenceConfidenceModel noSoftClipModel = new ReferenceConfidenceModel(samples, header, MAX_INDEL_SIZE, -1, (byte) 30, false, false);
        final int start = SPAN.getStart() + 10;
        final GATKRead read = matchingRead("reverted", start, 30);
        read.setAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG, start + 5);
        read.setAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_END_TAG, start + 20);
        final List<ReferenceConfidenceResult> sites = siteResults(noSoftClipModel, Collections.singletonList(read));
        for (int position = start; position < start + 30; position++) {
            final boolean originallyAligned = position >= start + 5 && position <= start + 20;
            Assert.assertEquals(sites.get(position - SPAN.getStart()).getDP(), originallyAligned ? 1 : 0, "position " + position);
        }
    }

    // ---------------------------------------------------------------------------------------------------------------
    // expected values from the per-position pileup path
    // ---------------------------------------------------------------------------------------------------------------

    /**
     * Asserts that the read-major site results equal those computed from a pileup at each position.
     *
     * @return the largest number of indel-informative reads at any position, before the model's cap
     */
    private int assertSiteResultsMatchPileups(final ReferenceConfidenceModel model, final List<GATKRead> reads, final int ploidy, final String label) {
        final AlleleLikelihoods<GATKRead, Haplotype> likelihoods = likelihoods(reads);
        final AssemblyRegion region = region();
        final byte[] paddedRef = paddedRef();
        final List<ReferenceConfidenceResult> actual = model.calculateSiteResults(likelihoods, region, paddedRef, ploidy);

        final List<ReadPileup> pileups = AssemblyBasedCallerUtils.getPileupsOverReference(header, SPAN, likelihoods, samples);
        Assert.assertEquals(actual.size(), pileups.size(), label);
        final int globalRefOffset = SPAN.getStart() - region.getPaddedSpan().getStart();
        final List<ReferenceConfidenceResult> expected = new ArrayList<>();
        int maxInformative = 0;
        for (int i = 0; i < pileups.size(); i++) {
            final ReadPileup pileup = pileups.get(i);
            final int refOffset = i + globalRefOffset;
            final RefVsAnyResult site = (RefVsAnyResult) model.calcGenotypeLikelihoodsOfRefVsAny(ploidy, pileup, paddedRef[refOffset], ReferenceConfidenceModel.BASE_QUAL_THRESHOLD, null, true);
            int informative = 0;
            for (final PileupElement element : pileup) {
                if (element.isBeforeDeletionStart() || element.isBeforeInsertion() || element.isDeletion()) {
                    continue;
                }
                final int alignedOffset = ReferenceConfidenceTestUtils.referenceAlignedOffset(element);
                final BitSet bits = ReferenceConfidenceModel.indelInformativeBases(element.getRead(), alignedOffset, paddedRef, refOffset, MAX_INDEL_SIZE);
                if (bits.get(alignedOffset)) {
                    informative++;
                }
            }
            model.applyIndelRefConfidence(ploidy, informative, site);
            expected.add(site);
            maxInformative = Math.max(maxInformative, informative);
        }
        assertSameSites(actual, expected, label);
        return maxInformative;
    }

    private static void assertSameSites(final List<ReferenceConfidenceResult> actual, final List<ReferenceConfidenceResult> expected, final String label) {
        Assert.assertEquals(actual.size(), expected.size(), label);
        for (int i = 0; i < actual.size(); i++) {
            final RefVsAnyResult a = (RefVsAnyResult) actual.get(i);
            final RefVsAnyResult e = (RefVsAnyResult) expected.get(i);
            final String site = label + ", site " + i;
            Assert.assertEquals(a.refDepth, e.refDepth, site + " refDepth");
            Assert.assertEquals(a.nonRefDepth, e.nonRefDepth, site + " nonRefDepth");
            Assert.assertEquals(a.genotypeLikelihoods.length, e.genotypeLikelihoods.length, site);
            for (int k = 0; k < a.genotypeLikelihoods.length; k++) {
                Assert.assertTrue(Double.doubleToLongBits(a.genotypeLikelihoods[k]) == Double.doubleToLongBits(e.genotypeLikelihoods[k]),
                        site + " likelihood " + k + ": " + a.genotypeLikelihoods[k] + " vs " + e.genotypeLikelihoods[k]);
            }
            Assert.assertEquals(a.finalPhredScaledGenotypeLikelihoods, e.finalPhredScaledGenotypeLikelihoods, site + " PLs");
        }
    }

    // ---------------------------------------------------------------------------------------------------------------
    // fixtures
    // ---------------------------------------------------------------------------------------------------------------

    private List<ReferenceConfidenceResult> siteResults(final ReferenceConfidenceModel model, final List<GATKRead> reads) {
        return model.calculateSiteResults(likelihoods(reads), region(), paddedRef(), PLOIDY);
    }

    private AssemblyRegion region() {
        return new AssemblyRegion(SPAN, PADDING, header);
    }

    private byte[] paddedRef() {
        final SimpleInterval padded = region().getPaddedSpan();
        return Arrays.copyOfRange(contigBases, padded.getStart() - 1, padded.getEnd());
    }

    private AlleleLikelihoods<GATKRead, Haplotype> likelihoods(final List<GATKRead> reads) {
        final AssemblyRegion region = region();
        final Haplotype refHaplotype = ReferenceConfidenceModel.createReferenceHaplotype(region, paddedRef(), region.getPaddedSpan());
        return new AlleleLikelihoods<>(samples, new IndexedAlleleList<>(refHaplotype), Map.of(sample, new ArrayList<>(reads)));
    }

    private GATKRead read(final String name, final int start, final byte[] bases, final byte[] quals, final String cigar) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, start, bases, quals, cigar);
        read.setReadGroup(readGroup.getId());
        return read;
    }

    private GATKRead matchingRead(final String name, final int start, final int length) {
        return read(name, start, Arrays.copyOfRange(contigBases, start - 1, start - 1 + length), quals(length, (byte) 30), length + "M");
    }

    private static byte[] quals(final int length, final byte qual) {
        final byte[] quals = new byte[length];
        Arrays.fill(quals, qual);
        return quals;
    }

    private static byte randomBase(final Random rng) {
        return (byte) "ACGT".charAt(rng.nextInt(4));
    }

    /**
     * A read with a random mix of aligned blocks (as M or as runs of = and X), insertions, deletions, reference skips,
     * soft and hard clips, random base qualities on both sides of the quality threshold, and sometimes a short fragment
     * on either strand that puts one end of the read inside the adaptor.
     */
    private GATKRead randomRead(final Random rng, final String name) {
        final int start = SPAN.getStart() - 80 + rng.nextInt(SPAN.size() + 100);
        final List<CigarElement> elements = new ArrayList<>();
        final List<Byte> bases = new ArrayList<>();
        int refIndex = start - 1;

        if (rng.nextInt(6) == 0) {
            addToCigar(elements, 1 + rng.nextInt(5), CigarOperator.H);
        }
        if (rng.nextInt(4) == 0) {
            addToCigar(elements, addRandomBases(bases, rng, 1 + rng.nextInt(4)), CigarOperator.S);
        }
        final int blocks = 1 + rng.nextInt(4);
        for (int b = 0; b <= blocks; b++) {
            final int blockLength = b < blocks ? 5 + rng.nextInt(36) : 5 + rng.nextInt(26);
            final boolean explicitMatches = rng.nextInt(4) == 0;
            for (int k = 0; k < blockLength; k++) {
                final byte refBase = contigBases[refIndex++];
                final boolean mismatch = rng.nextInt(20) == 0;
                bases.add(mismatch ? differentBase(rng, refBase) : refBase);
                addToCigar(elements, 1, explicitMatches ? (mismatch ? CigarOperator.X : CigarOperator.EQ) : CigarOperator.M);
            }
            if (b < blocks) {
                switch (rng.nextInt(5)) {
                    case 1:
                        addToCigar(elements, addRandomBases(bases, rng, 1 + rng.nextInt(3)), CigarOperator.I);
                        break;
                    case 2:
                        final int deleted = 1 + rng.nextInt(3);
                        refIndex += deleted;
                        addToCigar(elements, deleted, CigarOperator.D);
                        break;
                    case 3:
                        final int skipped = 10 + rng.nextInt(30);
                        refIndex += skipped;
                        addToCigar(elements, skipped, CigarOperator.N);
                        break;
                    default:
                        break;
                }
            }
        }
        if (rng.nextInt(4) == 0) {
            addToCigar(elements, addRandomBases(bases, rng, 1 + rng.nextInt(4)), CigarOperator.S);
        }
        if (rng.nextInt(6) == 0) {
            addToCigar(elements, 1 + rng.nextInt(5), CigarOperator.H);
        }

        final byte[] readBases = new byte[bases.size()];
        final byte[] readQuals = new byte[bases.size()];
        for (int k = 0; k < readBases.length; k++) {
            readBases[k] = bases.get(k);
            readQuals[k] = QUALS_TO_DRAW[rng.nextInt(QUALS_TO_DRAW.length)];
        }
        final GATKRead read = read(name, start, readBases, readQuals, TextCigarCodec.encode(new Cigar(elements)));
        final int fragment = rng.nextInt(5);
        if (fragment == 0) {
            read.setIsPaired(true);
            read.setMatePosition(CONTIG, start + 10);
            read.setMateIsReverseStrand(true);
            read.setFragmentLength(20 + rng.nextInt(80));
        } else if (fragment == 1) {
            read.setIsPaired(true);
            read.setIsReverseStrand(true);
            read.setMatePosition(CONTIG, start + 5 + rng.nextInt(20));
            read.setMateIsReverseStrand(false);
            read.setFragmentLength(-(20 + rng.nextInt(80)));
        }
        return read;
    }

    /** Appends an element, merging it into the previous one when the operators match. */
    private static void addToCigar(final List<CigarElement> elements, final int length, final CigarOperator op) {
        final int last = elements.size() - 1;
        if (last >= 0 && elements.get(last).getOperator() == op) {
            elements.set(last, new CigarElement(elements.get(last).getLength() + length, op));
        } else {
            elements.add(new CigarElement(length, op));
        }
    }

    private static int addRandomBases(final List<Byte> bases, final Random rng, final int count) {
        for (int k = 0; k < count; k++) {
            bases.add(randomBase(rng));
        }
        return count;
    }

    private static byte differentBase(final Random rng, final byte refBase) {
        byte base = randomBase(rng);
        while (base == refBase) {
            base = randomBase(rng);
        }
        return base;
    }
}
