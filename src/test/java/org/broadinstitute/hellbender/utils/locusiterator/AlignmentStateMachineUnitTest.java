package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public final class AlignmentStateMachineUnitTest extends LocusIteratorByStateBaseTest {
    @DataProvider(name = "AlignmentStateMachineTest")
    public Object[][] makeAlignmentStateMachineTest() {
        return createLIBSTests(
                Arrays.asList(1, 2),
                Arrays.asList(1, 2, 3, 4));
    }

    @Test(dataProvider = "AlignmentStateMachineTest")
    public void testAlignmentStateMachineTest(LIBSTest params) {
        final GATKRead read = params.makeRead();
        final AlignmentStateMachine state = new AlignmentStateMachine(read);
        final LIBS_position tester = new LIBS_position(read);

        // min is one because always visit something, even for 10I reads
        final int expectedBpToVisit = read.getEnd() - read.getStart() + 1;

        Assert.assertSame(state.getRead(), read);
        Assert.assertNotNull(state.toString());

        int bpVisited = 0;
        int lastOffset = -1;

        // TODO -- more tests about test state machine state before first step?
        Assert.assertTrue(state.isLeftEdge());
        Assert.assertNull(state.getCigarOperator());
        Assert.assertNotNull(state.toString());
        Assert.assertEquals(state.getReadOffset(), -1);
        Assert.assertEquals(state.getGenomeOffset(), -1);
        Assert.assertEquals(state.getCurrentCigarElementOffset(), -1);
        Assert.assertEquals(state.getCurrentCigarElement(), null);

        while ( state.stepForwardOnGenome() != null ) {
            Assert.assertNotNull(state.toString());

            tester.stepForwardOnGenome();

            Assert.assertTrue(state.getReadOffset() >= lastOffset, "Somehow read offsets are decreasing: lastOffset " + lastOffset + " current " + state.getReadOffset());
            Assert.assertEquals(state.getReadOffset(), tester.getCurrentReadOffset(), "Read offsets are wrong at " + bpVisited);

            Assert.assertFalse(state.isLeftEdge());

            Assert.assertEquals(state.getCurrentCigarElement(), read.getCigar().getCigarElement(tester.currentOperatorIndex), "CigarElement index failure");
            Assert.assertEquals(state.getOffsetIntoCurrentCigarElement(), tester.getCurrentPositionOnOperatorBase0(), "CigarElement index failure");

            Assert.assertEquals(read.getCigar().getCigarElement(state.getCurrentCigarElementOffset()), state.getCurrentCigarElement(), "Current cigar element isn't what we'd get from the read itself");

            Assert.assertTrue(state.getOffsetIntoCurrentCigarElement() >= 0, "Offset into current cigar too small");
            Assert.assertTrue(state.getOffsetIntoCurrentCigarElement() < state.getCurrentCigarElement().getLength(), "Offset into current cigar too big");

            Assert.assertEquals(state.getGenomeOffset(), tester.getCurrentGenomeOffsetBase0(), "Offset from alignment start is bad");
            Assert.assertEquals(state.getGenomePosition(), tester.getCurrentGenomeOffsetBase0() + read.getStart(), "GenomePosition start is bad");
            Assert.assertEquals(state.getLocation().size(), 1, "GenomeLoc position should have size == 1");
            Assert.assertEquals(state.getLocation().getStart(), state.getGenomePosition(), "GenomeLoc position is bad");
            // most tests of this functionality are in LIBS
            Assert.assertNotNull(state.makePileupElement());

            lastOffset = state.getReadOffset();
            bpVisited++;
        }

        Assert.assertEquals(bpVisited, expectedBpToVisit, "Didn't visit the expected number of bp");
        Assert.assertEquals(state.getReadOffset(), read.getLength());
        Assert.assertEquals(state.getCurrentCigarElementOffset(), read.numCigarElements());
        Assert.assertEquals(state.getCurrentCigarElement(), null);
        Assert.assertNotNull(state.toString());
    }

    @Test
    public void stepsOncePerReferenceBaseForBamReadWithCigarInCgTag() throws IOException {
        // A BAM cigar field holds at most 65,535 operators; the writer moves a longer cigar into the CG tag and
        // stores a two-element placeholder in its place, which the reader swaps back (and removes the tag) on decode.
        final int matchInsertionPairs = 35_000;
        final StringBuilder cigar = new StringBuilder();
        for (int i = 0; i < matchInsertionPairs; i++) {
            cigar.append("1M1I");
        }
        cigar.append("1M");
        final int readLength = 2 * matchInsertionPairs + 1;
        final int referenceLength = matchInsertionPairs + 1;
        final byte[] bases = new byte[readLength];
        Arrays.fill(bases, (byte) 'A');
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, (byte) 30);
        final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(ArtificialReadUtils.READ_GROUP_ID);
        readGroup.setSample("sample");
        final SAMFileHeader longHeader = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(readGroup);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(longHeader, "longCigar", 0, 1, bases, quals, cigar.toString());

        final File bam = createTempFile("longCigar", ".bam");
        try (final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(longHeader, true, bam)) {
            writer.addAlignment(read.convertToSAMRecord(longHeader));
        }
        final SAMRecord fromBam;
        try (final SamReader reader = SamReaderFactory.makeDefault().open(bam)) {
            fromBam = reader.iterator().next();
        }
        final GATKRead bamRead = new SAMRecordToGATKReadAdapter(fromBam);
        Assert.assertEquals(bamRead.getCigar().numCigarElements(), 2 * matchInsertionPairs + 1, "the cigar should survive the round trip");

        final AlignmentStateMachine state = new AlignmentStateMachine(bamRead);
        int steps = 0;
        while (state.stepForwardOnGenome() != null) {
            steps++;
        }
        Assert.assertEquals(steps, referenceLength);
        Assert.assertEquals(state.getReadOffset(), readLength);
    }
}
