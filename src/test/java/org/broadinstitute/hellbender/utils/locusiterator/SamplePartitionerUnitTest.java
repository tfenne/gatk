package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public final class SamplePartitionerUnitTest extends GATKBaseTest {
    private static final String READ_GROUP = "rg";
    private static final String SAMPLE = "sample";

    private static SAMFileHeader headerWithSample() {
        final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(READ_GROUP);
        readGroup.setSample(SAMPLE);
        return ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(readGroup);
    }

    /** Batches of reads, one batch per alignment start, all in the one sample. */
    private static List<List<GATKRead>> batches(final SAMFileHeader header, final int batchCount, final int readsPerBatch) {
        final List<List<GATKRead>> batches = new ArrayList<>();
        for (int b = 0; b < batchCount; b++) {
            final List<GATKRead> batch = new ArrayList<>();
            for (int i = 0; i < readsPerBatch; i++) {
                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "b" + b + "r" + i, 0, 100 * (b + 1), 50);
                read.setReadGroup(READ_GROUP);
                batch.add(read);
            }
            batches.add(batch);
        }
        return batches;
    }

    /**
     * Runs one partitioner cycle per batch and returns what each cycle yields for the sample. Between batches,
     * runs the given number of cycles in which nothing is submitted, asserting that each yields no reads.
     */
    private static List<List<GATKRead>> collect(final SamplePartitioner partitioner, final List<List<GATKRead>> batches, final int emptyCyclesBetweenBatches) {
        final List<List<GATKRead>> collected = new ArrayList<>();
        for (final List<GATKRead> batch : batches) {
            for (int cycle = 0; cycle < emptyCyclesBetweenBatches; cycle++) {
                partitioner.doneSubmittingReads();
                Assert.assertTrue(partitioner.getReadsForSample(SAMPLE).isEmpty());
                partitioner.reset();
            }
            batch.forEach(partitioner::submitRead);
            partitioner.doneSubmittingReads();
            collected.add(new ArrayList<>(partitioner.getReadsForSample(SAMPLE)));
            partitioner.reset();
        }
        return collected;
    }

    @Test
    public void emptyCyclesDoNotChangeWhatReservoirDownsamplingKeeps() {
        final SAMFileHeader header = headerWithSample();
        final List<List<GATKRead>> batches = batches(header, 4, 30);
        final LIBSDownsamplingInfo downsampleToFive = new LIBSDownsamplingInfo(true, 5);

        Utils.resetRandomGenerator();
        final List<List<GATKRead>> withoutEmptyCycles = collect(new SamplePartitioner(downsampleToFive, Collections.singletonList(SAMPLE), header), batches, 0);
        Utils.resetRandomGenerator();
        final List<List<GATKRead>> withEmptyCycles = collect(new SamplePartitioner(downsampleToFive, Collections.singletonList(SAMPLE), header), batches, 3);

        for (final List<GATKRead> kept : withoutEmptyCycles) {
            Assert.assertEquals(kept.size(), 5);
        }
        Assert.assertEquals(withEmptyCycles, withoutEmptyCycles);
    }

    @Test
    public void emptyCyclesDoNotChangeWhatPassThroughPartitioningYields() {
        final SAMFileHeader header = headerWithSample();
        final List<List<GATKRead>> batches = batches(header, 3, 4);
        final SamplePartitioner partitioner = new SamplePartitioner(LocusIteratorByState.NO_DOWNSAMPLING, Collections.singletonList(SAMPLE), header);
        Assert.assertEquals(collect(partitioner, batches, 2), batches);
    }
}
