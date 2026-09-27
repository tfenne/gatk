package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/** Reference definitions the reference confidence tests check the model's incremental computations against. */
final class ReferenceConfidenceTestUtils {
    private ReferenceConfidenceTestUtils() {}

    /**
     * The offset of a pileup element's position into its read aligned one-to-one with the reference: insertions
     * collapsed, deletions padded and soft clips counted.
     */
    static int referenceAlignedOffset(final PileupElement element) {
        final GATKRead read = element.getRead();
        int offset = alignedLength(element.getCurrentCigarElement()) > 0 ? element.getOffsetInCurrentCigar() : 0;
        for (int i = 0; i < element.getCurrentCigarOffset(); i++) {
            offset += alignedLength(read.getCigarElement(i));
        }
        return offset;
    }

    private static int alignedLength(final CigarElement element) {
        final CigarOperator op = element.getOperator();
        return op.consumesReferenceBases() || op == CigarOperator.S ? element.getLength() : 0;
    }
}
