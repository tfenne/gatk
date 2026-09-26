package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

public final class CalledHaplotypesUnitTest extends GATKBaseTest {
    private static final VariantContext CALL = new VariantContextBuilder("test", "1", 100, 100, Arrays.asList(Allele.create("A", true), Allele.create("C"))).make();
    private static final Haplotype HAPLOTYPE = new Haplotype("ACGT".getBytes());

    @Test
    public void callsWithHaplotypesAreAccepted() {
        final List<VariantContext> calls = Collections.singletonList(CALL);
        final Set<Haplotype> haplotypes = Collections.singleton(HAPLOTYPE);
        final CalledHaplotypes called = new CalledHaplotypes(calls, haplotypes);
        Assert.assertEquals(called.getCalls(), calls);
        Assert.assertEquals(called.getCalledHaplotypes(), haplotypes);
    }

    @Test
    public void noCallsAndNoHaplotypesAreAccepted() {
        final CalledHaplotypes called = new CalledHaplotypes(Collections.emptyList(), Collections.emptySet());
        Assert.assertTrue(called.getCalls().isEmpty());
        Assert.assertTrue(called.getCalledHaplotypes().isEmpty());
    }

    @Test
    public void callsWithoutHaplotypesAreRejectedWithBothInTheMessage() {
        try {
            new CalledHaplotypes(Collections.singletonList(CALL), Collections.emptySet());
            Assert.fail("expected an IllegalArgumentException");
        } catch (final IllegalArgumentException e) {
            Assert.assertTrue(e.getMessage().startsWith("Calls and calledHaplotypes should both be empty or both not but got calls="), e.getMessage());
            Assert.assertTrue(e.getMessage().contains(CALL.toString()), e.getMessage());
            Assert.assertTrue(e.getMessage().contains("calledHaplotypes=[]"), e.getMessage());
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void haplotypesWithoutCallsAreRejected() {
        new CalledHaplotypes(Collections.emptyList(), Collections.singleton(HAPLOTYPE));
    }
}
