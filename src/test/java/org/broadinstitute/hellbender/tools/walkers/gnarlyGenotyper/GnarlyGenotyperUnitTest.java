package org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GnarlyGenotyperUnitTest extends GATKBaseTest {

    @Test
    public void testGenomicsDBSkipsNonVariantIntervalsWhenOnlyVariantSitesAreEmitted() {
        Assert.assertTrue(genomicsDBOptionsAfterParsing().skipNonVariantIntervals());
    }

    @Test
    public void testGenomicsDBSkipsNoIntervalsWhenAllSitesAreKept() {
        Assert.assertFalse(genomicsDBOptionsAfterParsing("--keep-all-sites").skipNonVariantIntervals());
    }

    /** @return the GenomicsDB options of a GnarlyGenotyper whose command line is parsed but which is not run */
    private static GenomicsDBOptions genomicsDBOptionsAfterParsing(final String... extraArgs) {
        final List<String> args = new ArrayList<>(Arrays.asList("-R", b38_reference_20_21, "-V", "gendb://unused", "-O", "unused.vcf"));
        args.addAll(Arrays.asList(extraArgs));
        final GnarlyGenotyper tool = new GnarlyGenotyper();
        tool.getCommandLineParser().parseArguments(System.err, args.toArray(new String[0]));
        return tool.getGenomicsDBOptions();
    }
}
