package org.broadinstitute.hellbender.tools.genomicsdb;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.genomicsdb.model.GenomicsDBExportConfiguration;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

public class GenomicsDBOptionsUnitTest extends GATKBaseTest {

    @Test
    public void testSkipsNonVariantIntervalsByDefaultWhenOnlyVariantSitesAreNeeded() {
        Assert.assertTrue(optionsFor(new GenomicsDBArgumentCollection(), true).skipNonVariantIntervals());
    }

    @Test
    public void testDoesNotSkipNonVariantIntervalsWhenTheArgumentIsFalse() {
        final GenomicsDBArgumentCollection args = new GenomicsDBArgumentCollection();
        args.skipNonVariantIntervals = false;
        Assert.assertFalse(optionsFor(args, true).skipNonVariantIntervals());
    }

    @Test
    public void testDoesNotSkipNonVariantIntervalsWhenNonVariantSitesAreNeeded() {
        Assert.assertFalse(optionsFor(new GenomicsDBArgumentCollection(), false).skipNonVariantIntervals());
    }

    @Test
    public void testDoesNotSkipNonVariantIntervalsWhenTheToolDoesNotSayOnlyVariantSitesAreNeeded() {
        Assert.assertFalse(new GenomicsDBOptions().skipNonVariantIntervals());
        Assert.assertFalse(new GenomicsDBOptions(reference(), new GenomicsDBArgumentCollection()).skipNonVariantIntervals());
        Assert.assertFalse(new GenomicsDBOptions(reference(), new GenomicsDBArgumentCollection(),
                new GenotypeCalculationArgumentCollection()).skipNonVariantIntervals());
    }

    @Test
    public void testExportConfigurationSkipsBothKindsOfNonVariantIntervalWhenSkipping() {
        final GenomicsDBExportConfiguration.ExportConfiguration config =
                exportConfigurationFor(optionsFor(new GenomicsDBArgumentCollection(), true));

        Assert.assertTrue(config.getSkipReferenceOnlyIntervals());
        Assert.assertTrue(config.getSkipSpanningDeletionOnlyIntervals());
    }

    @Test
    public void testExportConfigurationSkipsNoIntervalsWhenNotSkipping() {
        final GenomicsDBExportConfiguration.ExportConfiguration config =
                exportConfigurationFor(optionsFor(new GenomicsDBArgumentCollection(), false));

        Assert.assertFalse(config.getSkipReferenceOnlyIntervals());
        Assert.assertFalse(config.getSkipSpanningDeletionOnlyIntervals());
    }

    @Test
    public void testSkippingLeavesOnlyTheOtherRecordsOfAGenomicsDBQuery() {
        final SimpleInterval interval = new SimpleInterval("chr20", 17960187, 17981445);
        final List<File> gvcfs = Arrays.asList(new File(largeFileTestDir + "gvcfs/HG00096.g.vcf.gz"),
                new File(largeFileTestDir + "gvcfs/HG00268.g.vcf.gz"), new File(largeFileTestDir + "gvcfs/NA19625.g.vcf.gz"));
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(GenomicsDBTestUtils.createTempGenomicsDB(gvcfs, interval));

        final List<VariantContext> allRecords = GenomicsDBTestUtils.readGenomicsDBRecords(genomicsDBUri, interval,
                new GenomicsDBOptions(reference()));
        final List<VariantContext> skippingRecords = GenomicsDBTestUtils.readGenomicsDBRecords(genomicsDBUri, interval,
                optionsFor(new GenomicsDBArgumentCollection(), true));

        Assert.assertTrue(allRecords.stream().anyMatch(GenomicsDBTestUtils::isReferenceOnly));
        Assert.assertTrue(allRecords.stream().anyMatch(GenomicsDBTestUtils::isSpanningDeletionOnly));
        final List<String> otherRecords = allRecords.stream()
                .filter(vc -> !GenomicsDBTestUtils.isReferenceOnly(vc) && !GenomicsDBTestUtils.isSpanningDeletionOnly(vc))
                .map(GenomicsDBOptionsUnitTest::siteAndAlleles).toList();
        Assert.assertFalse(otherRecords.isEmpty());
        Assert.assertEquals(skippingRecords.stream().map(GenomicsDBOptionsUnitTest::siteAndAlleles).toList(), otherRecords);
    }

    private static String siteAndAlleles(final VariantContext vc) {
        return vc.getContig() + ":" + vc.getStart() + "-" + vc.getEnd() + " " + vc.getAlleles();
    }

    private static GenomicsDBOptions optionsFor(final GenomicsDBArgumentCollection args, final boolean onlyVariantSitesNeeded) {
        return new GenomicsDBOptions(reference(), args, new GenotypeCalculationArgumentCollection(), onlyVariantSitesNeeded);
    }

    private static Path reference() {
        return new File(b38_reference_20_21).toPath();
    }

    private GenomicsDBExportConfiguration.ExportConfiguration exportConfigurationFor(final GenomicsDBOptions options) {
        final File workspace = createTempDir("GenomicsDBOptionsUnitTest");
        return GATKGenomicsDBUtils.createExportConfiguration(workspace.getAbsolutePath(),
                new File(workspace, "callset.json").getAbsolutePath(), new File(workspace, "vidmap.json").getAbsolutePath(),
                new File(workspace, "vcfheader.vcf").getAbsolutePath(), options);
    }
}
