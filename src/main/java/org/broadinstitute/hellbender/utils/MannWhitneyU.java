package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


import java.util.*;

/**
 * Imported with changes from Picard private.
 *
 * @author Tim Fennell
 */
public class MannWhitneyU {

    protected static Logger logger = LogManager.getLogger(MannWhitneyU.class);

    private static final class Rank implements Comparable<Rank> {
        final double value;
        float rank;
        final int series;

        private Rank(double value, float rank, int series) {
            this.value = value;
            this.rank = rank;
            this.series = series;
        }

        @Override
        public int compareTo(Rank that) {
            return (int) (this.value - that.value);
        }

        @Override
        public String toString() {
            return "Rank{" +
                    "value=" + value +
                    ", rank=" + rank +
                    ", series=" + series +
                    '}';
        }
    }

    /**
     * The results of performing a rank sum test.
     */
    public static class Result {
        private final double u;
        private final double z;
        private final double p;
        private final double medianShift;

        public Result(double u, double z, double p, double medianShift) {
            this.u = u;
            this.z = z;
            this.p = p;
            this.medianShift = medianShift;
        }

        public double getU() {
            return u;
        }

        public double getZ() {
            return z;
        }

        public double getP() {
            return p;
        }

        public double getMedianShift() {
            return medianShift;
        }
    }

    /**
     * The values of U1, U2 and the transformed number of ties needed for the calculation of sigma
     * in the normal approximation.
     */
    public static class TestStatistic {
        private final double u1;
        private final double u2;
        private final double trueU;
        private final double numOfTiesTransformed;

        public TestStatistic(double u1, double u2, double numOfTiesTransformed) {
            this.u1 = u1;
            this.u2 = u2;
            this.numOfTiesTransformed = numOfTiesTransformed;
            this.trueU = Double.NaN;
        }

        public TestStatistic(double trueU, double numOfTiesTransformed) {
            this.trueU = trueU;
            this.numOfTiesTransformed = numOfTiesTransformed;
            this.u1 = Double.NaN;
            this.u2 = Double.NaN;
        }

        public double getU1() {
            return u1;
        }

        public double getU2() {
            return u2;
        }

        public double getTies() {
            return numOfTiesTransformed;
        }

        public double getTrueU() {
            return trueU;
        }
    }

    /**
     * The ranked data in one list and a list of the number of ties.
     */
    public static class RankedData {
        private final Rank[] rank;
        private final ArrayList<Integer> numOfTies;

        public RankedData(Rank[] rank, ArrayList<Integer> numOfTies) {
            this.rank = rank;
            this.numOfTies = numOfTies;
        }

        public Rank[] getRank() {
            return rank;
        }

        public ArrayList<Integer> getNumOfTies() {
            return numOfTies;
        }
    }

    // Constructs a normal distribution; this needs to be a standard normal in order to get a Z-score in the exact case
    private static final double NORMAL_MEAN = 0;
    private static final double NORMAL_SD = 1;
    private static final NormalDistribution NORMAL = new NormalDistribution(NORMAL_MEAN, NORMAL_SD);

    /**
     * The minimum length for both data series in order to use a normal distribution
     * to calculate Z and p. If both series are shorter than this value then a permutation test
     * will be used.
     */
    private static final int MINIMUM_NORMAL_N = 10;

    /**
     * A variable that indicates if the test is one sided or two sided and if it's one sided
     * which group is the dominator in the null hypothesis.
     */
    public enum TestType {
        FIRST_DOMINATES,
        SECOND_DOMINATES,
        TWO_SIDED
    }

    public RankedData calculateRank(final double[] series1, final double[] series2) {
        Arrays.sort(series1);
        Arrays.sort(series2);

        // Make a merged ranks array
        final Rank[] ranks = new Rank[series1.length + series2.length];
        {
            int i = 0, j = 0, r = 0;
            while (r < ranks.length) {
                if (i >= series1.length) {
                    ranks[r++] = new Rank(series2[j++], r, 2);
                } else if (j >= series2.length) {
                    ranks[r++] = new Rank(series1[i++], r, 1);
                } else if (series1[i] <= series2[j]) {
                    ranks[r++] = new Rank(series1[i++], r, 1);
                } else {
                    ranks[r++] = new Rank(series2[j++], r, 2);
                }
            }
        }

        ArrayList<Integer> numOfTies = new ArrayList<>();

        // Now sort out any tie bands
        for (int i = 0; i < ranks.length; ) {
            float rank = ranks[i].rank;
            int count = 1;

            for (int j = i + 1; j < ranks.length && ranks[j].value == ranks[i].value; ++j) {
                rank += ranks[j].rank;
                ++count;
            }

            if (count > 1) {
                rank /= count;
                for (int j = i; j < i + count; ++j) {
                    ranks[j].rank = rank;
                }
                numOfTies.add(count);
            }

            // Skip forward the right number of items
            i += count;
        }

        return new RankedData(ranks, numOfTies);
    }

    /**
     * Rank both groups together and return a TestStatistic object that includes U1, U2 and number of ties for sigma
     */
    public TestStatistic calculateU1andU2(final double[] series1, final double[] series2) {
        RankedData ranked = calculateRank(series1, series2);
        Rank[] ranks = ranked.getRank();
        ArrayList<Integer> numOfTies = ranked.getNumOfTies();
        int lengthOfRanks = ranks.length;

        double numOfTiesForSigma = transformTies(lengthOfRanks, numOfTies);

        // Calculate R1 and R2 and U.
        float r1 = 0, r2 = 0;
        for (Rank rank : ranks) {
            if (rank.series == 1) r1 += rank.rank;
            else r2 += rank.rank;
        }

        double n1 = series1.length;
        double n2 = series2.length;
        double u1 = r1 - ((n1 * (n1 + 1)) / 2);
        double u2 = r2 - ((n2 * (n2 + 1)) / 2);

        TestStatistic result = new TestStatistic(u1, u2, numOfTiesForSigma);
        return result;
    }

    public double transformTies(int numOfRanks, ArrayList<Integer> numOfTies) {
        //Calculate number of ties transformed for formula for Sigma to calculate Z-score
        ArrayList<Double> transformedTies = new ArrayList<>();
        for (int count : numOfTies) {
            //If every single datapoint is tied then we want to return a p-value of .5 and
            //the formula for sigma that includes the number of ties breaks down. Setting
            //the number of ties to 0 in this case gets the desired result in the normal
            //approximation case.
            if (count != numOfRanks) {
                transformedTies.add((Math.pow(count, 3)) - count);
            }
        }

        double numOfTiesForSigma = 0.0;
        for (double count : transformedTies) {
            numOfTiesForSigma += count;
        }

        return(numOfTiesForSigma);
    }
    /**
     * Calculates the rank-sum test statisic U (sometimes W) from two sets of input data for a one-sided test
     * with an int indicating which group is the dominator. Returns a test statistic object with trueU and number of
     * ties for sigma.
     */
    public TestStatistic calculateOneSidedU(final double[] series1, final double[] series2, final TestType whichSeriesDominates) {
        TestStatistic stat = calculateU1andU2(series1, series2);
        TestStatistic result;
        if (whichSeriesDominates == TestType.FIRST_DOMINATES) {
            result = new TestStatistic(stat.getU1(), stat.getTies());
        } else {
            result = new TestStatistic(stat.getU2(), stat.getTies());
        }
        return result;
    }

    /**
     * Calculates the two-sided rank-sum test statisic U (sometimes W) from two sets of input data.
     * Returns a test statistic object with trueU and number of ties for sigma.
     */
    public TestStatistic calculateTwoSidedU(final double[] series1, final double[] series2) {
        TestStatistic u1AndU2 = calculateU1andU2(series1, series2);
        double u = Math.min(u1AndU2.getU1(), u1AndU2.getU2());
        TestStatistic result = new TestStatistic(u, u1AndU2.getTies());
        return result;
    }

    /**
     * Calculates the Z score (i.e. standard deviations from the mean) of the rank sum
     * test statistics given input data of lengths n1 and n2 respectively, as well as the number of ties, for normal
     * approximation only.
     */
    public double calculateZ(final double u, final int n1, final int n2, final double nties, final TestType whichSide) {
        double m = (n1 * n2) / 2d;

        //Adds a continuity correction
        double correction;
        if (whichSide == TestType.TWO_SIDED) {
            correction = (u - m) >= 0 ? .5 : -.5;
        } else {
            correction = whichSide == TestType.FIRST_DOMINATES ? -.5 : .5;
        }

        //If all the data is tied, the number of ties for sigma is set to 0. In order to get a p-value of .5 we need to
        //remove the continuity correction.
        if (nties == 0) {
            correction = 0;
        }

        double sigma = Math.sqrt((n1 * n2 / 12d) * ((n1 + n2 + 1) - nties / ((n1 + n2) * (n1 + n2 - 1))));
        return (u - m - correction) / sigma;
    }

    /**
     * Finds or calculates the median value of a sorted array of double.
     */
    public double median(final double[] data) {
        final int len = data.length;
        final int mid = len / 2;
        if (data.length % 2 == 0) {
            return (data[mid] + data[mid - 1]) / 2d;
        } else {
            return data[mid];
        }
    }


    /**
     * Constructs a new rank sum test with the given data.
     *
     * @param series1   group 1 data
     * @param series2   group 2 data
     * @param whichSide indicator of two sided test, 0 for two sided, 1 for series1 as dominator, 2 for series2 as dominator
     * @return Result including U statistic, Z score, p-value, and difference in medians.
     */
    public Result test(final double[] series1, final double[] series2, final TestType whichSide) {
        final int n1 = series1.length;
        final int n2 = series2.length;

        //If one of the groups is empty we return NaN
        if (n1 == 0 || n2 == 0) {
            return new Result(Float.NaN, Float.NaN, Float.NaN, Float.NaN);
        }

        double u;
        double nties;

        if (whichSide == TestType.TWO_SIDED) {
            TestStatistic result = calculateTwoSidedU(series1, series2);
            u = result.getTrueU();
            nties = result.getTies();
        } else {
            TestStatistic result = calculateOneSidedU(series1, series2, whichSide);
            u = result.getTrueU();
            nties = result.getTies();
        }

        double z;
        double p;

        if (n1 >= MINIMUM_NORMAL_N || n2 >= MINIMUM_NORMAL_N) {
            z = calculateZ(u, n1, n2, nties, whichSide);
            p = 2 * NORMAL.cumulativeProbability(NORMAL_MEAN + z * NORMAL_SD);
            if (whichSide != TestType.TWO_SIDED) {
                p = p / 2;
            }
        } else {
            // TODO -- This exact test is only implemented for the one sided test, but we currently don't call the two sided version
            if (whichSide != TestType.FIRST_DOMINATES) {
                logger.warn("An exact two-sided MannWhitneyU test was called. Only the one-sided exact test is implemented, use the approximation instead by setting minimumNormalN to 0.");
            }
            p = permutationTest(series1, series2, u);
            z = NORMAL.inverseCumulativeProbability(p);
        }

        return new Result(u, z, p, Math.abs(median(series1) - median(series2)));
    }

    /**
     * Computes the exact one-sided p-value of the observed U over every assignment of the pooled ranks to groups of
     * the observed sizes: the fraction of assignments whose group-1 U is smaller than the observed U, counting
     * assignments with the observed U by half. An observed U that no assignment gives group 1 is counted the same
     * way, so it never needs a matching bin. Counting the observed bin by half deals with the edge case where the observed value
     * is also the most extreme value; the plain cumulative distribution gives a p-value of 1 there, which doesn't
     * result in a usable z-score.
     *
     * @param series1 Data from group 1
     * @param series2 Data from group 2
     * @param testStatU Test statistic U from observed data
     * @return P-value of the observed U over all assignments of the ranks to the two groups
     */
    private double permutationTest(final double[] series1, final double[] series2, final double testStatU) {
        final int n1 = series1.length;
        final Rank[] ranks = calculateRank(series1, series2).getRank();

        // Tied ranks are averaged within their band, so every rank is a multiple of 1/2 and twice it is an exact
        // integer; so is twice U, which is twice the group-1 rank sum less n1 * (n1 + 1).
        final int[] doubledRanks = new int[ranks.length];
        int maxDoubledRankSum = 0;
        for (int i = 0; i < ranks.length; i++) {
            doubledRanks[i] = Math.round(2 * ranks[i].rank);
            maxDoubledRankSum += doubledRanks[i];
        }

        // assignments[k][s]: the number of ways to choose k of the ranks seen so far with doubled ranks summing to s
        final long[][] assignments = new long[n1 + 1][maxDoubledRankSum + 1];
        assignments[0][0] = 1;
        for (final int doubledRank : doubledRanks) {
            for (int k = n1; k >= 1; k--) {
                for (int sum = maxDoubledRankSum; sum >= doubledRank; sum--) {
                    assignments[k][sum] += assignments[k - 1][sum - doubledRank];
                }
            }
        }

        final long observedDoubledU = Math.round(2 * testStatU);
        final long doubledUOffset = (long) n1 * (n1 + 1);
        long smaller = 0;
        long equal = 0;
        long total = 0;
        for (int sum = 0; sum <= maxDoubledRankSum; sum++) {
            final long count = assignments[n1][sum];
            final long doubledU = sum - doubledUOffset;
            total += count;
            if (doubledU < observedDoubledU) {
                smaller += count;
            } else if (doubledU == observedDoubledU) {
                equal += count;
            }
        }
        return (equal / 2.0 + smaller) / total;
    }

}
