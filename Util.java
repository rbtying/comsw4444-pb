package pb.g3;

import pb.sim.Asteroid;

import java.util.function.Function;

/**
 * Created by rbtying on 10/6/15.
 */
public class Util {
    public static long findArgMin(long start, long end, Function<Long, Double> f) {
        long idx = start;
        double minval = Double.MAX_VALUE;
        for (long l = start; l < end; ++l) {
            double test = f.apply(l);
            if (test < minval) {
                minval = test;
                idx = l;
            }
        }
        return idx;
    }

    public static double findMin(long start, long end, Function<Long, Double> f) {
        long idx = findArgMin(start, end, f);
        return f.apply(idx);
    }

    public static long findArgMax(long start, long end, Function<Long, Double> f) {
        long idx = start;
        double maxval = Double.MIN_VALUE;
        for (long l = start; l < end; ++l) {
            double test = f.apply(l);
            if (test > maxval) {
                maxval = test;
                idx = l;
            }
        }
        return idx;
    }

    public static double findMax(long start, long end, Function<Long, Double> f) {
        long idx = findArgMax(start, end, f);
        return f.apply(idx);
    }

    public static double energy(Asteroid a, long t) {
        double vel = a.orbit.velocityAt(t - a.epoch).magnitude();
        return 0.5 * a.mass * vel * vel;
    }

}
