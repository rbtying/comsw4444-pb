package pb.g3;

import pb.sim.Asteroid;
import pb.sim.Point;

import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;

/**
 * Created by rbtying on 10/6/15.
 */
public class Util {

    public static long nextAfterTime(long t_offset, Asteroid a, long current_time) {
        long t = current_time % a.orbit.period() + t_offset;
        if (t < current_time) {
            t += a.orbit.period();
        }
        return t;
    }

    /**
     * Returns the time of the collision between asteroids between [initial_time, max_time)
     * <p>
     * runs in O(shorter_orbit_period + num_intersection_points)
     *
     * @param a            the first asteroid to collide
     * @param b            the second asteroid to collide
     * @param initial_time
     * @param max_time
     * @return -1 if collision not found, otherwise, time of collision
     */
    public static long findCollision(Asteroid a, Asteroid b, long initial_time, long max_time) {
        Asteroid shorter_orbit = a;
        Asteroid longer_orbit = b;
        if (b.orbit.period() < a.orbit.period()) {
            shorter_orbit = b;
            longer_orbit = a;
        }

        // if more than 5 points intersect per real intersection, something is very wrong
        ArrayList<Long> potential_intersection_points = new ArrayList<>(512);

        Point p1 = new Point(), p2 = new Point(), focus = new Point();
        longer_orbit.orbit.center(focus);
        focus.x *= 2;
        focus.y *= 2;

        double collision_distance = shorter_orbit.radius() + longer_orbit.radius();

        long T = shorter_orbit.orbit.period();

        for (long time = 0; time < T && initial_time + time < max_time; ++time) {
            long t = time + initial_time;
            positionAt(shorter_orbit, t, p1);
            double d1 = Point.distance(p1, focus);
            double r = p1.magnitude();
            if (Math.abs(d1 + r - 2 * longer_orbit.orbit.a) < collision_distance) {
                potential_intersection_points.add(time);
            }
        }

        // potential_intersection_points now contains up to num_intersections points to check
        // points in potential_intersection_points are epoch-adjusted
        long base_offset = initial_time;

        while (base_offset < max_time) {
            for (Long intersection_point : potential_intersection_points) {
                long t = base_offset + intersection_point;
                positionAt(shorter_orbit, t, p1);
                positionAt(longer_orbit, t, p2);
                double dist = Point.distance(p1, p2);
                if (dist < collision_distance) {
                    return t;
                }
            }
            // check in the next period
            base_offset += T;
        }

        return -1;
    }

    /**
     * Computes the eccentricity of an asteroid
     *
     * @param a
     * @return
     */
    public static double eccentricity(Asteroid a) {
        return Math.sqrt(1 - Math.pow(a.orbit.b / a.orbit.a, 2));
    }

    /**
     * Computes the argmax of f from i=start to i=end
     * @param start initial value of long
     * @param end one more than last value of long
     * @param f f: long -> double
     * @return
     */
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

    /**
     * Computes the argmin of f from i=start to i=end
     *
     * @param start initial value of long
     * @param end   one more than last value of long
     * @param f     f: long -> double
     * @return
     */
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

    /**
     * Computes the min of f from i=start to i=end
     *
     * @param start initial value of long
     * @param end   one more than last value of long
     * @param f     f: long -> double
     * @return
     */
    public static double findMin(long start, long end, Function<Long, Double> f) {
        long idx = findArgMin(start, end, f);
        return f.apply(idx);
    }


    /**
     * Computes the max of f from i=start to i=end
     * @param start initial value of long
     * @param end one more than last value of long
     * @param f f: long -> double
     * @return
     */
    public static double findMax(long start, long end, Function<Long, Double> f) {
        long idx = findArgMax(start, end, f);
        return f.apply(idx);
    }

    /**
     * Computes the epoch-adjusted velocity of an asteroid at a known time
     *
     * @param a the asteroid
     * @param t the time
     * @return epoch-adjusted velocity
     */
    public static Point velocityAt(Asteroid a, long t) {
        return a.orbit.velocityAt(t - a.epoch);
    }

    public static void velocityAt(Asteroid a, long t, Point p) {
        a.orbit.velocityAt(t - a.epoch, p);
    }

    /**
     * Computes the epoch-adjusted position of an asteroid at a known time
     *
     * @param a the asteroid
     * @param t the time
     * @return epoch-adjusted position
     */
    public static Point positionAt(Asteroid a, long t) {
        return a.orbit.positionAt(t - a.epoch);
    }

    public static void positionAt(Asteroid a, long t, Point p) {
        a.orbit.positionAt(t - a.epoch, p);
    }

    /**
     * Computes the epoch-adjusted kinetic energy at a time
     *
     * @param a the asteroid
     * @param t the time
     * @return epoch-adjusted kinetic energy
     */
    public static double energyAtTime(Asteroid a, long t) {
        return energyFromVelocity(a, velocityAt(a, t).magnitude());
    }

    /**
     * Computes the amount of energy necessary to introduce delta V of v on a
     *
     * @param a the asteroid
     * @param v the velocity
     */
    public static double energyFromVelocity(Asteroid a, double v) {
        return 0.5 * a.mass * v * v;
    }

    public static class Key {
        private static Map<Asteroid, Key> cache = new ConcurrentHashMap<>();
        private final double mass;
        private final long epoch;
        private final double orbit_a;
        private final double orbit_b;
        private final double orbit_A;
        private final double orbit_Mo;
        private int hash;

        private Key(Asteroid a) {
            mass = a.mass;
            epoch = a.epoch;
            orbit_a = a.orbit.a;
            orbit_b = a.orbit.b;
            orbit_A = a.orbit.A;
            orbit_Mo = a.orbit.Mo;
            hash = computeHashCode();
        }

        public static Key factory(Asteroid a) {
            return cache.computeIfAbsent(a, Key::new);
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof Key)) {
                return false;
            }

            Key ok = (Key) o;
            return mass == ok.mass && epoch == ok.epoch && orbit_a == ok.orbit_a && orbit_b == ok.orbit_b &&
                    orbit_A == ok.orbit_A && orbit_Mo == ok.orbit_Mo;
        }

        private int computeHashCode() {
            int hash = 17;
            hash = hash * 23 + Double.hashCode(mass);
            hash = hash * 23 + Double.hashCode(epoch);
            hash = hash * 23 + Double.hashCode(orbit_a);
            hash = hash * 23 + Double.hashCode(orbit_b);
            hash = hash * 23 + Double.hashCode(orbit_A);
            hash = hash * 23 + Double.hashCode(orbit_Mo);
            return hash;
        }

        @Override
        public int hashCode() {
            return hash;
        }
    }

    public static class Push implements Comparable<Push> {
        public int asteroid_idx;
        public long push_time;
        public long expected_collision_time;
        public double energy;
        public double direction;
        public double mass;
        private Asteroid _simulated;

        public Push(Asteroid a, int aidx, long pt, double e) {
            double dir = a.orbit.velocityAt(pt - a.epoch).direction();
            if (e < 0) {
                dir += Math.PI;
            }
            init(aidx, pt, Math.abs(e), dir, a.mass);
        }

        public Push(int aidx, long pt, double e, double d, double m) {
            init(aidx, pt, e, d, m);
        }

        public static Push add(Push a, Push b) {
            if (a.asteroid_idx != b.asteroid_idx || a.push_time != b.push_time || a.mass != b.mass) {
                throw new RuntimeException("Adding incompatible pushes");
            }
            double double_root_m_times_v_a = Math.sqrt(a.energy);
            double double_root_m_times_v_b = Math.sqrt(b.energy);
            double dxa = Math.cos(a.direction) * double_root_m_times_v_a;
            double dya = Math.sin(a.direction) * double_root_m_times_v_a;
            double dxb = Math.cos(b.direction) * double_root_m_times_v_b;
            double dyb = Math.sin(b.direction) * double_root_m_times_v_b;
            double dx = dxa + dxb;
            double dy = dya + dyb;

            double E = dx * dx + dy * dy;
            double direction = Math.atan2(dy, dx);

            return new Push(a.asteroid_idx, a.push_time, E, direction, a.mass);
        }

        private void init(int aidx, long pt, double e, double d, double m) {
            asteroid_idx = aidx;
            push_time = pt;
            energy = e;
            direction = d;
            mass = m;
            _simulated = null;
            expected_collision_time = -1;
        }

        public Asteroid simulatedAsteroid(Asteroid a[]) {
            if (_simulated == null) {
                _simulated = Asteroid.push(a[asteroid_idx], push_time, energy, direction);
            }
            return _simulated;
        }

        public String toString() {
            return String.format("Push: Asteroid %d at time %d with energyAtTime %f in direction %f", asteroid_idx,
                    push_time, energy, direction);
        }

        @Override
        public int compareTo(Push o) {
            return (int) (push_time - o.push_time);
        }
    }

}
