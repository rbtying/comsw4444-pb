package pb.g3;

import pb.sim.Asteroid;
import pb.sim.Orbit;
import pb.sim.Point;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;

/**
 * Created by rbtying on 10/6/15.
 */
public class Util {

    public static final double TWO_PI = 2 * Math.PI;

    public static double normalizedAngle(double angle) {
        return angle - TWO_PI * Math.floor((angle + Math.PI) / TWO_PI);
    }

    /**
     * Creates a string with a time in years and days
     *
     * @param t time in days
     * @return stringified form
     */
    public static String toYearString(long t) {
        return "day " + (t % 365) + " of year " + ((t / 365) + 1);
    }

    /**
     * Finds the next t_offset in the cycle after current_time
     *
     * @param t_offset     the time from the original orbit
     * @param a            the asteroid whose period we are inspecting
     * @param current_time the minimum time
     * @return
     */
    public static long nextAfterTime(long t_offset, Asteroid a, long current_time) {
        long t = current_time % a.orbit.period() + t_offset;
        while (t < current_time) {
            t += a.orbit.period();
        }
        return t;
    }

    /**
     * Returns the time of the collision between asteroids between [initial_time, max_time)
     * runs in O(shorter_orbit_period + num_intersection_points)
     *
     * @param a            the first asteroid to collide
     * @param b            the second asteroid to collide
     * @param initial_time time to start searching from
     * @param max_time     time to stop search at
     * @return -1 if collision not found, otherwise, time of collision
     */
    public static long findCollision(Asteroid a, Asteroid b, long initial_time, long max_time) {
        OrbitPair op = OrbitPair.factory(a, b);

        Asteroid shorter_orbit = a;
        Asteroid longer_orbit = b;
        if (b.orbit.period() < a.orbit.period()) {
            shorter_orbit = b;
            longer_orbit = a;
        }

        Point p1 = new Point(), p2 = new Point();

        if (max_time - initial_time < shorter_orbit.orbit.period()) {
            for (long t = initial_time; t < max_time; ++t) {
                positionAt(shorter_orbit, t, p1);
                positionAt(longer_orbit, t, p2);
                double dist = Point.distance(p1, p2);
                if (dist < op.collision_distance) {
                    return t;
                }
            }
        }

        List<Long> potential_intersection_points;

        try {
            potential_intersection_points = op.intersections();
        } catch (OutOfMemoryError e) {
            OrbitPair.intersection_cache.clear();
            potential_intersection_points = op.intersections();
        }

        // potential_intersection_points now contains up to num_intersections points to check
        // points in potential_intersection_points are not epoch-adjusted
        long base_offset = initial_time + shorter_orbit.epoch;
        long T = op.shorter_orbit.period();

        while (base_offset < max_time) {
            for (Long intersection_point : potential_intersection_points) {
                long t = base_offset + intersection_point;
                positionAt(shorter_orbit, t, p1);
                positionAt(longer_orbit, t, p2);
                double dist = Point.distance(p1, p2);
                if (dist < op.collision_distance) {
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
     * @param a the asteroid
     * @return its eccentricity
     */
    public static double eccentricity(Asteroid a) {
        return Math.sqrt(1 - Math.pow(a.orbit.b / a.orbit.a, 2));
    }

    /**
     * Computes the argmax of f from i=start to i=end
     *
     * @param start initial value of long
     * @param end   one more than last value of long
     * @param f     f: long -> double
     * @return value of long that maximizes f
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
     * @return value of long which minimizes f
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
     * @return minimum value of f
     */
    public static double findMin(long start, long end, Function<Long, Double> f) {
        long idx = findArgMin(start, end, f);
        return f.apply(idx);
    }


    /**
     * Computes the max of f from i=start to i=end
     *
     * @param start initial value of long
     * @param end   one more than last value of long
     * @param f     f: long -> double
     * @return maximum value of f
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

    public static boolean eq(double a, double b) {
        return Double.isNaN(a) && Double.isNaN(b) || a == b;
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

    public static class OrbitPair {
        private static Map<OrbitPair, List<Long>> intersection_cache = new ConcurrentHashMap<>();
        public Orbit shorter_orbit, longer_orbit;
        public double collision_distance;

        private OrbitPair(Orbit a, Orbit b, double collision_distance) {
            if (a.period() > b.period()) {
                shorter_orbit = b;
                longer_orbit = a;
            } else {
                shorter_orbit = a;
                longer_orbit = b;
            }
            this.collision_distance = collision_distance;
        }

        public static void clearCache() {
            System.out.println("Deleting " + intersection_cache.size() + " items from intersection cache");
            intersection_cache.clear();
        }

        public static boolean orbitsEqual(Orbit a, Orbit b) {
            return eq(a.a, b.a) && eq(a.b, b.b) && eq(a.A, b.A) && eq(a.Mo, b.Mo);
        }

        public static OrbitPair factory(Asteroid a, Asteroid b) {
            return new OrbitPair(a.orbit, b.orbit, a.radius() + b.radius());
        }

        @Override
        public int hashCode() {
            int hash = 17;
            hash = hash * 23 + Double.hashCode(longer_orbit.a);
            hash = hash * 23 + Double.hashCode(longer_orbit.b);
            hash = hash * 23 + Double.hashCode(longer_orbit.A);
            hash = hash * 23 + Double.hashCode(longer_orbit.Mo);
            hash = hash * 23 + Double.hashCode(shorter_orbit.a);
            hash = hash * 23 + Double.hashCode(shorter_orbit.b);
            hash = hash * 23 + Double.hashCode(shorter_orbit.A);
            hash = hash * 23 + Double.hashCode(shorter_orbit.Mo);
            hash = hash * 23 + Double.hashCode(collision_distance);
            return hash;
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof OrbitPair)) {
                return false;
            }
            OrbitPair op = (OrbitPair) o;
            boolean eq_s = orbitsEqual(shorter_orbit, shorter_orbit);
            boolean eq_l = orbitsEqual(longer_orbit, longer_orbit);
            boolean coll_eq = collision_distance == op.collision_distance;

            return eq_s && eq_l && coll_eq;
        }

        /**
         * Computes intersections between orbits. Cached.
         *
         * @return List of intersections, denoted in positions on the shorter orbit
         */
        public List<Long> intersections() {
            return intersection_cache.computeIfAbsent(this, (orbitPair) -> {
                List<Long> intersections = new ArrayList<>(2000);
                Point focus = orbitPair.longer_orbit.center();
                focus.x *= 2;
                focus.y *= 2;

                Point p = new Point();
                for (long i = 0; i < orbitPair.shorter_orbit.period(); ++i) {
                    orbitPair.shorter_orbit.positionAt(i, p);
                    double d1 = p.magnitude();
                    double d2 = Point.distance(p, focus);

                    if (d1 + d2 - 2 * longer_orbit.a < orbitPair.collision_distance) {
                        intersections.add(i);
                    }
                }
                return intersections;
            });
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
            double v_a_mag = Math.sqrt((2 * a.energy) / a.mass);
            double v_b_mag = Math.sqrt((2 * b.energy) / b.mass);

            double dx = Math.cos(a.direction) * v_a_mag + Math.cos(b.direction) * v_b_mag;
            double dy = Math.sin(a.direction) * v_a_mag + Math.sin(b.direction) * v_b_mag;

            double E = 0.5 * a.mass * (dx * dx + dy * dy);
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
            return String.format("Push: asteroid %d on %s with energyAtTime %6.3e in direction %.02f deg",
                    asteroid_idx, toYearString(push_time), energy, Math.toDegrees(direction));
        }

        @Override
        public int compareTo(Push o) {
            return (int) (push_time - o.push_time);
        }
    }

}
