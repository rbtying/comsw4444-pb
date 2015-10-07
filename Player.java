package pb.g3;

import org.omg.CORBA.DynAnyPackage.Invalid;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;
import pb.sim.Orbit;
import pb.sim.Point;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

public class Player implements pb.sim.Player {

    private static class Key {
        private final double mass;
        private final long epoch;
        private final double orbit_a;
        private final double orbit_b;
        private final double orbit_A;
        private final double orbit_Mo;
        private int hash;

        private static Map<Asteroid, Key> cache = new ConcurrentHashMap<>();

        public static Key factory(Asteroid a) {
            return cache.computeIfAbsent(a, Key::new);
        }

        private Key(Asteroid a) {
            mass = a.mass;
            epoch = a.epoch;
            orbit_a = a.orbit.a;
            orbit_b = a.orbit.b;
            orbit_A = a.orbit.A;
            orbit_Mo = a.orbit.Mo;
            hash = computeHashCode();
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof Key)) {
                return false;
            }

            Key ok = (Key) o;
            return mass == ok.mass && epoch == ok.epoch && orbit_a == ok.orbit_a && orbit_b == ok.orbit_b && orbit_A == ok.orbit_A && orbit_Mo == ok.orbit_Mo;
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

    private class Push implements Comparable<Push>{
        public int asteroid_idx;
        public long push_time;
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

        private void init(int aidx, long pt, double e, double d, double m) {
            asteroid_idx = aidx;
            push_time = pt;
            energy = e;
            direction = d;
            mass = m;
            _simulated = null;
        }

        public Asteroid simulatedAsteroid(Asteroid a[]) {
            if (_simulated == null) {
                _simulated = Asteroid.push(a[asteroid_idx], push_time, energy, direction);
            }
            return _simulated;
        }

        public String toString() {
            return String.format("Push: Asteroid %d at time %d with energy %f in direction %f", asteroid_idx, push_time, energy, direction);
        }

        @Override
        public int compareTo(Push o) {
            return (int) (push_time - o.push_time);
        }
    }

    public static double hohmannTransfer(Asteroid a, double r_b) {
        // transfers Asteroid a to the aphelion distance of b
        // and hopes for the best
        double r_a = a.orbit.positionAt(findPerihelion(a)).magnitude();

        double dv = Math.sqrt(Orbit.GM / r_a) * (Math.sqrt(2 * r_b / (r_a + r_b)) - 1);
        return Math.copySign(a.mass * dv * dv * 0.5, dv);
    }

    public static double reverseHohmannTransfer(Asteroid a, double r_b) {
        double r_a = a.orbit.positionAt(findPerihelion(a)).magnitude();

        double dv = Math.sqrt(Orbit.GM / r_b) * (1 - Math.sqrt(2 * r_a / (r_a + r_b)));
        return Math.copySign(a.mass * dv * dv * 0.5, dv);
    }

	// used to pick asteroid and velocity boost randomly
	private Random random = new Random();

	// current time, time limit
	private long time = -1;
	private long time_limit = -1;
    private long time_skip = -1;

    private PriorityQueue<Push> next_pushes = new PriorityQueue<>();

    private static double min_distance = Double.MAX_VALUE;

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit)
	{
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
	}

    private static Map<Key, Long> cache_perihelion = new ConcurrentHashMap<>();
    private static Map<Key, Long> cache_aphelion = new ConcurrentHashMap<>();

    private static long findPerihelion(Asteroid a) {
        Key k = Key.factory(a);
        Point p = new Point();
        return (cache_perihelion.computeIfAbsent(k, (x) -> Util.findArgMin(0, a.orbit.period(), (y) -> {
            a.orbit.positionAt(y, p);
            return p.magnitude();
        })) + a.epoch) % a.orbit.period();
    }

    private static long findAphelion(Asteroid a) {
        // maximum distance point is half the orbit away from the minimum distance point
        return (findPerihelion(a) + a.orbit.period() / 2) % a.orbit.period();
    }

    private Push moveToRadius(long push_time, int aidx, Asteroid a, double radius) {
        double r_a = a.orbit.positionAt(push_time - a.epoch).magnitude();
        double normalization_energy = reverseHohmannTransfer(a, r_a); // make it a circle?
        double hohmann_energy = hohmannTransfer(a, radius); // instantaneously ellipsize it
        return new Push(a, aidx, push_time, hohmann_energy + normalization_energy);
    }

	// try to push asteroid
	public void play(Asteroid[] asteroids,
	                 double[] energy, double[] direction) {
        ++time;

        if (time_skip > time) {
            return;
        } else {
            time_skip = -1;
        }

        try {
            if (!next_pushes.isEmpty()) {
                Push next_push = next_pushes.peek();
                if (time < next_push.push_time) {
                    return;
                }
                if (time == next_push.push_time) {
                    // apply push
                    energy[next_push.asteroid_idx] = next_push.energy;
                    direction[next_push.asteroid_idx] = next_push.direction;

                    System.out.println("Making push " + next_push);

                    next_pushes.remove();
                    return;
                }
            } else {
                System.out.println("no next push, computing more");

                Integer indexes[] = new Integer[asteroids.length];
                for (int i = 0; i < indexes.length; ++i) {
                    indexes[i] = i;
                }

                Arrays.sort(indexes, (o1, o2) -> {
                    double d = asteroids[o2].mass - asteroids[o1].mass;
                    if (d > 0) {
                        return 1;
                    } else if (d < 0) {
                        return -1;
                    } else {
                        return 0;
                    }
                });

                // grab largest asteroid
                Asteroid largest = asteroids[indexes[0]];

                double r_largest_ph = largest.orbit.positionAt(findPerihelion(largest)).magnitude();
                double r_largest_ap = largest.orbit.positionAt(findAphelion(largest)).magnitude();
                double r_largest = (r_largest_ap + r_largest_ph) / 2;

                Push potential_pushes[] = new Push[indexes.length];
                for (int i = 0; i < potential_pushes.length; ++i) {
                    potential_pushes[i] = null;
                }

                for (int i = 1; i < indexes.length; ++i) {
                    Asteroid a = asteroids[indexes[i]];

                    for (long push_time = time + 1; push_time < time + largest.orbit.period() * 2; ++push_time) {
                        Push next_push = moveToRadius(push_time, indexes[i], a, r_largest);
                        Asteroid a2 = next_push.simulatedAsteroid(asteroids);
                        long t = a2.orbit.period();
                        long collision_time = push_time + t / 2;

                        Point potential_collision_point = a2.orbit.positionAt(collision_time - a2.epoch);
                        Point largest_location = largest.orbit.positionAt(collision_time - largest.epoch);

                        if (Point.distance(potential_collision_point, largest_location) < a.radius() + largest.radius()) {
                            if (potential_pushes[i] == null || next_push.energy < potential_pushes[i].energy) {
                                potential_pushes[i] = next_push;
                            }
                        }
                    }

                    if (potential_pushes[i] != null) {
                        next_pushes.add(potential_pushes[i]);
                        System.out.println("computed a next push: " + next_pushes.peek());
                        return;
                    }
                }
                for (int i = 0; i < asteroids.length; ++i) {
                    double ecc = Math.sqrt(1 - Math.pow(asteroids[i].orbit.b / asteroids[i].orbit.a, 2));
                    double r_ph = largest.orbit.positionAt(findPerihelion(largest)).magnitude();
                    double r_ap = largest.orbit.positionAt(findAphelion(largest)).magnitude();
                    double r = (r_ap + r_ph) / 2;
                    if (ecc > 0.05) {
                        Push reversePush = new Push(asteroids[i], i, time + 1, reverseHohmannTransfer(largest, r));
                        next_pushes.add(reversePush);
                    } else {
                        time_skip = largest.orbit.period();
                    }
                }
            }
        } catch (ArrayIndexOutOfBoundsException ex) {
            next_pushes.clear();
        }
    }
}
