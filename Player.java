package pb.g3;

import pb.sim.Asteroid;
import pb.sim.Orbit;
import pb.sim.Point;

import java.util.Arrays;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

public class Player implements pb.sim.Player {

    private static double min_distance = Double.MAX_VALUE;
    private static Map<Util.Key, Long> cache_perihelion = new ConcurrentHashMap<>();
    private static Map<Util.Key, Long> cache_aphelion = new ConcurrentHashMap<>();
    // used to pick asteroid and velocity boost randomly
	private Random random = new Random();
	// current time, time limit
	private long time = -1;
	private long time_limit = -1;
    private long time_skip = -1;
    private int num_asteroids = -1;
    private PriorityQueue<Util.Push> next_pushes = new PriorityQueue<>();

    /**
     * Computes the energyAtTime needed to transfer asteroid a to radius b when
     * launching from its apoapsis
     *
     * @param a   Asteroid to transfer
     * @param r_b Radius to transfer to
     * @return Signed energyAtTime requirement
     */
    public static double hohmannTransfer(Asteroid a, double r_b) {
        // transfers Asteroid a to the aphelion distance of b
        // and hopes for the best
        double r_a = a.orbit.positionAt(findApoapsis(a)).magnitude();

        double dv = Math.sqrt(Orbit.GM / r_a) * (Math.sqrt(2 * r_b / (r_a + r_b)) - 1);
        return Math.copySign(a.mass * dv * dv * 0.5, dv);
    }

    /**
     * Computes the energyAtTime needed to transfer asteroid a from radius b
     *
     * @param a
     * @param r_b
     * @return
     */
    public static double reverseHohmannTransfer(Asteroid a, double r_b) {
        double r_a = a.orbit.positionAt(findPeriapsis(a)).magnitude();

        double dv = Math.sqrt(Orbit.GM / r_b) * (1 - Math.sqrt(2 * r_a / (r_a + r_b)));
        return Math.copySign(a.mass * dv * dv * 0.5, dv);
    }

    private static long findPeriapsis(Asteroid a) {
        Util.Key k = Util.Key.factory(a);
        Point p = new Point();
        return (cache_perihelion.computeIfAbsent(k, (x) -> Util.findArgMin(0, a.orbit.period(), (y) -> {
            a.orbit.positionAt(y, p);
            return p.magnitude();
        })) + a.epoch) % a.orbit.period();
    }

    private static long findApoapsis(Asteroid a) {
        // maximum distance point is half the orbit away from the minimum distance point
        return (findPeriapsis(a) + a.orbit.period() / 2) % a.orbit.period();
    }

    private static double findPeriapsisDistance(Asteroid a) {
        return Util.positionAt(a, findPeriapsis(a)).magnitude();
    }

    private static double findApoapsisDistance(Asteroid a) {
        return Util.positionAt(a, findApoapsis(a)).magnitude();
    }

    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit) {
        if (Orbit.dt() != 24 * 60 * 60) {
            throw new IllegalStateException("Time quantum is not a day");
        }
        this.time_limit = time_limit;
        this.num_asteroids = asteroids.length;
    }

    private Util.Push moveToRadius(long push_time, int aidx, Asteroid a, double radius) {
        double r_a = a.orbit.positionAt(push_time - a.epoch).magnitude();
        double normalization_energy = reverseHohmannTransfer(a, r_a); // make it a circle?
        double hohmann_energy = hohmannTransfer(a, radius); // instantaneously ellipsize it
        return new Util.Push(a, aidx, push_time, hohmann_energy + normalization_energy);
    }

	// try to push asteroid
	public void play(Asteroid[] asteroids,
	                 double[] energy, double[] direction) {
        ++time;

        if (num_asteroids != asteroids.length) {
            next_pushes.clear();
            System.out.println(String.format("#asteroids changed from %d to %d, dropping queued pushes", num_asteroids,
                    asteroids.length));
            num_asteroids = asteroids.length;
        }

        if (time_skip > time) {
            return;
        } else {
            time_skip = -1;
        }
        if (!next_pushes.isEmpty()) {
            Util.Push next_push = next_pushes.peek();
            while (time > next_push.push_time) {
                next_pushes.remove();
                next_push = next_pushes.peek();
                if (next_push == null) {
                    return;
                }
            }

            if (time < next_push.push_time) {
                return;
            }
            if (time == next_push.push_time) {
                // apply push
                energy[next_push.asteroid_idx] = next_push.energy;
                direction[next_push.asteroid_idx] = next_push.direction;

                System.out.println("Making push " + next_push);

                if (time_skip < next_push.expected_collision_time) {
                    time_skip = next_push.expected_collision_time;
                    System.out.println("Waiting until time " + time_skip);
                }

                next_pushes.remove();
                return;
            }
        } else {
            // no pushes computed

            System.out.println("no next push, computing more");

            Integer indexes[] = new Integer[asteroids.length];
            Integer radius_indexes[] = new Integer[asteroids.length];
            for (int i = 0; i < indexes.length; ++i) {
                indexes[i] = i;
                radius_indexes[i] = i;
            }

            Arrays.sort(radius_indexes, (o1, o2) -> (int) Math.signum(Util.positionAt(asteroids[o2], time).magnitude() -
                    Util.positionAt(asteroids[o1], time).magnitude()));

            Arrays.sort(indexes, (o1, o2) -> (int) Math.signum(asteroids[o2].mass - asteroids[o1].mass));

            // grab largest asteroid
            Asteroid largest = asteroids[indexes[0]];

            double r_largest_ph = findPeriapsisDistance(largest);
            double r_largest_ap = findApoapsisDistance(largest);

            PriorityQueue<Util.Push> best_next_push_heap = new PriorityQueue<Util.Push>((p1, p2) -> (int) Math.signum(p1.energy - p2.energy));

            for (int i = 1; i < indexes.length; ++i) {
                Asteroid a = asteroids[indexes[i]];
                long push_time;
                double r_peri = findPeriapsisDistance(a);
                double r_ap = findApoapsisDistance(a);

                if (Math.abs(r_ap - r_peri) < a.radius()) {
                    // this thing is basically a circle
                    for (long dpush_time = 0; dpush_time < a.orbit.period(); ++dpush_time) {
                        push_time = time + dpush_time;
                        Util.Push push_to_r_largest_ph = moveToRadius(push_time, indexes[i], a, r_largest_ph);
                        Asteroid r_largest_ph_a = push_to_r_largest_ph.simulatedAsteroid(asteroids);
                        push_to_r_largest_ph.expected_collision_time = Util.findCollision(largest, r_largest_ph_a, push_time, push_time + 365 * 5);

                        Util.Push push_to_r_largest_ap = moveToRadius(push_time, indexes[i], a, r_largest_ap);
                        Asteroid r_largest_ap_a = push_to_r_largest_ap.simulatedAsteroid(asteroids);
                        push_to_r_largest_ap.expected_collision_time = Util.findCollision(largest, r_largest_ap_a, push_time, push_time + 365 * 5);

                        if (push_to_r_largest_ap.expected_collision_time >= 0) {
                            best_next_push_heap.add(push_to_r_largest_ap);
                        }
                        if (push_to_r_largest_ph.expected_collision_time >= 0) {
                            best_next_push_heap.add(push_to_r_largest_ph);
                        }
                    }
                } else {
                    // only check the apoapsis
                    push_time = Util.nextAfterTime(findApoapsis(a), a, time);

                    Util.Push push_to_r_largest_ph = moveToRadius(push_time, indexes[i], a, r_largest_ph);
                    Asteroid r_largest_ph_a = push_to_r_largest_ph.simulatedAsteroid(asteroids);
                    long collision_time_ph = Util.findCollision(largest, r_largest_ph_a, push_time, push_time + 365 * 5);

                    Util.Push push_to_r_largest_ap = moveToRadius(push_time, indexes[i], a, r_largest_ap);
                    Asteroid r_largest_ap_a = push_to_r_largest_ap.simulatedAsteroid(asteroids);
                    long collision_time_ap = Util.findCollision(largest, r_largest_ap_a, push_time, push_time + 365 * 5);

                    if (collision_time_ap >= 0) {
                        best_next_push_heap.add(push_to_r_largest_ap);
                    }
                    if (collision_time_ph >= 0) {
                        best_next_push_heap.add(push_to_r_largest_ph);
                    }
                }
            }

            if (!best_next_push_heap.isEmpty()) {
                next_pushes.add(best_next_push_heap.remove());
                System.out.println("Next: " + next_pushes.peek());
                System.out.println("Expected collision time: " + next_pushes.peek().expected_collision_time);
                return;
            }

            System.err.println("Couldn't find a good move :(");
            System.err.println("Considering circularizing orbits");

            Util.Push next_push = null;

            for (int i = 1; i < indexes.length; ++i) {
                Asteroid a = asteroids[indexes[i]];
                long apoapsis_time = Util.nextAfterTime(findApoapsis(a), a, time);
                double r_ap = Util.positionAt(a, apoapsis_time).magnitude();
                double r_ph = findPeriapsisDistance(a);
                double E = reverseHohmannTransfer(a, r_ap);
                if (Math.abs(r_ap - r_ph) > a.radius()) {
                    // make it a circles!!
                    Util.Push circularize = new Util.Push(a, indexes[i], apoapsis_time, E);
                    if (next_push == null || circularize.push_time < next_push.push_time) {
                        next_push = circularize;
                    }
                }
            }

            if (next_push != null) {
                System.out.println("Adding circularization " + next_push);
                next_pushes.add(next_push);
            }
        }
    }
}
