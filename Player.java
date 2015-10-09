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
            Integer radius_indexes[] = new Integer[asteroids.length];
            for (int i = 0; i < indexes.length; ++i) {
                indexes[i] = i;
                radius_indexes[i] = i;
            }

            Arrays.sort(radius_indexes, (o1, o2) -> {
                double d = Util.positionAt(asteroids[o2], time).magnitude() -
                        Util.positionAt(asteroids[o1], time).magnitude();
                return (int) Math.signum(d);
            });

            Arrays.sort(indexes, (o1, o2) -> {
                double d = asteroids[o2].mass - asteroids[o1].mass;
                return (int) Math.signum(d);
            });

            // grab largest asteroid
            Asteroid largest = asteroids[indexes[0]];

            double r_largest_ph = largest.orbit.positionAt(findPeriapsis(largest)).magnitude();
            double r_largest_ap = largest.orbit.positionAt(findApoapsis(largest)).magnitude();
            double r_largest = (r_largest_ap + r_largest_ph) / 2;

            Util.Push potential_pushes[] = new Util.Push[indexes.length];
            for (int i = 0; i < potential_pushes.length; ++i) {
                potential_pushes[i] = null;
            }

            for (int i = 1; i < indexes.length; ++i) {
                Asteroid a = asteroids[indexes[i]];
                long push_time = time + 1;
                Util.Push next_push = moveToRadius(push_time, indexes[i], a, r_largest);

                Asteroid a2 = next_push.simulatedAsteroid(asteroids);
                long t = a2.orbit.period();
                long collision_time = push_time + t / 2;

                if (Util.findCollision(a2, largest, time, time + 10000) >= 0) {
                    if (potential_pushes[i] == null || next_push.energy < potential_pushes[i].energy) {
                        potential_pushes[i] = next_push;
                    }
                }

                if (potential_pushes[i] != null) {
                    next_pushes.add(potential_pushes[i]);
                    System.out.println("computed a next push: " + next_pushes.peek());
                    return;
                }
                break;
            }

            // no pushes computed

            for (int i = 0; i < asteroids.length; ++i) {
                double ecc = Math.sqrt(1 - Math.pow(asteroids[i].orbit.b / asteroids[i].orbit.a, 2));
//                double r_ph = largest.orbit.positionAt(findPeriapsis(largest)).magnitude();
                double r_ap = asteroids[i].orbit.positionAt(findApoapsis(asteroids[i])).magnitude();
//                double r = (r_ap + r_ph) / 2;
                if (ecc > 0.01) {
                    Util.Push reversePush = new Util.Push(asteroids[i], i, findApoapsis(asteroids[i]),
                            reverseHohmannTransfer(asteroids[i], r_ap));
                    next_pushes.add(reversePush);
                } else {
                    time_skip = largest.orbit.period();
                }
            }
        }
    }
}
