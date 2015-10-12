package pb.g3;

import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;
import pb.sim.Orbit;
import pb.sim.Point;

import java.util.Arrays;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

public class Player implements pb.sim.Player {

    private static Map<Util.Key, Long> cache_perihelion = new ConcurrentHashMap<>();

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
     * @param a asteroid to transfer (elliptical at tangent)
     * @param r_b radius to transfer it to
     * @return reverse Hohmann transfer energy
     */
    public static double reverseHohmannTransfer(Asteroid a, double r_b) {
        double r_a = a.orbit.positionAt(findPeriapsis(a)).magnitude();

        double dv = Math.sqrt(Orbit.GM / r_b) * (1 - Math.sqrt(2 * r_a / (r_a + r_b)));
        return Math.copySign(a.mass * dv * dv * 0.5, dv);
    }

    /**
     * Finds the periapsis time of a
     *
     * @param a the asteroid
     * @return the time in the period which is a's periapsis
     */
    private static long findPeriapsis(Asteroid a) {
        Util.Key k = Util.Key.factory(a);
        Point p = new Point();
        return (cache_perihelion.computeIfAbsent(k, (x) -> Util.findArgMin(0, a.orbit.period(), (y) -> {
            a.orbit.positionAt(y, p);
            return p.magnitude();
        })) + a.epoch) % a.orbit.period();
    }

    /**
     * Finds the apoapsis time of a
     * @param a the asteroid
     * @return the time in the period which is a's apoapsis
     */
    private static long findApoapsis(Asteroid a) {
        // maximum distance point is half the orbit away from the minimum distance point
        return (findPeriapsis(a) + a.orbit.period() / 2) % a.orbit.period();
    }

    /**
     * Finds the periapsis distance of a
     * @param a an asteroid
     * @return periapsis distance
     */
    private static double findPeriapsisDistance(Asteroid a) {
        return Util.positionAt(a, findPeriapsis(a)).magnitude();
    }

    /**
     * Finds the apoapsis distance of a
     * @param a an asteroid
     * @return apoapsis distance
     */
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
        double r_a = Util.positionAt(a, push_time).magnitude();
        double hohmann_energy = hohmannTransfer(a, radius); // instantaneously ellipsize it
        Util.Push hohmann_push = new Util.Push(a, aidx, push_time, hohmann_energy);
        return Util.Push.add(hohmann_push, circularize(push_time, aidx, a));
    }

    private Util.Push circularize(long push_time, int aidx, Asteroid a) {
        // circularizes at current radius
        Point p = new Point();
        Util.positionAt(a, push_time, p);
        double r = p.magnitude();
        double ortho_dir = p.direction();

        Util.velocityAt(a, push_time, p);
        double tangent_dir = p.direction();
        double current_velocity = p.magnitude();


        double bad_vel = Math.cos(ortho_dir - tangent_dir) * current_velocity;
        double good_vel = Math.sin(ortho_dir - tangent_dir) * current_velocity;

        double energy_to_remove_tangent_vel = 0.5 * a.mass * bad_vel * bad_vel;
        double reverse_bad_direction = ortho_dir + Math.PI;

        double target_v = Math.sqrt(Orbit.GM / r);

        Util.Push counteract_bad = new Util.Push(aidx, push_time, energy_to_remove_tangent_vel, reverse_bad_direction, a.mass);
        Util.Push add_good = new Util.Push(aidx, push_time, 0.5 * a.mass * Math.pow(target_v - Math.abs(good_vel), 2), ortho_dir + Math.PI / 2, a.mass);

        return Util.Push.add(counteract_bad, add_good);
    }

    public double hohmannAngularOffset(double r1, double r2) {
        return Util.normalizedAngle(Math.PI * (1 - ((1 / (2 * Math.sqrt(2))) * Math.sqrt(Math.pow(r1 / r2 + 1, 3)))));
    }

    /**
     * Evaluates whether the asteroid at a_idx will collide with target. If so, it adds it to the best_next_push_heap
     *
     * @param target
     * @param asteroids
     * @param a_idx
     * @param best_next_push_heap
     */
    public void evaluateAsteroid(Asteroid target, Asteroid[] asteroids, int a_idx, PriorityQueue<Util.Push>
            best_next_push_heap, long max_relative_start_time) {
        long push_time;
        Asteroid a = asteroids[a_idx];
        double r_peri = findPeriapsisDistance(a);
        double r_ap = findApoapsisDistance(a);
        double r_target_ph = findPeriapsisDistance(target);
        double r_target_ap = findApoapsisDistance(target);

        // this thing is basically a circle
        Point p = new Point();
        for (long dpush_time = 1; dpush_time < max_relative_start_time; ++dpush_time) {
            push_time = time + dpush_time;
            Util.positionAt(a, push_time, p);
            double local_angle = p.direction();
            double local_radius = p.magnitude();
            Util.positionAt(target, push_time, p);
            double target_angle = p.direction();

            double reqd_ang_offset = hohmannAngularOffset(local_radius, r_target_ph);
            double local_ang_offset = Util.normalizedAngle(target_angle - local_angle);

            if (Math.abs(reqd_ang_offset - local_ang_offset) > Math.PI / 72) {
                continue;
            }


            Util.Push push_to_r_largest_ph = moveToRadius(push_time, a_idx, a, r_target_ph);
            if (push_to_r_largest_ph.energy > 0) {
                Asteroid r_largest_ph_a = push_to_r_largest_ph.simulatedAsteroid(asteroids);
                push_to_r_largest_ph.expected_collision_time = Util.findCollision(target, r_largest_ph_a,
                        push_time + r_largest_ph_a.orbit.period() / 2 - 5, push_time + +r_largest_ph_a.orbit
                                .period() / 2 + 5);

                if (push_to_r_largest_ph.expected_collision_time >= 0) {
                    best_next_push_heap.add(push_to_r_largest_ph);
                }
            }

            if (Math.abs(r_target_ap - r_target_ph) < a.radius() + target.radius()) {
                // if the largest is sufficiently eccentric, try also the other intersection

                Util.Push push_to_r_largest_ap = moveToRadius(push_time, a_idx, a, r_target_ap);
                if (push_to_r_largest_ap.energy > 0) {
                    Asteroid r_largest_ap_a = push_to_r_largest_ap.simulatedAsteroid(asteroids);
                    push_to_r_largest_ap.expected_collision_time = Util.findCollision(target, r_largest_ap_a,
                            push_time, push_time + a.orbit.period() * 2);

                    if (push_to_r_largest_ap.expected_collision_time >= 0) {
                        best_next_push_heap.add(push_to_r_largest_ap);
                    }
                }
            }
        }
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids,
                     double[] energy, double[] direction) {
        ++time;

        if (num_asteroids != asteroids.length) {
            next_pushes.clear();
            Util.OrbitPair.clearCache();
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
                System.out.println("Removing " + next_push + " because it is now " + Util.toYearString(time));
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
                    time_skip = next_push.expected_collision_time + 1;
                    System.out.println("Waiting until " + Util.toYearString(time_skip));
                }

                next_pushes.remove();
                return;
            }
        } else {
        	// no pushes computed
        	System.out.println("no next push, computing more");
        	
        	long startTime = System.nanoTime();

    		Integer indexes[] = new Integer[asteroids.length];
    		Integer radius_indexes[] = new Integer[asteroids.length];
    		for (int i = 0; i < indexes.length; ++i) {
    			indexes[i] = i;
    			radius_indexes[i] = i;
            }

            Arrays.sort(radius_indexes, (o1, o2) -> (int) Math.signum(findPeriapsisDistance(asteroids[o2]) -
                    findPeriapsisDistance(asteroids[o1])));

    		Arrays.sort(indexes, (o1, o2) -> (int) Math.signum(asteroids[o2].mass - asteroids[o1].mass));

            long search_space = (long) Math.ceil((double) (time_limit - time) / (2 * asteroids.length));

            // if less than half the time limit has passed
            if ((float)time / time_limit < 0.5) {
                // grab outer asteroid, or most massive if it's much more massive
                Asteroid target = asteroids[radius_indexes[0]];
                int tgt_idx = radius_indexes[0];
                if (asteroids[indexes[0]].mass - asteroids[indexes[1]].mass > asteroids[indexes[1]].mass) {
                    target = asteroids[indexes[0]];
                    tgt_idx = indexes[0];
                }

                PriorityQueue<Util.Push> best_next_push_heap = new PriorityQueue<>((p1, p2) -> (int) Math.signum
                        (p1.energy - p2.energy));

                for (int i = 0; i < radius_indexes.length; ++i) {
                    if (tgt_idx == radius_indexes[i]) {
                        continue;
                    }
                    evaluateAsteroid(target, asteroids, radius_indexes[i], best_next_push_heap, search_space);
                }

        		if (!best_next_push_heap.isEmpty()) {
        			next_pushes.add(best_next_push_heap.remove());
        			System.out.println("Next: " + next_pushes.peek());
        			System.out.println("Expected collision time: " + Util.toYearString(next_pushes.peek()
        					.expected_collision_time));
        			System.out.println("Elapsed wall time: " + (System.nanoTime() - startTime) / 1e9);
        			return;
        		}

        		System.err.println("Couldn't find a good move :(");
        		System.err.println("Considering circularizing orbits");

        		Util.Push next_push = null;

                for (int i = 0; i < indexes.length; ++i) {
                    Asteroid a = asteroids[i];
                    long apoapsis_time = Util.nextAfterTime(findApoapsis(a), a, time);
                    double r_ap = Util.positionAt(a, apoapsis_time).magnitude();
        			double r_ph = findPeriapsisDistance(a);
        			double E = reverseHohmannTransfer(a, r_ap);
        			if (Math.abs(E) == 0) {
        				continue;
                    }
                    if (Math.abs(r_ap - r_ph) > a.radius()) {
                        // make it a circles!!
                        Util.Push circularize = circularize(apoapsis_time, i, a);
                        if (next_push == null || circularize.push_time < next_push.push_time) {
                            next_push = circularize;
        				}
        			}
        		}

        		if (next_push != null) {
        			System.out.println("Adding circularization " + next_push);
        			next_pushes.add(next_push);
        			System.out.println("Elapsed wall time: " + (System.nanoTime() - startTime) / 1e9);
        			return;
        		} else {
                    System.out.println("Skipping " + search_space + "  days");
                    time_skip = time + search_space;
                }
            }
        	// more than half the the time limit has passed
        	else {
        		Random random = new Random();
        		Asteroid a1 = null;
                int i = radius_indexes[0];
                Point v = Util.velocityAt(asteroids[i], time);
                double v1 = v.magnitude();
                double v2 = v1 * (random.nextDouble() * 0.25 + 0.05);
                double d1 = v.direction();
                double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
        		double E = 0.5 * asteroids[i].mass * v2 * v2;

                Util.Push test_push = new Util.Push(asteroids[i], i, time, E);

                try {
                    a1 = test_push.simulatedAsteroid(asteroids);
                    for (int j = 0; j < asteroids.length; ++j) {
                        if (i == j) {
                            continue;
                        }
                        long collision_time = Util.findCollision(a1, asteroids[j], time, time + 1825);
                        if (collision_time >= 0) {
                            next_pushes.add(test_push);
                        }
                    }
                } catch (InvalidOrbitException e) {
                    System.out.println("Invalid Orbit: " + e.getMessage());
                    return;
                }
            }

        }
    }
}
