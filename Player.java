package pb.g3;

import pb.sim.Asteroid;
import pb.sim.Orbit;
import pb.sim.Point;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

public class Player implements pb.sim.Player {

    private static Map<Util.Key, Long> cache_perihelion = new ConcurrentHashMap<>();

    // current time, time limit
    private long time = -1;
    private long time_limit = -1;
    private long time_skip = -1;
    private int num_asteroids = -1;
    private long num_pushes = 0;
    private PriorityQueue<Util.Push> next_pushes = new PriorityQueue<>((p1, p2) -> Long.compare(p1.push_time, p2
            .push_time));

    private Asteroid target;
    private double target_mass;
    private boolean consider_circularize = true;

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
     * @param a   asteroid to transfer (elliptical at tangent)
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
     *
     * @param a the asteroid
     * @return the time in the period which is a's apoapsis
     */
    private static long findApoapsis(Asteroid a) {
        // maximum distance point is half the orbit away from the minimum distance point
        return (findPeriapsis(a) + a.orbit.period() / 2) % a.orbit.period();
    }

    /**
     * Finds the periapsis distance of a
     *
     * @param a an asteroid
     * @return periapsis distance
     */
    private static double findPeriapsisDistance(Asteroid a) {
        return Util.positionAt(a, findPeriapsis(a)).magnitude();
    }

    /**
     * Finds the apoapsis distance of a
     *
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

    private Util.Push moveToRadius(long push_time, Asteroid a, double radius) {
        double r_a = Util.positionAt(a, push_time).magnitude();
        double hohmann_energy = hohmannTransfer(a, radius); // instantaneously ellipsize it
        Util.Push hohmann_push = new Util.Push(a, push_time, hohmann_energy);
        return Util.Push.add(hohmann_push, circularize(push_time, a));
    }

    private Util.Push circularize(long push_time, Asteroid a) {
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

        Util.Push counteract_bad = new Util.Push(a, push_time, energy_to_remove_tangent_vel, reverse_bad_direction, a
                .mass);
        Util.Push add_good = new Util.Push(a, push_time, 0.5 * a.mass * Math.pow(target_v - Math.abs(good_vel), 2),
                ortho_dir + Math.PI / 2, a.mass);

        return Util.Push.add(counteract_bad, add_good);
    }

    public double hohmannAngularOffset(double r1, double r2) {
        return Util.normalizedAngle(Math.PI * (1 - ((1 / (2 * Math.sqrt(2))) * Math.sqrt(Math.pow(r1 / r2 + 1, 3)))));
    }

    private long findCollisionGivenHohmannTransfer(Asteroid target, Util.Push p) {
        Asteroid pushed = p.simulatedAsteroid();
        long t = pushed.orbit.period() / 2;
        long t_col = p.push_time + t;
        Point p1 = Util.positionAt(target, t_col);
        Point p2 = Util.positionAt(pushed, t_col);
        if (Point.distance(p1, p2) < target.radius() + pushed.radius()) {
            return t_col;
        } else {
            return -1;
        }
    }

    /**
     * Evaluates whether the asteroid at a_idx will collide with target. If so, it adds it to the best_next_push_heap
     *
     * @param target              target asteroid
     * @param a                   asteroid to attempt to collide
     * @param best_next_push_heap heap to keep track of the best next push
     */
    public void evaluateAsteroid(Asteroid target, Asteroid a, PriorityQueue<Util.Push>
            best_next_push_heap, long time, long max_time) {
        double r_target_ph = findPeriapsisDistance(target);
        double r_target_ap = findApoapsisDistance(target);

        Point p = new Point();

        long must_collide_by = Math.min(max_time + (max_time - time) / 2, time_limit - 1);

        for (long push_time = time; push_time < max_time; ++push_time) {
            Util.positionAt(a, push_time, p);
            double local_angle = p.direction();
            double local_radius = p.magnitude();
            Util.positionAt(target, push_time, p);
            double target_angle = p.direction();

            double reqd_ang_offset = hohmannAngularOffset(local_radius, r_target_ph);
            double local_ang_offset = Util.normalizedAngle(target_angle - local_angle);

            if (Math.abs(Util.normalizedAngle(reqd_ang_offset - local_ang_offset)) > Math.PI / 72) {
                continue;
            }

            Util.Push push_target_ph = moveToRadius(push_time, a, r_target_ph);
            if (push_target_ph.energy > 0) {
                push_target_ph.expected_collision_time = findCollisionGivenHohmannTransfer(target, push_target_ph);

                if (push_target_ph.expected_collision_time >= 0 && push_target_ph.expected_collision_time <=
                        must_collide_by) {
                    best_next_push_heap.add(push_target_ph);
                }
            }

            Util.Push push_target_ap = moveToRadius(push_time, a, r_target_ap);
            if (push_target_ap.energy > 0) {
                push_target_ap.expected_collision_time = findCollisionGivenHohmannTransfer(target, push_target_ap);

                if (push_target_ap.expected_collision_time >= 0 && push_target_ap.expected_collision_time <= must_collide_by) {
                    best_next_push_heap.add(push_target_ap);
                }
            }
        }
    }

    public Util.Push circularizeAsteroids(Asteroid asteroids[], long time, long max_time) {
        Util.Push next_push = null;

        for (Asteroid a : asteroids) {
            long apoapsis_time = Util.nextAfterTime(findApoapsis(a), a, time);
            if (apoapsis_time > max_time) {
                continue;
            }
            double r_ap = Util.positionAt(a, apoapsis_time).magnitude();
            double r_ph = findPeriapsisDistance(a);
            double E = reverseHohmannTransfer(a, r_ap);
            if (Math.abs(E) == 0) {
                continue;
            }
            if (Math.abs(r_ap - r_ph) > a.radius()) {
                // make it a circles!!
                Util.Push circularize = circularize(apoapsis_time, a);
                if (next_push == null || circularize.energy < next_push.energy) {
                    next_push = circularize;
                }
            }
        }
        return next_push;
    }

    public void computePushes(Asteroid asteroids[], long localtime, long max_time) {

        // sorted ascending
        PriorityQueue<Util.Push> best_next_push_heap = new PriorityQueue<>((p1, p2) -> Double.compare(p1.energy, p2
                .energy));

        if (target == null) {
            Asteroid asteroids_by_radius[] = new Asteroid[asteroids.length];
            System.arraycopy(asteroids, 0, asteroids_by_radius, 0, asteroids.length);
            // sorted descending
            Arrays.sort(asteroids_by_radius, (a1, a2) -> Double.compare(findPeriapsisDistance(a2),
                    findPeriapsisDistance(a1)));

            double costs[] = new double[asteroids.length];

            for (int i = 0; i < asteroids_by_radius.length; ++i) {
                double r = findPeriapsisDistance(asteroids_by_radius[i]);
                for (Asteroid a : asteroids_by_radius) {
                    costs[i] += Math.abs(moveToRadius(0, a, r).energy);
                }
            }

            target = asteroids_by_radius[Util.findArgMinI(0, costs.length, (x) -> costs[x])];
            target_mass = target.mass;

            if (Math.abs(findPeriapsisDistance(asteroids_by_radius[0]) - findPeriapsisDistance
                    (asteroids_by_radius[asteroids_by_radius.length - 1])) < asteroids_by_radius[0].radius()) {
                System.out.println("Same-orbit detected, initiate random push");
                consider_circularize = false;
                for (int i = 0; i < asteroids.length; ++i) {
                    next_pushes.add(moveToRadius(time + i, asteroids_by_radius[i], findPeriapsisDistance
                            (asteroids_by_radius[i]) * (i * 0.5) + 1));
                }
                return;
            }
        }

        for (Asteroid a : asteroids) {
            if (a.id == target.id) {
                continue;
            }
            evaluateAsteroid(target, a, best_next_push_heap, localtime, max_time);
        }

        if (!best_next_push_heap.isEmpty()) {
            next_pushes.add(best_next_push_heap.remove());
            System.out.println("Next: " + next_pushes.peek());
            System.out.println("Expected collision time: " + Util.toYearString(next_pushes.peek()
                    .expected_collision_time));
        } else {
            System.err.println("Couldn't find a good move :(");
            if (consider_circularize) {
                System.err.println("Considering circularizing orbits");

                Util.Push next_push = circularizeAsteroids(asteroids, localtime + 1, max_time);

                if (next_push != null) {
                    System.out.println("Adding circularization " + next_push);
                    next_pushes.add(next_push);
                } else {
                    System.out.println("Skipping until " + Util.toYearString(max_time));
                    time_skip = max_time;
                }
            }
        }
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids,
                     double[] energy, double[] direction) {
        ++time;


        Asteroid asteroids_by_mass[] = new Asteroid[asteroids.length];
        System.arraycopy(asteroids, 0, asteroids_by_mass, 0, asteroids.length);
        // sorted ascending
        Arrays.sort(asteroids_by_mass, (a1, a2) -> Double.compare(a1.mass, a2.mass));

        Set<Long> extant_asteroids = new HashSet<>();
        double max_mass = -1;
        double total_mass = 0;
        for (Asteroid a : asteroids) {
            extant_asteroids.add(a.id);
            max_mass = Math.max(max_mass, a.mass);
            total_mass += a.mass;
        }

        double frac_mass = max_mass / total_mass;
        double remaining_frac_mass = 0.5 - frac_mass;
        int max_asteroids_needed = 0;
        for (Asteroid a : asteroids_by_mass) {
            ++max_asteroids_needed;
            remaining_frac_mass -= a.mass / total_mass;
            if (remaining_frac_mass < 0) {
                break;
            }
        }

        long search_space = (long) (((time_limit - time) * 1.0) / max_asteroids_needed);


        if (num_asteroids != asteroids.length) {
            next_pushes.clear();
            Util.OrbitPair.clearCache();
            System.out.println(String.format("#asteroids changed from %d to %d, dropping queued pushes", num_asteroids,
                    asteroids.length));
            num_asteroids = asteroids.length;

            // finding new target
            if (!extant_asteroids.contains(target.id)) {
                for (Asteroid a : asteroids) {
                    if (a.mass == target_mass) {
                        target = a;
                    }
                }
            }
        }

        if (time_skip > time) {
            return;
        } else {
            time_skip = -1;
        }

        if (!next_pushes.isEmpty()) {
            Util.Push next_push = next_pushes.peek();
            while (time > next_push.push_time || !extant_asteroids.contains(next_push.asteroid.id)) {
                System.out.println("Removing " + next_push + " because it is now invalid");
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
                if (extant_asteroids.contains(next_push.asteroid.id)) {
                    // apply push
                    for (int i = 0; i < asteroids.length; ++i) {
                        if (asteroids[i].id == next_push.asteroid.id) {
                            energy[i] = next_push.energy;
                            direction[i] = next_push.direction;
                            ++num_pushes;
                        }
                    }
                    System.out.println("Making push " + next_push);

                    if (time_skip < next_push.expected_collision_time) {
                        time_skip = next_push.expected_collision_time + 1;
                        target_mass = target.mass + next_push.asteroid.mass;
                        System.out.println("Waiting until " + Util.toYearString(time_skip));
                    }

                    next_pushes.remove();
                } else {
                    System.err.println("Could not find asteroid for scheduled push, dropping push");
                    next_pushes.remove();
                }
            }
        } else {
            // no pushes computed
            System.out.println("no next push, computing more");

            long startTime = System.nanoTime();

            computePushes(asteroids, time, time + search_space);

            System.out.println("Elapsed wall time: " + (System.nanoTime() - startTime) / 1e9);
        }
    }
}
