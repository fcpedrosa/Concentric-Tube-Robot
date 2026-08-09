// Hand-rolled benchmark harness for the CTR kinematics hot paths.
//
// Deliberately framework-free: warm-up + repeated timed runs, reporting
// min / median / p90. Build with the `bench` preset (release-native) and run
// on a quiet machine.

#include "TestRobot.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace
{

template <typename T> inline void doNotOptimize(const T &value)
{
    asm volatile("" : : "r,m"(value) : "memory");
}

struct Stats
{
    double min_us;
    double median_us;
    double p90_us;
};

template <typename F> Stats timeIt(F &&fn, std::size_t reps, std::size_t warmup)
{
    using clock = std::chrono::steady_clock;

    for (std::size_t i = 0; i < warmup; ++i)
        fn();

    std::vector<double> samples;
    samples.reserve(reps);
    for (std::size_t i = 0; i < reps; ++i)
    {
        const auto t0 = clock::now();
        fn();
        const auto t1 = clock::now();
        samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }

    std::sort(samples.begin(), samples.end());
    return {samples.front(), samples[samples.size() / 2UL], samples[(samples.size() * 9UL) / 10UL]};
}

void report(const std::string &name, const Stats &s)
{
    std::printf("%-28s min %12.1f us   median %12.1f us   p90 %12.1f us\n", name.c_str(), s.min_us, s.median_us,
                s.p90_us);
}

} // namespace

int main()
{
    using namespace ctr;

    std::printf("CTR kinematics benchmarks (times in microseconds)\n");
    std::printf("--------------------------------------------------\n");

    // ── FK: cold start (zero shooting guess) ─────────────────────────────────
    {
        CTR robot = testing::makeReferenceRobot();
        blaze::StaticVector<double, 6UL> q = testing::kHomeConfig;
        q[4UL] = math::deg2Rad(50.0); // non-trivial configuration: torsion active
        q[5UL] = math::deg2Rad(110.0);

        report("actuate_cold", timeIt(
                                   [&]
                                   {
                                       bvp_type guess{};
                                       const FKResult fk = robot.actuate(q, guess);
                                       doNotOptimize(fk);
                                       doNotOptimize(robot.tipPosition());
                                   },
                                   50UL, 5UL));
    }

    // ── FK: warm start (converged guess, 1 mm perturbation in beta) ─────────
    {
        CTR robot = testing::makeReferenceRobot();
        blaze::StaticVector<double, 6UL> q = testing::kHomeConfig;
        q[4UL] = math::deg2Rad(50.0);
        q[5UL] = math::deg2Rad(110.0);
        bvp_type guess{};
        std::ignore = robot.actuate(q, guess);

        blaze::StaticVector<double, 6UL> q2 = q;
        report("actuate_warm_1mm", timeIt(
                                       [&]
                                       {
                                           q2[0UL] = q[0UL] + 1.0e-3;
                                           bvp_type g = guess;
                                           const FKResult fk = robot.actuate(q2, g);
                                           doNotOptimize(fk);
                                       },
                                       50UL, 5UL));
    }

    // ── Kinematic Jacobian (3×6, finite differences) ─────────────────────────
    {
        CTR robot = testing::makeReferenceRobot();
        blaze::StaticVector<double, 6UL> q = testing::kHomeConfig;
        q[4UL] = math::deg2Rad(50.0);
        q[5UL] = math::deg2Rad(110.0);
        bvp_type guess{};
        std::ignore = robot.actuate(q, guess);
        const auto tip = robot.tipPosition();

        report("kinematic_jacobian", timeIt(
                                         [&]
                                         {
                                             const auto J = robot.kinematicJacobian(guess, tip);
                                             doNotOptimize(J);
                                         },
                                         50UL, 5UL));
    }

    // ── IK: cold solve to the demo target ────────────────────────────────────
    {
        report("ik_cold_demo_target", timeIt(
                                          [&]
                                          {
                                              CTR robot = testing::makeReferenceRobot();
                                              bvp_type guess{};
                                              std::ignore = robot.actuate(testing::kHomeConfig, guess);
                                              const blaze::StaticVector<double, 3UL> target = {-0.053210, 0.043606,
                                                                                               0.179527};
                                              const IKResult ik = robot.solveIK(target, 5.0e-4, guess);
                                              doNotOptimize(ik);
                                          },
                                          3UL, 1UL));
    }

    // ── IK: tracking 10 nearby targets (2 mm apart, warm-started) ────────────
    {
        CTR robot = testing::makeReferenceRobot();
        bvp_type guess{};
        std::ignore = robot.actuate(testing::kHomeConfig, guess);
        const blaze::StaticVector<double, 3UL> base = robot.tipPosition();

        report("ik_tracking_10x2mm", timeIt(
                                         [&]
                                         {
                                             std::size_t successes = 0UL;
                                             for (std::size_t i = 1UL; i <= 10UL; ++i)
                                             {
                                                 blaze::StaticVector<double, 3UL> target = base;
                                                 target[1UL] += 2.0e-3 * static_cast<double>(i);
                                                 const IKResult ik = robot.solveIK(target, 5.0e-4, guess);
                                                 successes += ik.converged ? 1UL : 0UL;
                                             }
                                             doNotOptimize(successes);
                                             std::printf("    (tracking successes: %zu/10)\n", successes);
                                         },
                                         1UL, 0UL));
    }

    return 0;
}
