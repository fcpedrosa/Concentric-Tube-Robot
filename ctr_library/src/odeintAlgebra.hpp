#pragma once

// PRIVATE header — Boost.Odeint glue for Blaze state vectors. Not installed;
// scheduled for removal when the integrator moves to a hand-driven
// fixed-step loop.

#include <boost/numeric/odeint.hpp>
#include "ctr/Types.hpp"

namespace ctr
{

/**
 * @brief Custom Boost.Odeint algebra that delegates element-wise operations
 *        and norm computation to Blaze.
 *
 * Defined in namespace ctr (not boost::numeric::odeint) to avoid injecting
 * user-defined types into a third-party namespace, which is undefined behaviour.
 * Pass it as an explicit template argument to the stepper type.
 */
struct BlazeBVPAlgebra
{
    template <typename StateType> static double norm_inf(const StateType &s) { return blaze::linfNorm(s); }

    template <typename S1, typename S2, typename Op> static void for_each2(S1 &s1, S2 &s2, Op op)
    {
        for (std::size_t i = 0UL; i < s1.size(); ++i)
            op(s1[i], s2[i]);
    }

    template <typename S1, typename S2, typename S3, typename Op> static void for_each3(S1 &s1, S2 &s2, S3 &s3, Op op)
    {
        for (std::size_t i = 0UL; i < s1.size(); ++i)
            op(s1[i], s2[i], s3[i]);
    }

    template <typename S1, typename S2, typename S3, typename S4, typename Op>
    static void for_each4(S1 &s1, S2 &s2, S3 &s3, S4 &s4, Op op)
    {
        for (std::size_t i = 0UL; i < s1.size(); ++i)
            op(s1[i], s2[i], s3[i], s4[i]);
    }
};

} // namespace ctr
