#pragma once

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
// #include <boost/numeric/odeint/external/blaze/blaze_algebra_dispatcher.hpp>

// Define the state type
using State = blaze::StaticVector<double, 15UL>;

namespace boost
{
    namespace numeric
    {
        namespace odeint
        {
            // Custom algebra to handle pow and abs operations
            struct custom_algebra
            {
                // Implement the norm_inf method
                template<typename StateType>
                static double norm_inf(const StateType &s)
                {
                    return blaze::linfNorm(s);
                }
                
                template<typename S1, typename S2, typename Op>
                static void for_each2(S1 &s1, S2 &s2, Op op)
                {
                    for (size_t i = 0UL; i < s1.size(); ++i)
                        op(s1[i], s2[i]);
                }

                template<typename S1, typename S2, typename S3, typename Op>
                static void for_each3(S1 &s1, S2 &s2, S3 &s3, Op op)
                {
                    for (size_t i = 0UL; i < s1.size(); ++i)
                        op(s1[i], s2[i], s3[i]);
                }

                template<typename S1, typename S2, typename S3, typename S4, typename Op>
                static void for_each4(S1 &s1, S2 &s2, S3 &s3, S4 &s4, Op op)
                {
                    for (size_t i = 0UL; i < s1.size(); ++i)
                        op(s1[i], s2[i], s3[i], s4[i]);
                }
            };
        } // namespace odeint
    } // namespace numeric
} // namespace boost