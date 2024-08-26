#pragma once

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
// #include <boost/numeric/odeint/external/blaze/blaze_algebra_dispatcher.hpp>

// Define the state type using blaze::StaticVector with 15 elements.
using State = blaze::StaticVector<double, 15UL>;

namespace boost
{
    namespace numeric
    {
        namespace odeint
        {
            /**
             * @brief Custom algebra to handle operations for blaze::StaticVector.
             *
             * This structure defines the algebra needed for integration when using Blaze vectors,
             * including operations like `pow`, `abs`, and the infinity norm.
             */
            struct custom_algebra
            {
                /**
                 * @brief Compute the infinity norm (maximum absolute value) of the state vector.
                 * 
                 * @tparam StateType Type of the state vector.
                 * @param s The state vector for which the norm is calculated.
                 * @return The infinity norm of the vector.
                 */
                template<typename StateType>
                static double norm_inf(const StateType &s)
                {
                    return blaze::linfNorm(s);  // Use Blaze's infinity norm function.
                }

                /**
                 * @brief Apply an operation element-wise to three state vectors.
                 * 
                 * This function iterates over three vectors (`s1`, `s2`, `s3`) and applies the operation `op`
                 * to their respective elements.
                 *
                 * @tparam S1 Type of the first state vector.
                 * @tparam S2 Type of the second state vector.
                 * @tparam S3 Type of the third state vector.
                 * @tparam Op Type of the operation.
                 * @param s1 First state vector.
                 * @param s2 Second state vector.
                 * @param s3 Third state vector.
                 * @param op Operation to be applied element-wise.
                 */
                template<typename S1, typename S2, typename S3, typename Op>
                static void for_each3(S1 &s1, S2 &s2, S3 &s3, Op op)
                {
                    for (size_t i = 0UL; i < s1.size(); ++i)
                        op(s1[i], s2[i], s3[i]); // Apply operation to elements at index `i`.
                }

                /**
                 * @brief Apply an operation element-wise to four state vectors.
                 * 
                 * This function iterates over four vectors and applies the operation `op`
                 * to their respective elements.
                 *
                 * @tparam S1 Type of the first state vector.
                 * @tparam S2 Type of the second state vector.
                 * @tparam S3 Type of the third state vector.
                 * @tparam S4 Type of the fourth state vector.
                 * @tparam Op Type of the operation.
                 * @param s1 First state vector.
                 * @param s2 Second state vector.
                 * @param s3 Third state vector.
                 * @param s4 Fourth state vector.
                 * @param op Operation to be applied element-wise.
                 */
                template<typename S1, typename S2, typename S3, typename S4, typename Op>
                static void for_each4(S1 &s1, S2 &s2, S3 &s3, S4 &s4, Op op)
                {
                    for (size_t i = 0UL; i < s1.size(); ++i)
                        op(s1[i], s2[i], s3[i], s4[i]); // Apply operation to elements at index `i`.
                }

                /**
                 * @brief Apply an operation element-wise to two state vectors.
                 * 
                 * This function iterates over two vectors (`s1`, `s2`) and applies the operation `op`
                 * to their respective elements.
                 *
                 * @tparam S1 Type of the first state vector.
                 * @tparam S2 Type of the second state vector.
                 * @tparam Op Type of the operation.
                 * @param s1 First state vector.
                 * @param s2 Second state vector.
                 * @param op Operation to be applied element-wise.
                 */
                template<typename S1, typename S2, typename Op>
                static void for_each2(S1 &s1, S2 &s2, Op op)
                {
                    for (size_t i = 0UL; i < s1.size(); ++i)
                        op(s1[i], s2[i]); // Apply operation to elements at index `i`.
                }
            };
        } // namespace odeint
    } // namespace numeric
} // namespace boost
