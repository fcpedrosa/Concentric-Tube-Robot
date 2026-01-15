#pragma once

// #define _USE_MATH_DEFINES
#include <cmath>
#include <blaze/Math.h>
#include <iostream>

/**
 * @file mathOperations.hpp
 * @brief Math utilities for kinematic computations.
 * @ingroup ctr_math
 * @details Provides numeric helpers, rotations, and linear algebra utilities.
 */

/**
 * @brief Math helpers and utility functions.
 * @ingroup ctr_math
 */
namespace mathOp
{
	/**
	 * @brief Enumeration for defining which root finding algorithm to be used by the shooting method to solve the CTR's Boundary Value Problem (BVP)
	 */
	enum class rootFindingMethod
	{
		NEWTON_RAPHSON,
		LEVENBERG_MARQUARDT,
		POWELL_DOG_LEG,
		MODIFIED_NEWTON_RAPHSON,
		BROYDEN,
		BROYDEN_II
	};

	/**
	 * @brief Converts an angle in degrees to radians
	 *
	 * @param degree Angle in degrees.
	 * @return The corresponding angle in radians.
	 */
	inline constexpr double deg2Rad(double degree)
	{
		constexpr double pi_180 = M_PI / 180.00;
		return degree * pi_180;
	}

	/**
	 * @brief Computes the closest congruent of an angle within the [0, 3Pi] interval
	 *
	 * @param[in] angle Angle in radians.
	 * @return The corresponding closest congruent angle in radians within [0, 3Pi].
	 */
	inline double congruentAngle(double angle)
	{
		constexpr double THREE_PI = 3.00 * M_PI;
		return std::fmod(std::fabs(angle), THREE_PI) * (angle < 0.00 ? -1.00 : 1.00);
	}

	/**
	 * @brief Computes and returns a vector that is orthogonal to a vector provided by the user
	 *
	 * @param[in] v The user provided 3-dimensional static Blaze vector.
	 * @return A 3-dimensional Blaze vector orthogonal to the input vector 'v'.
	 */
	inline blaze::StaticVector<double, 3UL> orthogonal(const blaze::StaticVector<double, 3UL> &v)
	{
		blaze::StaticVector<double, 3UL> aux;
		if (std::fabs(v[0UL]) < std::fabs(v[1UL]))
		{
			aux[std::fabs(v[0UL]) < std::fabs(v[2UL]) ? 0UL : 2UL] = 1.00;
		}
		else
		{
			aux[std::fabs(v[1UL]) < std::fabs(v[2UL]) ? 1UL : 2UL] = 1.00;
		}

		return blaze::cross(v, aux);
	}

	/**
	 * @brief Computes and returns the least quaternion rotation between two user provided vectors vectors a and b. The rotation is from vector a to vector b. Both vectors are 3-dimensional.
	 *
	 * @param[in] a The input 3-dimensional static Blaze vector 'a'.
	 * @param[in] b The input 3-dimensional static Blaze vector 'b'.
	 * @return The static Blaze vector corresponding to the least quaternion rotation from 'a' to 'b'.
	 */
	inline blaze::StaticVector<double, 4UL> getRotationBetween(const blaze::StaticVector<double, 3UL> &a, const blaze::StaticVector<double, 3UL> &b)
	{
		blaze::StaticVector<double, 3UL> from = blaze::normalize(a);
		blaze::StaticVector<double, 3UL> to = blaze::normalize(b);
		blaze::StaticVector<double, 4UL> quaternion;

		// Handles the case of parallel vectors & opposing directions
		if (from == -to)
		{
			quaternion[blaze::argmax(from)] = 1.00;
		}
		else
		{
			double cos_theta = blaze::dot(from, to);
			double half_cos = std::sqrt((1.00 + cos_theta) * 0.50);
			blaze::StaticVector<double, 3UL> axis = blaze::cross(from, to) / (2.00 * half_cos);
			quaternion = {half_cos, axis[0UL], axis[1UL], axis[2UL]};
		}

		return quaternion;
	}

	/**
	 * @brief Computes pseudo-inverse of a MxN matrix, where max(M,N) = 6, via SVD decomposition.
	 *
	 * @param[in] M The input MxN Hybrid Blaze matrix M.
	 * @return The corresponding pseudo-inverse of M as a NxM Blaze matrix.
	 */
	inline blaze::HybridMatrix<double, 6UL, 6UL> pInv(const blaze::HybridMatrix<double, 6UL, 6UL> &M)
	{
		// declaring the auxiliary matrices for pInv computation
		blaze::HybridMatrix<double, 6UL, 6UL> U; // The matrix for the left singular vectors
		blaze::HybridVector<double, 6UL> s;		 // The vector for the singular values
		blaze::HybridMatrix<double, 6UL, 6UL> V; // The matrix for the right singular vectors

		// SVD decomposition of the matrix M
		try
		{
			blaze::svd(M, U, s, V);
		}
		catch (std::exception &e)
		{
			std::cerr << "Blaze SVD has failed: " << e.what() << std::endl;
			std::cerr << M << std::endl;
		}

		// computing the pseudoinverse of M via SVD decomposition
		blaze::DiagonalMatrix<blaze::HybridMatrix<double, 6UL, 6UL>> S_inv(s.size(), s.size());

		// Creating a reference to the diagonal of matrix S
		auto diag = blaze::diagonal(S_inv);
		// applies a "damping factor" to the zero singular values of the matrix M
		diag = blaze::map(s, [](double d)
						  { return (d <= 1.00E-25) ? 0.00 : d / ((d * d) + 1.00E-25); }); // Damped least squares -- SVD pseudo inverse

		return blaze::trans(U * S_inv * V);
	}

	/**
	 * @brief Converts a rotation matrix in SO(3) into the corresponding orientation in quaternions. The quaternion is passed by reference and is modified by the function.
	 *
	 * @param[in,out] h The input 4-dimensional static Blaze vector corresponding to the quaternion to be computed.
	 * @param[in] R The input static rotation Blaze matrix in SO(3).
	 */
	inline void SO3_To_Quaternion(blaze::StaticVector<double, 4UL> &h, const blaze::StaticMatrix<double, 3UL, 3UL> &R)
	{
		double n4;						   // the norm of quaternion multiplied by 4
		const double tr = blaze::trace(R); // trace of matrix

		if (tr > 0.00)
		{
			h[0UL] = tr + 1.00;
			h[1UL] = R(1UL, 2UL) - R(2UL, 1UL);
			h[2UL] = R(2UL, 0UL) - R(0UL, 2UL);
			h[3UL] = R(0UL, 1UL) - R(1UL, 0UL);
			n4 = h[0UL];
		}
		else
		{
			size_t i = 0UL;
			if (R(1UL, 1UL) > R(0UL, 0UL))
				i = 1UL;
			if (R(2UL, 2UL) > R(i, i))
				i = 2UL;

			switch (i)
			{
			case 0UL:
				h[0UL] = R(1UL, 2UL) - R(2UL, 1UL);
				h[1UL] = 1.00 + R(0UL, 0UL) - R(1UL, 1UL) - R(2UL, 2UL);
				h[2UL] = R(1UL, 0UL) + R(0UL, 1UL);
				h[3UL] = R(2UL, 0UL) + R(0UL, 2UL);
				n4 = h[1UL];
				break;
			case 1UL:
				h[0UL] = R(2UL, 0UL) - R(0UL, 2UL);
				h[1UL] = R(1UL, 0UL) + R(0UL, 1UL);
				h[2UL] = 1.00 + R(1UL, 1UL) - R(0UL, 0UL) - R(2UL, 2UL);
				h[3UL] = R(2UL, 1UL) + R(1UL, 2UL);
				n4 = h[2UL];
				break;
			case 2UL:
				h[0UL] = R(0UL, 1UL) - R(1UL, 0UL);
				h[1UL] = R(2UL, 0UL) + R(0UL, 2UL);
				h[2UL] = R(2UL, 1UL) + R(1UL, 2UL);
				h[3UL] = 1.00 + R(2UL, 2UL) - R(0UL, 0UL) - R(1UL, 1UL);
				n4 = h[3UL];
				break;
			}
		}

		h *= 1.00 / (2.00 * std::sqrt(n4));
	}

	/**
	 * @brief Converts an Euler-angle representation of 3D rotations into the corresponding orientation in quaternions.
	 * The quaternion is passed by reference and is modified by the function.
	 *
	 * @param[in] heading The heading angle in the Euler-angles representation
	 * @param[in] attitude The attitude angle in the Euler-angles representation
	 * @param[in] bank The bank angle in the Euler-angles representation
	 * @param[in,out] h The input 4-dimensional static Blaze vector corresponding to the quaternion to be computed.
	 */
	inline void euler2Quaternion(const double heading, const double attitude, const double bank, blaze::StaticVector<double, 4UL> &h)
	{
		// Gotta convert the angles to radians first.
		double theta(0.50 * heading), phi(0.50 * attitude), psi(0.50 * bank);
		double c1(cos(theta)), s1(sin(theta)), c2(cos(phi)), s2(sin(phi)), c3(cos(psi)), s3(sin(psi)), c1c2, s1s2;

		c1c2 = c1 * c2;
		s1s2 = s1 * s2;

		h[0UL] = c1c2 * c3 - s1s2 * s3;
		h[1UL] = c1c2 * s3 + s1s2 * c3;
		h[2UL] = s1 * c2 * c3 + c1 * s2 * s3;
		h[3UL] = c1 * s2 * c3 - s1 * c2 * s3;
	}

	/**
	 * @brief Converts a 3D rotation expressed in quaternions into the corresponding rotation matrix in SO(3).
	 * The rotation matrix is passed by reference and is modified by the function.
	 *
	 * @param h The input 4-dimensional static Blaze vector corresponding to the quaternion.
	 * @param R The input static rotation Blaze matrix in SO(3) to be computed.
	 */
	inline void getSO3(const blaze::StaticVector<double, 4UL> &h, blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> &R)
	{
		const double scale = 2.00 / blaze::sqrNorm(h);

		R(0UL, 0UL) = 1.00 + scale * (-h[2UL] * h[2UL] - h[3UL] * h[3UL]);
		R(0UL, 1UL) = scale * (h[1UL] * h[2UL] - h[3UL] * h[0UL]);
		R(0UL, 2UL) = scale * (h[1UL] * h[3UL] + h[2UL] * h[0UL]);

		R(1UL, 0UL) = scale * (h[1UL] * h[2UL] + h[3UL] * h[0UL]);
		R(1UL, 1UL) = 1.00 + scale * (-h[1UL] * h[1UL] - h[3UL] * h[3UL]);
		R(1UL, 2UL) = scale * (h[2UL] * h[3UL] - h[1UL] * h[0UL]);

		R(2UL, 0UL) = scale * (h[1UL] * h[3UL] - h[2UL] * h[0UL]);
		R(2UL, 1UL) = scale * (h[2UL] * h[3UL] + h[1UL] * h[0UL]);
		R(2UL, 2UL) = 1.00 + scale * (-h[1UL] * h[1UL] - h[2UL] * h[2UL]);
	}

	/**
	 * @brief Computes and returns the rotation matrix in SO(3) around the z-axis R_z for an angle 'theta' in radians.
	 *
	 * @param theta The input angle in radians.
	 * @return The corresponding rotation Blaze matrix aroudn the z-axis.
	 */
	inline blaze::StaticMatrix<double, 3UL, 3UL> rotz(const double &theta)
	{
		double c(cos(theta)), s(sin(theta));
		blaze::StaticMatrix<double, 3UL, 3UL> R{{c, -s, 0.00},
												{s, c, 0.00},
												{0.00, 0.00, 1.00}};

		return R;
	}

	/**
	 * @brief Computes and returns the time derivative of the rotation matrix in SO(3) around the z-axis R_z for an angle 'theta' in radians.
	 *
	 * @param theta The input angle in radians.
	 * @return The corresponding time derivative of the rotation Blaze matrix aroudn the z-axis.
	 */
	inline blaze::StaticMatrix<double, 3UL, 3UL> rotz_dot_transpose(const double &theta)
	{
		double c(cos(theta)), s(sin(theta));
		blaze::StaticMatrix<double, 3UL, 3UL> dR{{-s, c, 0.00},
												 {-c, -s, 0.00},
												 {0.00, 0.00, 0.00}};

		return dR;
	}

	/**
	 * @brief Computes and returns hat operator applied on a 3-dimensional vector. This operation maps a 3-dimensional vector in R^3 into the corresponding element in the space of skew-symmetric matrices in so(3).
	 *
	 * @param v The input 3-dimensional static Blaze vector.
	 * @return The corresponding skew-symmetric matrix in so(3).
	 */
	inline blaze::StaticMatrix<double, 3UL, 3UL> hatOperator(const blaze::StaticVector<double, 3UL> &v)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> Res{{0.00, -v[2UL], v[1UL]},
												  {v[2UL], 0.00, -v[0UL]},
												  {-v[1UL], v[0UL], 0.00}};

		return Res;
	}

	/**
	 * @brief Computes and returns the pre-multiplication of a 3-dimensional square matrix by the hat operator applied on a 3-dimensional vector.
	 * That is to say, the operation: \f$\hat{v} M\f$.
	 *
	 * @param[in] v The input 3-dimensional static Blaze vector.
	 * @param[in] M The corresponding 3-dimensional square Blaze matrix to be premultiplied.
	 * @return The 3-dimensional square Blaze matrix containing the result of the operation.
	 */
	inline blaze::StaticMatrix<double, 3UL, 3UL> hatPreMultiply(const blaze::StaticVector<double, 3UL> &v, const blaze::StaticMatrix<double, 3UL, 3UL> &M)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> Res{{-M(1UL, 0UL) * v[2UL] + M(2UL, 0UL) * v[1UL], -M(1UL, 1UL) * v[2UL] + M(2UL, 1UL) * v[1UL], -M(1UL, 2UL) * v[2UL] + M(2UL, 2UL) * v[1UL]},
												  {M(0UL, 0UL) * v[2UL] - M(2UL, 0UL) * v[0UL], M(0UL, 1UL) * v[2UL] - M(2UL, 1UL) * v[0UL], M(0UL, 2UL) * v[2UL] - M(2UL, 2UL) * v[0UL]},
												  {-M(0UL, 0UL) * v[1UL] + M(1UL, 0UL) * v[0UL], -M(0UL, 1UL) * v[1UL] + M(1UL, 1UL) * v[0UL], -M(0UL, 2UL) * v[1UL] + M(1UL, 2UL) * v[0UL]}};

		return Res;
	}

	/**
	 * @brief Computes and returns the pos-multiplication of a 3-dimensional square matrix by the hat operator applied on a 3-dimensional vector.
	 * That is to say, the operation: \f$M \hat{v}\f$.
	 *
	 * @param[in] v The input 3-dimensional static Blaze vector.
	 * @param[in] M The corresponding 3-dimensional square Blaze matrix to be premultiplied.
	 * @return The 3-dimensional square Blaze matrix containing the result of the operation.
	 */
	inline blaze::StaticMatrix<double, 3UL, 3UL> hatPostMultiply(const blaze::StaticMatrix<double, 3UL, 3UL> &M, const blaze::StaticVector<double, 3UL> &v)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> Res{{M(0UL, 1UL) * v[2UL] - M(0UL, 2UL) * v[1UL], -M(0UL, 0UL) * v[2UL] + M(0UL, 2UL) * v[0UL], M(0UL, 0UL) * v[1UL] - M(0UL, 1UL) * v[0UL]},
												  {M(1UL, 1UL) * v[2UL] - M(1UL, 2UL) * v[1UL], -M(1UL, 0UL) * v[2UL] + M(1UL, 2UL) * v[0UL], M(1UL, 0UL) * v[1UL] - M(1UL, 1UL) * v[0UL]},
												  {M(2UL, 1UL) * v[2UL] - M(2UL, 2UL) * v[1UL], -M(2UL, 0UL) * v[2UL] + M(2UL, 2UL) * v[0UL], M(2UL, 0UL) * v[1UL] - M(2UL, 1UL) * v[0UL]}};

		return Res;
	}

	/**
	 * @brief Computes and returns product between a vector and the transpose of a matrix.
	 * That is to say, the operation: \f$R^\top v\f$.
	 *
	 * @param[in] R The input 3-dimensional static square Blaze matrix.
	 * @param[in] v The corresponding 3-dimensional square Blaze vector.
	 * @return The 3-dimensional square Blaze vector containing the result of the operation.
	 */
	inline blaze::StaticVector<double, 3UL> transposePreMultiply(const blaze::StaticMatrix<double, 3UL, 3UL> &R, const blaze::StaticVector<double, 3UL> &v)
	{
		blaze::StaticVector<double, 3UL> res{R(0UL, 0UL) * v[0UL] + R(1UL, 0UL) * v[1UL] + R(2UL, 0UL) * v[2UL],
											 R(0UL, 1UL) * v[0UL] + R(1UL, 1UL) * v[1UL] + R(2UL, 1UL) * v[2UL],
											 R(0UL, 2UL) * v[0UL] + R(1UL, 2UL) * v[1UL] + R(2UL, 2UL) * v[2UL]};

		return res;
	}

	/**
	 * @brief Computes and returns the derivative of the quaternion orientation with respect to arc-length
	 *
	 * @param[in] u The input 3-dimensional static Blaze vector corresponding to the local curvature of the backbone.
	 * @param[in] h The input 4-dimensional square Blaze vector corresponding to the orientation of the local body frame of the backbone.
	 * @return The 4-dimensional static Blaze vector corresponding to the local rate of change of the orientation of the backbone with respect to arc-length.
	 */
	inline blaze::StaticVector<double, 4UL> quaternionDiff(const blaze::StaticVector<double, 3UL> &u, const blaze::StaticVector<double, 4UL> &h)
	{
		blaze::StaticVector<double, 4UL> hs{0.50 * (-u[0UL] * h[1UL] - u[1UL] * h[2UL] - u[2UL] * h[3UL]),
											0.50 * (u[0UL] * h[0UL] + u[2UL] * h[2UL] - u[1UL] * h[3UL]),
											0.50 * (u[1UL] * h[0UL] - u[2UL] * h[1UL] + u[0UL] * h[3UL]),
											0.50 * (u[2UL] * h[0UL] + u[1UL] * h[1UL] - u[0UL] * h[2UL])};

		return hs;
	}
}