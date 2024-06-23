#pragma once

// #define _USE_MATH_DEFINES
#include <cmath>
#include <blaze/Math.h>
#include <iostream>

namespace mathOp
{
	// enumeration for defining which root finding method to be used
	enum class rootFindingMethod
	{
		NEWTON_RAPHSON,
		LEVENBERG_MARQUARDT,
		POWELL_DOG_LEG,
		MODIFIED_NEWTON_RAPHSON,
		BROYDEN,
		BROYDEN_II
	};

	// Function for converting degrees to radians
	inline constexpr double deg2Rad(double degree)
	{
		constexpr double pi_180 = M_PI / 180.00;
		return degree * pi_180;
	}

	// Function to reduce any angle to its congruent in the [0, 2*Pi) interval
	inline double wrapTo2Pi(double theta)
	{
		constexpr double TWO_PI = 2.00 * M_PI;
		return std::fmod(std::fmod(theta, TWO_PI) + TWO_PI, TWO_PI);
	}

	// Function that computes the congruent angle of theta in [0, 3Pi]
	inline double congruentAngle(double angle)
	{
		constexpr double THREE_PI = 3.00 * M_PI;
		return std::fmod(std::fabs(angle), THREE_PI) * (angle < 0.00 ? -1.00 : 1.00);
	}

	// Function that returns an orthogonal vector to v
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

	// Function that returns the least quaternion rotation between vectors a and b
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

	// method for computing the Penrose Pseudoinverse via SVD decomposition
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

	// function that returns a rotation matrix in SO(3) from a set of non-unity quaternions
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

	// function that returns the rotation matrix Rz of any angle theta
	inline blaze::StaticMatrix<double, 3UL, 3UL> rotz(const double &theta)
	{
		double c(cos(theta)), s(sin(theta));
		blaze::StaticMatrix<double, 3UL, 3UL> R{{c, -s, 0.00},
												{s, c, 0.00},
												{0.00, 0.00, 1.00}};

		return R;
	}

	// function that returns the derivative of the rotation matrix Rz
	inline blaze::StaticMatrix<double, 3UL, 3UL> rotz_dot_transpose(const double &theta)
	{
		double c(cos(theta)), s(sin(theta));
		blaze::StaticMatrix<double, 3UL, 3UL> dR{{-s, c, 0.00},
												 {-c, -s, 0.00},
												 {0.00, 0.00, 0.00}};

		return dR;
	}

	// function that computes the premultiplication of a matrix M by v^, i.e., v^ * M
	inline blaze::StaticMatrix<double, 3UL, 3UL> hatOperator(const blaze::StaticVector<double, 3UL> &v)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> Res{{0.00, -v[2UL], v[1UL]},
												  {v[2UL], 0.00, -v[0UL]},
												  {-v[1UL], v[0UL], 0.00}};

		return Res;
	}

	// function that computes the premultiplication of a matrix M by v^, i.e., v^ * M
	inline blaze::StaticMatrix<double, 3UL, 3UL> hatPreMultiply(const blaze::StaticVector<double, 3UL> &v, const blaze::StaticMatrix<double, 3UL, 3UL> &M)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> Res{{-M(1UL, 0UL) * v[2UL] + M(2UL, 0UL) * v[1UL], -M(1UL, 1UL) * v[2UL] + M(2UL, 1UL) * v[1UL], -M(1UL, 2UL) * v[2UL] + M(2UL, 2UL) * v[1UL]},
												  {M(0UL, 0UL) * v[2UL] - M(2UL, 0UL) * v[0UL], M(0UL, 1UL) * v[2UL] - M(2UL, 1UL) * v[0UL], M(0UL, 2UL) * v[2UL] - M(2UL, 2UL) * v[0UL]},
												  {-M(0UL, 0UL) * v[1UL] + M(1UL, 0UL) * v[0UL], -M(0UL, 1UL) * v[1UL] + M(1UL, 1UL) * v[0UL], -M(0UL, 2UL) * v[1UL] + M(1UL, 2UL) * v[0UL]}};

		return Res;
	}

	// function that computes the posmultiplication of a matrix M by v^, i.e., M * v^
	inline blaze::StaticMatrix<double, 3UL, 3UL> hatPostMultiply(const blaze::StaticMatrix<double, 3UL, 3UL> &M, const blaze::StaticVector<double, 3UL> &v)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> Res{{M(0UL, 1UL) * v[2UL] - M(0UL, 2UL) * v[1UL], -M(0UL, 0UL) * v[2UL] + M(0UL, 2UL) * v[0UL], M(0UL, 0UL) * v[1UL] - M(0UL, 1UL) * v[0UL]},
												  {M(1UL, 1UL) * v[2UL] - M(1UL, 2UL) * v[1UL], -M(1UL, 0UL) * v[2UL] + M(1UL, 2UL) * v[0UL], M(1UL, 0UL) * v[1UL] - M(1UL, 1UL) * v[0UL]},
												  {M(2UL, 1UL) * v[2UL] - M(2UL, 2UL) * v[1UL], -M(2UL, 0UL) * v[2UL] + M(2UL, 2UL) * v[0UL], M(2UL, 0UL) * v[1UL] - M(2UL, 1UL) * v[0UL]}};

		return Res;
	}

	// efficiently computes the product between a vector and the transpose of a matrix (R^T * v)
	inline blaze::StaticVector<double, 3UL> transposePreMultiply(const blaze::StaticMatrix<double, 3UL, 3UL> &R, const blaze::StaticVector<double, 3UL> &v)
	{
		blaze::StaticVector<double, 3UL> res{R(0UL, 0UL) * v[0UL] + R(1UL, 0UL) * v[1UL] + R(2UL, 0UL) * v[2UL],
											 R(0UL, 1UL) * v[0UL] + R(1UL, 1UL) * v[1UL] + R(2UL, 1UL) * v[2UL],
											 R(0UL, 2UL) * v[0UL] + R(1UL, 2UL) * v[1UL] + R(2UL, 2UL) * v[2UL]};

		return res;
	}

	// function that computes the differential quaternion evolution
	inline blaze::StaticVector<double, 4UL> quaternionDiff(const blaze::StaticVector<double, 3UL> &u, const blaze::StaticVector<double, 4UL> &h)
	{
		blaze::StaticVector<double, 4UL> hs{0.50 * (-u[0UL] * h[1UL] - u[1UL] * h[2UL] - u[2UL] * h[3UL]),
											0.50 * (u[0UL] * h[0UL] + u[2UL] * h[2UL] - u[1UL] * h[3UL]),
											0.50 * (u[1UL] * h[0UL] - u[2UL] * h[1UL] + u[0UL] * h[3UL]),
											0.50 * (u[2UL] * h[0UL] + u[1UL] * h[1UL] - u[0UL] * h[2UL])};

		return hs;
	}
}