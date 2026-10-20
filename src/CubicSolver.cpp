#include "thermo/CubicSolver.hpp"
#include "thermo/Constants.hpp"

#include <cmath>
#include <algorithm>

namespace thermo {

namespace {
inline constexpr double PI = 3.1415926535897932384626433832795;
inline constexpr double TWO_PI_OVER_3 = (2.0 * PI) / 3.0;
inline constexpr double FOUR_PI_OVER_3 = (4.0 * PI) / 3.0;
} // namespace

double CubicSolver::evaluate_poly(double Z, double c2, double c1, double c0) noexcept {
    // Horner's method: ((Z + c2)*Z + c1)*Z + c0
    return ((Z + c2) * Z + c1) * Z + c0;
}

double CubicSolver::evaluate_derivative(double Z, double c2, double c1) noexcept {
    // Horner's method for derivative: (3*Z + 2*c2)*Z + c1
    return (3.0 * Z + 2.0 * c2) * Z + c1;
}

CubicSolution CubicSolver::solve(double c2, double c1, double c0) {
    CubicSolution sol;

    // Depressed cubic substitution: Z = x - c2/3 -> x^3 + px + q = 0
    const double c2_sq = c2 * c2;
    const double p = c1 - (c2_sq / 3.0);
    const double q = c0 - ((c1 * c2) / 3.0) + ((2.0 * c2_sq * c2) / 27.0);

    const double p_over_3 = p / 3.0;
    const double q_over_2 = q / 2.0;

    const double delta = (q_over_2 * q_over_2) + (p_over_3 * p_over_3 * p_over_3);

    sol.discriminant = delta;
    sol.p = p;
    sol.q = q;

    std::vector<double> raw_roots;
    raw_roots.reserve(3);

    const double shift = c2 / 3.0;

    if (delta < -constants::CUBIC_DISCRIMINANT_ZERO_TOL) {
        // Case 1: Three distinct real roots (delta < 0, p < 0)
        // Use Vieta's trigonometric method to avoid complex arithmetic
        const double minus_p_over_3 = -p_over_3;
        const double sqrt_minus_p_over_3 = std::sqrt(minus_p_over_3);
        const double r_denom = minus_p_over_3 * sqrt_minus_p_over_3;

        // Clamp ratio to [-1.0, 1.0] to prevent NaN from tiny floating-point roundoff
        double ratio = -q_over_2 / r_denom;
        ratio = std::clamp(ratio, -1.0, 1.0);

        const double theta = (1.0 / 3.0) * std::acos(ratio);
        const double factor = 2.0 * sqrt_minus_p_over_3;

        raw_roots.push_back(factor * std::cos(theta) - shift);
        raw_roots.push_back(factor * std::cos(theta + TWO_PI_OVER_3) - shift);
        raw_roots.push_back(factor * std::cos(theta + FOUR_PI_OVER_3) - shift);
    }
    else if (std::abs(delta) <= constants::CUBIC_DISCRIMINANT_ZERO_TOL) {
        // Case 2: Discriminant is zero within numerical tolerance (repeated roots)
        if (std::abs(p) <= constants::CUBIC_DISCRIMINANT_ZERO_TOL &&
            std::abs(q) <= constants::CUBIC_DISCRIMINANT_ZERO_TOL) {
            // Subcase 2A: Triple root (x = 0 -> Z = -c2/3), occurs at critical point
            const double r = -shift;
            raw_roots.push_back(r);
            raw_roots.push_back(r);
            raw_roots.push_back(r);
        } else {
            // Subcase 2B: One single root and one double root
            const double u = std::cbrt(-q_over_2);
            raw_roots.push_back(2.0 * u - shift);
            raw_roots.push_back(-u - shift);
            raw_roots.push_back(-u - shift);
        }
    }
    else {
        // Case 3: Exactly one real root and two complex conjugate roots (delta > 0)
        // Use Cardano's analytical formulation with real cube roots
        const double sqrt_delta = std::sqrt(delta);
        const double u_arg = -q_over_2 + sqrt_delta;
        const double v_arg = -q_over_2 - sqrt_delta;

        const double u = std::cbrt(u_arg);
        const double v = std::cbrt(v_arg);

        raw_roots.push_back((u + v) - shift);
    }

    // Polish each root with a single Newton-Raphson iteration if needed
    for (double& r : raw_roots) {
        const double f_val = evaluate_poly(r, c2, c1, c0);
        const double f_prime = evaluate_derivative(r, c2, c1);
        if (std::abs(f_val) > 1.0e-14 && std::abs(f_prime) > 1.0e-10) {
            const double r_polished = r - (f_val / f_prime);
            if (std::abs(evaluate_poly(r_polished, c2, c1, c0)) < std::abs(f_val)) {
                r = r_polished;
            }
        }
    }

    // Sort roots in ascending order (Z_0 <= Z_1 <= Z_2)
    std::sort(raw_roots.begin(), raw_roots.end());
    sol.roots = raw_roots;
    sol.num_real_roots = static_cast<int>(sol.roots.size());

    // Compute polynomial residuals
    sol.residuals.reserve(sol.roots.size());
    for (double r : sol.roots) {
        sol.residuals.push_back(std::abs(evaluate_poly(r, c2, c1, c0)));
    }

    // Identify distinct roots within CUBIC_ROOT_MERGE_TOL
    sol.distinct_roots.reserve(sol.roots.size());
    for (double r : sol.roots) {
        if (sol.distinct_roots.empty() ||
            std::abs(r - sol.distinct_roots.back()) > constants::CUBIC_ROOT_MERGE_TOL) {
            sol.distinct_roots.push_back(r);
        }
    }
    sol.num_distinct_roots = static_cast<int>(sol.distinct_roots.size());

    return sol;
}

} // namespace thermo
