#ifndef THERMO_CUBIC_SOLVER_HPP
#define THERMO_CUBIC_SOLVER_HPP

#include <vector>

namespace thermo {

/**
 * @brief Output container for cubic polynomial solutions.
 * Pure mathematical structure representing all real roots and solver diagnostics.
 */
struct CubicSolution {
    int num_real_roots{0};             ///< Total count of real roots found (1, 2, or 3)
    int num_distinct_roots{0};         ///< Count of distinct real roots after merging within tolerance
    std::vector<double> roots;         ///< Real roots sorted in ascending order (Z_0 <= Z_1 <= Z_2)
    std::vector<double> distinct_roots;///< Distinct real roots sorted in ascending order
    std::vector<double> residuals;     ///< Polynomial residual |f(Z)| for each root
    double discriminant{0.0};          ///< Depressed cubic discriminant Delta = (q/2)^2 + (p/3)^3
    double p{0.0};                     ///< Depressed cubic coefficient p
    double q{0.0};                     ///< Depressed cubic coefficient q
};

/**
 * @brief Numerically robust analytical solver for monic real cubic polynomials:
 *        P(Z) = Z^3 + c2*Z^2 + c1*Z + c0 = 0
 * 
 * Features:
 * - Depressed cubic transformation: Z = x - c2/3 -> x^3 + px + q = 0.
 * - Vieta's trigonometric method for Delta < 0 (three real roots) to prevent complex roundoff.
 * - Cardano's formulation with signed cube roots for Delta > 0 (one real root).
 * - Exact handling of near-zero discriminants (repeated roots and critical point coalescence).
 * - Horner-scheme residual verification and optional Newton refinement.
 * - Separation of mathematical solving from EOS-specific physical root selection.
 */
class CubicSolver {
public:
    /**
     * @brief Solves the monic cubic polynomial Z^3 + c2*Z^2 + c1*Z + c0 = 0.
     * @param c2 Quadratic coefficient
     * @param c1 Linear coefficient
     * @param c0 Constant term
     * @return CubicSolution structure with sorted real roots and residuals
     */
    [[nodiscard]] static CubicSolution solve(double c2, double c1, double c0);

    /**
     * @brief Evaluates the monic polynomial f(Z) = Z^3 + c2*Z^2 + c1*Z + c0 using Horner's method.
     */
    [[nodiscard]] static double evaluate_poly(double Z, double c2, double c1, double c0) noexcept;

    /**
     * @brief Evaluates the first derivative f'(Z) = 3*Z^2 + 2*c2*Z + c1 using Horner's method.
     */
    [[nodiscard]] static double evaluate_derivative(double Z, double c2, double c1) noexcept;
};

} // namespace thermo

#endif // THERMO_CUBIC_SOLVER_HPP
