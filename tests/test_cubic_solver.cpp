#include "thermo/CubicSolver.hpp"
#include "thermo/Constants.hpp"

#include <iostream>
#include <cassert>
#include <cmath>
#include <string>

namespace {

int g_pass_count = 0;
int g_fail_count = 0;

void report_test(const std::string& test_name, bool condition) {
    if (condition) {
        std::cout << "[PASS] " << test_name << "\n";
        ++g_pass_count;
    } else {
        std::cerr << "[FAIL] " << test_name << "\n";
        ++g_fail_count;
    }
}

bool approx_equal(double a, double b, double tol = 1.0e-6) {
    return std::abs(a - b) <= tol;
}

} // namespace

int main() {
    std::cout << "==================================================\n";
    std::cout << "RUNNING CUBIC SOLVER UNIT TESTS (Session 2)\n";
    std::cout << "==================================================\n";

    // ------------------------------------------------------------------------
    // Test 1: Three distinct integer roots: (Z - 1)(Z - 2)(Z - 3) = 0
    // Expansion: Z^3 - 6*Z^2 + 11*Z - 6 = 0
    // Expected roots: 1.0, 2.0, 3.0
    // ------------------------------------------------------------------------
    {
        const auto sol = thermo::CubicSolver::solve(-6.0, 11.0, -6.0);
        report_test("Test 1: (Z-1)(Z-2)(Z-3) has 3 real roots", sol.num_real_roots == 3);
        report_test("Test 1: 3 distinct roots", sol.num_distinct_roots == 3);
        report_test("Test 1: Root 1 == 1.0", approx_equal(sol.roots[0], 1.0));
        report_test("Test 1: Root 2 == 2.0", approx_equal(sol.roots[1], 2.0));
        report_test("Test 1: Root 3 == 3.0", approx_equal(sol.roots[2], 3.0));
        for (size_t i = 0; i < sol.residuals.size(); ++i) {
            report_test("Test 1: Residual " + std::to_string(i) + " <= 1e-12",
                        sol.residuals[i] <= thermo::constants::CUBIC_POLY_RESIDUAL_TOL);
        }
    }

    // ------------------------------------------------------------------------
    // Test 2: Single real root: Z^3 + Z - 2 = 0 -> (Z - 1)(Z^2 + Z + 2) = 0
    // Expected root: 1.0 (complex roots: -0.5 +/- 1.3228i)
    // ------------------------------------------------------------------------
    {
        const auto sol = thermo::CubicSolver::solve(0.0, 1.0, -2.0);
        report_test("Test 2: Z^3 + Z - 2 has 1 real root", sol.num_real_roots == 1);
        report_test("Test 2: Distinct roots count == 1", sol.num_distinct_roots == 1);
        report_test("Test 2: Root == 1.0", approx_equal(sol.roots[0], 1.0));
        report_test("Test 2: Discriminant > 0", sol.discriminant > 0.0);
        report_test("Test 2: Residual <= 1e-12",
                    sol.residuals[0] <= thermo::constants::CUBIC_POLY_RESIDUAL_TOL);
    }

    // ------------------------------------------------------------------------
    // Test 3: Repeated double root: (Z - 1)^2 * (Z - 4) = 0
    // Expansion: (Z^2 - 2Z + 1)(Z - 4) = Z^3 - 6*Z^2 + 9*Z - 4 = 0
    // Expected roots: 1.0, 1.0, 4.0 (distinct: 1.0, 4.0)
    // ------------------------------------------------------------------------
    {
        const auto sol = thermo::CubicSolver::solve(-6.0, 9.0, -4.0);
        report_test("Test 3: (Z-1)^2*(Z-4) has 3 roots", sol.num_real_roots == 3);
        report_test("Test 3: 2 distinct roots", sol.num_distinct_roots == 2);
        report_test("Test 3: Root 0 == 1.0", approx_equal(sol.roots[0], 1.0));
        report_test("Test 3: Root 1 == 1.0", approx_equal(sol.roots[1], 1.0));
        report_test("Test 3: Root 2 == 4.0", approx_equal(sol.roots[2], 4.0));
        report_test("Test 3: Distinct[0] == 1.0", approx_equal(sol.distinct_roots[0], 1.0));
        report_test("Test 3: Distinct[1] == 4.0", approx_equal(sol.distinct_roots[1], 4.0));
    }

    // ------------------------------------------------------------------------
    // Test 4: Triple root: (Z - 2)^3 = 0
    // Expansion: Z^3 - 6*Z^2 + 12*Z - 8 = 0
    // Expected root: 2.0 (multiplicity 3)
    // ------------------------------------------------------------------------
    {
        const auto sol = thermo::CubicSolver::solve(-6.0, 12.0, -8.0);
        report_test("Test 4: (Z-2)^3 has 3 roots", sol.num_real_roots == 3);
        report_test("Test 4: 1 distinct root", sol.num_distinct_roots == 1);
        report_test("Test 4: Distinct root == 2.0", approx_equal(sol.distinct_roots[0], 2.0));
        report_test("Test 4: Residual <= 1e-12",
                    sol.residuals[0] <= thermo::constants::CUBIC_POLY_RESIDUAL_TOL);
    }

    // ------------------------------------------------------------------------
    // Test 5: Audited PR-1976 Case 1: Methane at 300 K, 5 MPa
    // Dimensionless parameters: A = 0.16248591, B = 0.05372479
    // Expected root: Z = 0.901844
    // ------------------------------------------------------------------------
    {
        const double A = 0.16248591;
        const double B = 0.05372479;
        const double c2 = -(1.0 - B);
        const double c1 = A - 3.0 * B * B - 2.0 * B;
        const double c0 = -(A * B - B * B - B * B * B);
        const auto sol = thermo::CubicSolver::solve(c2, c1, c0);
        report_test("Test 5 (PR Methane 5MPa): Exactly 1 real root", sol.num_real_roots == 1);
        report_test("Test 5 (PR Methane 5MPa): Root Z == 0.901844", approx_equal(sol.roots[0], 0.90184394, 1.0e-5));
        report_test("Test 5 (PR Methane 5MPa): Residual <= 1e-12",
                    sol.residuals[0] <= thermo::constants::CUBIC_POLY_RESIDUAL_TOL);
    }

    // ------------------------------------------------------------------------
    // Test 6: Audited PR-1976 Case 2: CO2 at 300 K, 4 MPa
    // Dimensionless parameters: A = 0.25725664, B = 0.04276376
    // Expected root: Z = 0.752334
    // ------------------------------------------------------------------------
    {
        const double A = 0.25725664;
        const double B = 0.04276376;
        const double c2 = -(1.0 - B);
        const double c1 = A - 3.0 * B * B - 2.0 * B;
        const double c0 = -(A * B - B * B - B * B * B);
        const auto sol = thermo::CubicSolver::solve(c2, c1, c0);
        report_test("Test 6 (PR CO2 4MPa): Exactly 1 real root", sol.num_real_roots == 1);
        report_test("Test 6 (PR CO2 4MPa): Root Z == 0.752334", approx_equal(sol.roots[0], 0.75233429, 1.0e-5));
        report_test("Test 6 (PR CO2 4MPa): Residual <= 1e-12",
                    sol.residuals[0] <= thermo::constants::CUBIC_POLY_RESIDUAL_TOL);
    }

    // ------------------------------------------------------------------------
    // Test 7: Audited PR-1976 Case 3: CO2 at 300 K, 7 MPa
    // Dimensionless parameters: A = 0.45019912, B = 0.07483658
    // Expected root: Z = 0.198057 (Single compressed liquid root, Delta > 0)
    // ------------------------------------------------------------------------
    {
        const double A = 0.45019912;
        const double B = 0.07483658;
        const double c2 = -(1.0 - B);
        const double c1 = A - 3.0 * B * B - 2.0 * B;
        const double c0 = -(A * B - B * B - B * B * B);
        const auto sol = thermo::CubicSolver::solve(c2, c1, c0);
        report_test("Test 7 (PR CO2 7MPa): Exactly 1 real root (Delta > 0)", sol.num_real_roots == 1);
        report_test("Test 7 (PR CO2 7MPa): Root Z == 0.198057", approx_equal(sol.roots[0], 0.19805743, 1.0e-5));
        report_test("Test 7 (PR CO2 7MPa): Residual <= 1e-12",
                    sol.residuals[0] <= thermo::constants::CUBIC_POLY_RESIDUAL_TOL);
    }

    // ------------------------------------------------------------------------
    // Test 8: Audited PR-1976 Case 4: Propane at 320 K, 1 MPa
    // Dimensionless parameters: A = 0.15613051, B = 0.02116534
    // Expected 3 real roots: Z_0 = 0.036453, Z_1 = 0.091826, Z_2 = 0.850555
    // ------------------------------------------------------------------------
    {
        const double A = 0.15613051;
        const double B = 0.02116534;
        const double c2 = -(1.0 - B);
        const double c1 = A - 3.0 * B * B - 2.0 * B;
        const double c0 = -(A * B - B * B - B * B * B);
        const auto sol = thermo::CubicSolver::solve(c2, c1, c0);
        report_test("Test 8 (PR Propane 1MPa): 3 real roots (Delta < 0)", sol.num_real_roots == 3);
        report_test("Test 8 (PR Propane 1MPa): 3 distinct roots", sol.num_distinct_roots == 3);
        report_test("Test 8 (PR Propane 1MPa): Root 0 == 0.036453", approx_equal(sol.roots[0], 0.03645310, 1.0e-5));
        report_test("Test 8 (PR Propane 1MPa): Root 1 == 0.091826", approx_equal(sol.roots[1], 0.09182609, 1.0e-5));
        report_test("Test 8 (PR Propane 1MPa): Root 2 == 0.850555", approx_equal(sol.roots[2], 0.85055548, 1.0e-5));
        for (size_t i = 0; i < sol.residuals.size(); ++i) {
            report_test("Test 8 (PR Propane 1MPa): Residual " + std::to_string(i) + " <= 1e-12",
                        sol.residuals[i] <= thermo::constants::CUBIC_POLY_RESIDUAL_TOL);
        }
    }

    // ------------------------------------------------------------------------
    // Test 9: Canonical PR Critical Point Triple Root (T = Tc, P = Pc)
    // Dimensionless parameters: A = Omega_a, B = Omega_b
    // Expected: p == 0, q == 0, Delta == 0, Triple root Z_c = -c2/3 ~ 0.30740131
    // ------------------------------------------------------------------------
    {
        const double A = thermo::constants::OMEGA_A;
        const double B = thermo::constants::OMEGA_B;
        const double c2 = -(1.0 - B);
        const double c1 = A - 3.0 * B * B - 2.0 * B;
        const double c0 = -(A * B - B * B - B * B * B);
        const auto sol = thermo::CubicSolver::solve(c2, c1, c0);

        report_test("Test 9 (PR Critical Triple Root): 3 real roots", sol.num_real_roots == 3);
        report_test("Test 9 (PR Critical Triple Root): Exactly 1 distinct root (multiplicity 3)",
                    sol.num_distinct_roots == 1);
        report_test("Test 9 (PR Critical Triple Root): Z_c ~ 0.30740131",
                    approx_equal(sol.distinct_roots[0], 0.30740131, 1.0e-6));
        report_test("Test 9 (PR Critical Triple Root): Residual <= 1e-12",
                    sol.residuals[0] <= thermo::constants::CUBIC_POLY_RESIDUAL_TOL);
    }

    // ------------------------------------------------------------------------
    // Test 10: Near-Critical State (T = 1.001*Tc, P = Pc)
    // Testing numerical stability of Cardano solver near discriminant cusp.
    // Physics note: near critical point, (Z - Z_c) ~ (Tr - 1)^(1/3) ~ (0.001)^(1/3) ~ 0.1,
    // so a 0.1% temperature increase shifts Z from 0.3074 to 0.3638.
    // ------------------------------------------------------------------------
    {
        const double Tr_near = 1.001;
        const double A_near = thermo::constants::OMEGA_A / (Tr_near * Tr_near);
        const double B_near = thermo::constants::OMEGA_B / Tr_near;
        const double c2 = -(1.0 - B_near);
        const double c1 = A_near - 3.0 * B_near * B_near - 2.0 * B_near;
        const double c0 = -(A_near * B_near - B_near * B_near - B_near * B_near * B_near);
        const auto sol = thermo::CubicSolver::solve(c2, c1, c0);

        report_test("Test 10 (Near-Critical Coalescence): 1 real root (Delta > 0)", sol.num_real_roots == 1);
        report_test("Test 10 (Near-Critical Coalescence): Root Z == 0.363838",
                    approx_equal(sol.roots[0], 0.363838, 1.0e-4));
        report_test("Test 10 (Near-Critical Coalescence): Residual <= 1e-12",
                    sol.residuals[0] <= thermo::constants::CUBIC_POLY_RESIDUAL_TOL);
    }

    std::cout << "==================================================\n";
    std::cout << "SUMMARY: Passed " << g_pass_count << ", Failed " << g_fail_count << "\n";
    std::cout << "==================================================\n";

    return (g_fail_count == 0) ? 0 : 1;
}
