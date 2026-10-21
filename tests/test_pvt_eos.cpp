#include "thermo/PengRobinsonEOS.hpp"
#include "thermo/ComponentDatabase.hpp"
#include "thermo/Constants.hpp"

#include <iostream>
#include <cassert>
#include <cmath>
#include <string>
#include <limits>

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

bool approx_equal(double a, double b, double tol = 1.0e-5) {
    return std::abs(a - b) <= tol;
}

bool approx_equal_rel(double a, double b, double rel_tol = 1.0e-7) {
    const double diff = std::abs(a - b);
    const double scale = std::max(std::abs(a), std::abs(b));
    if (scale == 0.0) return diff == 0.0;
    return (diff / scale) <= rel_tol;
}

} // namespace

int main() {
    std::cout << "==================================================\n";
    std::cout << "RUNNING PENG-ROBINSON EOS & PVT TESTS (Session 3)\n";
    std::cout << "==================================================\n";

    thermo::ComponentDatabase db;
    const auto& methane = db.get_component("CH4");
    const auto& co2 = db.get_component("CO2");
    const auto& propane = db.get_component("C3H8");

    thermo::PengRobinsonEOS pr_ch4(methane);
    thermo::PengRobinsonEOS pr_co2(co2);
    thermo::PengRobinsonEOS pr_c3h8(propane);

    // ========================================================================
    // Test A: Published-Form Parameter Consistency
    // ========================================================================
    {
        std::cout << "\n--- Test A: Published-Form Equation Consistency ---\n";
        // Check Methane parameters
        const double ac_ch4 = pr_ch4.ac();
        const double b_ch4 = pr_ch4.b();
        const double kappa_ch4 = pr_ch4.kappa();

        report_test("Test A: Methane a_c > 0", ac_ch4 > 0.0);
        report_test("Test A: Methane b > 0", b_ch4 > 0.0);
        report_test("Test A: Methane kappa > 0", kappa_ch4 > 0.0);
        report_test("Test A: Methane a_c ~ 0.24958 Pa*m^6/mol^2", approx_equal(ac_ch4, 0.249579, 1.0e-4));
        report_test("Test A: Methane b ~ 2.680e-5 m^3/mol", approx_equal(b_ch4, 2.680157e-5, 1.0e-8));

        // Dimensionless parameters A and B at T = 300 K, P = 5 MPa
        const double A_ch4 = pr_ch4.A_parameter(300.0, 5.0e6);
        const double B_ch4 = pr_ch4.B_parameter(300.0, 5.0e6);
        report_test("Test A: Methane A == 0.162486", approx_equal(A_ch4, 0.16248591, 1.0e-6));
        report_test("Test A: Methane B == 0.053725", approx_equal(B_ch4, 0.05372479, 1.0e-6));

        // Methane analytical derivative da/dT at T = 300 K
        const double da_dT_ch4 = pr_ch4.da_dT(300.0);
        report_test("Test A: da/dT is strictly negative (molecular attraction weakens with T)", da_dT_ch4 < 0.0);
        report_test("Test A: Methane da/dT == -3.686e-4", approx_equal(da_dT_ch4, -3.686115e-4, 1.0e-7));

        // EOS parameters container
        const auto prm_co2 = pr_co2.parameters(300.0);
        report_test("Test A: CO2 Tr(300 K) ~ 0.98642", approx_equal(prm_co2.Tr, 0.986420, 1.0e-5));
        report_test("Test A: CO2 alpha(300 K) ~ 1.00965", approx_equal(prm_co2.alpha, 1.009649, 1.0e-5));
        report_test("Test A: CO2 a(300 K) ~ 0.40015", approx_equal(prm_co2.a, 0.4001451, 1.0e-5));
    }

    // ========================================================================
    // Test B: Forward/Inverse Reciprocity P(T, v(T, P)) == P
    // ========================================================================
    {
        std::cout << "\n--- Test B: Forward/Inverse PVT Reciprocity ---\n";
        // Test multiple substances, temperatures, and pressures
        const struct ReciprocityCase {
            const thermo::PengRobinsonEOS& eos;
            std::string name;
            double T; // K
            double P; // Pa
        } cases[] = {
            {pr_ch4, "Methane", 200.0, 1.0e5},
            {pr_ch4, "Methane", 300.0, 1.0e6},
            {pr_ch4, "Methane", 300.0, 5.0e6},
            {pr_ch4, "Methane", 400.0, 2.0e7},
            {pr_co2, "CO2", 280.0, 2.0e6},
            {pr_co2, "CO2", 300.0, 4.0e6},
            {pr_co2, "CO2", 300.0, 7.0e6},
            {pr_co2, "CO2", 305.0, 7.4e6}, // Near critical
            {pr_co2, "CO2", 500.0, 1.0e7},
            {pr_c3h8, "Propane", 320.0, 1.0e6},
            {pr_c3h8, "Propane", 350.0, 3.0e6},
            {pr_c3h8, "Propane", 400.0, 5.0e6}
        };

        for (const auto& c : cases) {
            const auto pvt = c.eos.calculate_pvt(c.T, c.P);
            report_test("Test B: " + c.name + " at T=" + std::to_string(static_cast<int>(c.T)) +
                        "K, P=" + std::to_string(static_cast<int>(c.P/1e5)) + "bar has physical roots",
                        !pvt.roots.empty());

            for (size_t i = 0; i < pvt.roots.size(); ++i) {
                const double v = pvt.roots[i].molar_volume;
                const double P_calc = c.eos.calculate_pressure(c.T, v);
                const double rel_err = std::abs(P_calc - c.P) / c.P;
                report_test("Test B: " + c.name + " root " + std::to_string(i) + " reciprocity (|P_calc - P|/P < 1e-10)",
                            rel_err <= thermo::constants::PVT_RECIPROCITY_TOL);
            }
        }
    }

    // ========================================================================
    // Test C: Ideal-Gas Limit as P -> 0
    // ========================================================================
    {
        std::cout << "\n--- Test C: Ideal-Gas Limit as P -> 0 ---\n";
        const double T = 300.0;
        const double low_pressures[] = {100.0, 10.0, 1.0, 0.01, 1.0e-4};

        for (double P : low_pressures) {
            const auto pvt = pr_ch4.calculate_pvt(T, P);
            report_test("Test C: Exactly 1 root at P = " + std::to_string(P) + " Pa", pvt.roots.size() == 1);

            const double Z = pvt.roots[0].Z;
            const double v_real = pvt.roots[0].molar_volume;
            const double v_ideal = (thermo::constants::R * T) / P;

            // Z must approach 1.0
            const double Z_diff = std::abs(Z - 1.0);
            // v must approach v_ideal
            const double v_rel_err = std::abs(v_real - v_ideal) / v_ideal;

            report_test("Test C (P=" + std::to_string(P) + " Pa): |Z - 1| approaches 0", Z_diff <= (1.0e-4 * (P / 100.0) + 1.0e-7));
            report_test("Test C (P=" + std::to_string(P) + " Pa): v approaches v_ideal", v_rel_err <= (1.0e-4 * (P / 100.0) + 1.0e-7));

            // Pressure from ideal volume approaches P
            const double P_from_v_ideal = pr_ch4.calculate_pressure(T, v_ideal);
            const double P_rel_err = std::abs(P_from_v_ideal - P) / P;
            report_test("Test C (P=" + std::to_string(P) + " Pa): P(T, v_ideal) ~ P", P_rel_err <= (1.0e-4 * (P / 100.0) + 1.0e-6));
        }
    }

    // ========================================================================
    // Test D: Root Topology & Audited Benchmark States
    // ========================================================================
    {
        std::cout << "\n--- Test D: Audited Benchmark States & Root Topology ---\n";

        // Case 1: Methane at 300 K, 5 MPa (Single supercritical vapor root)
        {
            const auto pvt = pr_ch4.calculate_pvt(300.0, 5.0e6);
            report_test("Test D (BM-01 Methane 5MPa): Discriminant > 0", pvt.discriminant > 0.0);
            report_test("Test D (BM-01 Methane 5MPa): Exactly 1 physical root", pvt.roots.size() == 1);
            report_test("Test D (BM-01 Methane 5MPa): Z == 0.901844", approx_equal(pvt.roots[0].Z, 0.90184394, 1.0e-5));
            report_test("Test D (BM-01 Methane 5MPa): v == 4.4990e-4 m^3/mol",
                        approx_equal(pvt.roots[0].molar_volume, 4.499009e-4, 1.0e-7));
            report_test("Test D (BM-01 Methane 5MPa): Phase is SINGLE_PHASE",
                        pvt.roots[0].phase == thermo::PhaseTag::SINGLE_PHASE);
            report_test("Test D (BM-01 Methane 5MPa): Mechanically STABLE",
                        pvt.roots[0].stability == thermo::MechanicalStability::STABLE);
        }

        // Case 2: CO2 at 300 K, 4 MPa (Single subcritical vapor root)
        {
            const auto pvt = pr_co2.calculate_pvt(300.0, 4.0e6);
            report_test("Test D (BM-02 CO2 4MPa): Discriminant > 0", pvt.discriminant > 0.0);
            report_test("Test D (BM-02 CO2 4MPa): Exactly 1 physical root", pvt.roots.size() == 1);
            report_test("Test D (BM-02 CO2 4MPa): Z == 0.752334", approx_equal(pvt.roots[0].Z, 0.75233429, 1.0e-5));
            report_test("Test D (BM-02 CO2 4MPa): v == 4.6914e-4 m^3/mol",
                        approx_equal(pvt.roots[0].molar_volume, 4.691441e-4, 1.0e-7));
            report_test("Test D (BM-02 CO2 4MPa): Phase is SINGLE_PHASE",
                        pvt.roots[0].phase == thermo::PhaseTag::SINGLE_PHASE);
            report_test("Test D (BM-02 CO2 4MPa): Mechanically STABLE",
                        pvt.roots[0].stability == thermo::MechanicalStability::STABLE);
        }

        // Case 3: CO2 at 300 K, 7 MPa (Single compressed liquid root, Delta > 0)
        {
            const auto pvt = pr_co2.calculate_pvt(300.0, 7.0e6);
            report_test("Test D (BM-03 CO2 7MPa): Discriminant > 0", pvt.discriminant > 0.0);
            report_test("Test D (BM-03 CO2 7MPa): Exactly 1 physical root", pvt.roots.size() == 1);
            report_test("Test D (BM-03 CO2 7MPa): Z == 0.198057", approx_equal(pvt.roots[0].Z, 0.19805743, 1.0e-5));
            report_test("Test D (BM-03 CO2 7MPa): v == 7.0575e-5 m^3/mol",
                        approx_equal(pvt.roots[0].molar_volume, 7.057462e-5, 1.0e-7));
            report_test("Test D (BM-03 CO2 7MPa): Phase is SINGLE_PHASE",
                        pvt.roots[0].phase == thermo::PhaseTag::SINGLE_PHASE);
            report_test("Test D (BM-03 CO2 7MPa): Mechanically STABLE",
                        pvt.roots[0].stability == thermo::MechanicalStability::STABLE);
        }

        // Case 4: Propane at 320 K, 1 MPa (Three real physical roots inside spinodal)
        {
            const auto pvt = pr_c3h8.calculate_pvt(320.0, 1.0e6);
            report_test("Test D (BM-04 Propane 1MPa): Discriminant < 0", pvt.discriminant < 0.0);
            report_test("Test D (BM-04 Propane 1MPa): Exactly 3 physical roots", pvt.roots.size() == 3);
            report_test("Test D (BM-04 Propane 1MPa): Root 0 (Liquid) Z == 0.036453",
                        approx_equal(pvt.roots[0].Z, 0.03645310, 1.0e-5));
            report_test("Test D (BM-04 Propane 1MPa): Root 1 (Middle) Z == 0.091826",
                        approx_equal(pvt.roots[1].Z, 0.09182609, 1.0e-5));
            report_test("Test D (BM-04 Propane 1MPa): Root 2 (Vapor) Z == 0.850555",
                        approx_equal(pvt.roots[2].Z, 0.85055548, 1.0e-5));

            // Phase tags
            report_test("Test D (BM-04 Propane 1MPa): Root 0 tagged as LIQUID",
                        pvt.roots[0].phase == thermo::PhaseTag::LIQUID);
            report_test("Test D (BM-04 Propane 1MPa): Root 1 tagged as INTERMEDIATE",
                        pvt.roots[1].phase == thermo::PhaseTag::INTERMEDIATE);
            report_test("Test D (BM-04 Propane 1MPa): Root 2 tagged as VAPOR",
                        pvt.roots[2].phase == thermo::PhaseTag::VAPOR);

            // Stability tags
            report_test("Test D (BM-04 Propane 1MPa): Liquid root is STABLE",
                        pvt.roots[0].stability == thermo::MechanicalStability::STABLE);
            report_test("Test D (BM-04 Propane 1MPa): Intermediate root is UNSTABLE",
                        pvt.roots[1].stability == thermo::MechanicalStability::UNSTABLE);
            report_test("Test D (BM-04 Propane 1MPa): Vapor root is STABLE",
                        pvt.roots[2].stability == thermo::MechanicalStability::STABLE);

            // Accessors
            report_test("Test D (BM-04 Propane 1MPa): liquid_root() accessor works",
                        pvt.liquid_root() != nullptr && approx_equal(pvt.liquid_root()->Z, 0.036453, 1e-4));
            report_test("Test D (BM-04 Propane 1MPa): vapor_root() accessor works",
                        pvt.vapor_root() != nullptr && approx_equal(pvt.vapor_root()->Z, 0.850555, 1e-4));
            report_test("Test D (BM-04 Propane 1MPa): intermediate_root() accessor works",
                        pvt.intermediate_root() != nullptr && approx_equal(pvt.intermediate_root()->Z, 0.091826, 1e-4));
        }
    }

    // ========================================================================
    // Test E: Mechanical Stability Derivative Validation
    // ========================================================================
    {
        std::cout << "\n--- Test E: Mechanical Stability Derivative Cross-Validation ---\n";
        // Check analytical (dP/dv)_T against central finite difference
        const double T = 320.0;
        const auto pvt = pr_c3h8.calculate_pvt(T, 1.0e6);
        assert(pvt.roots.size() == 3);

        for (size_t i = 0; i < pvt.roots.size(); ++i) {
            const double v = pvt.roots[i].molar_volume;
            const double dP_dv_ana = pr_c3h8.calculate_dP_dv(T, v);
            report_test("Test E: Root " + std::to_string(i) + " stored dP_dv matches calculate_dP_dv",
                        approx_equal_rel(pvt.roots[i].dP_dv, dP_dv_ana, 1.0e-12));

            // Central difference with step delta_v = 1e-5 * v
            const double delta_v = 1.0e-5 * v;
            const double P_plus = pr_c3h8.calculate_pressure(T, v + delta_v);
            const double P_minus = pr_c3h8.calculate_pressure(T, v - delta_v);
            const double dP_dv_num = (P_plus - P_minus) / (2.0 * delta_v);

            report_test("Test E: Root " + std::to_string(i) + " analytical vs finite-difference dP/dv relative error <= 1e-5",
                        approx_equal_rel(dP_dv_ana, dP_dv_num, 1.0e-5));
        }

        // Verify sign: liquid < 0, middle > 0, vapor < 0
        report_test("Test E: Liquid root dP/dv < 0", pvt.roots[0].dP_dv < 0.0);
        report_test("Test E: Intermediate root dP/dv > 0", pvt.roots[1].dP_dv > 0.0);
        report_test("Test E: Vapor root dP/dv < 0", pvt.roots[2].dP_dv < 0.0);
    }

    // ========================================================================
    // Test F: Domain Validation & Invalid Inputs
    // ========================================================================
    {
        std::cout << "\n--- Test F: Domain Validation & Invalid Inputs ---\n";

        // Non-positive T in forward PVT
        bool threw = false;
        try { (void)pr_ch4.calculate_pvt(0.0, 1.0e5); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test F: T = 0 K in calculate_pvt throws invalid_argument", threw);

        threw = false;
        try { (void)pr_ch4.calculate_pvt(-100.0, 1.0e5); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test F: T < 0 K in calculate_pvt throws invalid_argument", threw);

        // Non-positive P in forward PVT
        threw = false;
        try { (void)pr_ch4.calculate_pvt(300.0, 0.0); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test F: P = 0 Pa in calculate_pvt throws invalid_argument", threw);

        threw = false;
        try { (void)pr_ch4.calculate_pvt(300.0, -1.0e5); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test F: P < 0 Pa in calculate_pvt throws invalid_argument", threw);

        // Non-positive T in calculate_pressure
        threw = false;
        try { (void)pr_ch4.calculate_pressure(0.0, 1.0e-3); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test F: T = 0 K in calculate_pressure throws invalid_argument", threw);

        // v <= b in calculate_pressure
        const double b = pr_ch4.b();
        threw = false;
        try { (void)pr_ch4.calculate_pressure(300.0, b); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test F: v == b in calculate_pressure throws invalid_argument", threw);

        threw = false;
        try { (void)pr_ch4.calculate_pressure(300.0, b * 0.9); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test F: v < b in calculate_pressure throws invalid_argument", threw);

        // Non-finite values
        const double nan_val = std::numeric_limits<double>::quiet_NaN();
        const double inf_val = std::numeric_limits<double>::infinity();

        threw = false;
        try { (void)pr_ch4.calculate_pvt(nan_val, 1.0e5); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test F: NaN T in calculate_pvt throws invalid_argument", threw);

        threw = false;
        try { (void)pr_ch4.calculate_pvt(300.0, inf_val); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test F: Inf P in calculate_pvt throws invalid_argument", threw);

        threw = false;
        try { (void)pr_ch4.calculate_pressure(300.0, nan_val); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test F: NaN v in calculate_pressure throws invalid_argument", threw);
    }

    std::cout << "==================================================\n";
    std::cout << "SUMMARY: Passed " << g_pass_count << ", Failed " << g_fail_count << "\n";
    std::cout << "==================================================\n";

    return (g_fail_count == 0) ? 0 : 1;
}
