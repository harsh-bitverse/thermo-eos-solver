#include "thermo/PengRobinsonEOS.hpp"
#include "thermo/ComponentDatabase.hpp"
#include "thermo/Constants.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>

namespace {

static int total_assertions = 0;
static int passed_assertions = 0;

#define TEST_ASSERT(cond, msg) \
    do { \
        total_assertions++; \
        if (cond) { \
            passed_assertions++; \
            std::cout << "[PASS] " << msg << "\n"; \
        } else { \
            std::cerr << "[FAIL] Line " << __LINE__ << ": " << msg << "\n"; \
            assert(false); \
        } \
    } while (0)

bool approx_equal(double val, double expected, double tol) {
    return std::abs(val - expected) <= tol;
}

bool approx_rel_equal(double val, double expected, double rel_tol) {
    if (std::abs(expected) < 1.0e-12) {
        return std::abs(val - expected) <= rel_tol;
    }
    return (std::abs(val - expected) / std::abs(expected)) <= rel_tol;
}

} // anonymous namespace

int main() {
    std::cout << "========================================================\n";
    std::cout << "RUNNING THERMODYNAMIC PROPERTIES & EQUILIBRIUM TESTS (S4)\n";
    std::cout << "========================================================\n\n";

    thermo::ComponentDatabase db;

    // ========================================================================
    // Part A: Residual-Property Identities & Maxwell Relation Verification
    // ========================================================================
    std::cout << "--- Part A: Residual-Property Identities & Maxwell Relation ---\n";
    {
        // 1. Methane BM-01: 300 K, 5 MPa
        const auto comp_ch4 = db.get_component("Methane");
        const thermo::PengRobinsonEOS eos_ch4(comp_ch4);
        const double T = 300.0;
        const double P = 5.0e6;
        const auto pvt = eos_ch4.calculate_pvt(T, P);

        TEST_ASSERT(pvt.roots.size() == 1, "Methane 300K, 5MPa has exactly 1 root");
        const auto& r = pvt.roots[0];

        TEST_ASSERT(std::isfinite(r.H_res), "Methane H_res is finite");
        TEST_ASSERT(std::isfinite(r.S_res), "Methane S_res is finite");
        TEST_ASSERT(std::isfinite(r.ln_phi), "Methane ln_phi is finite");
        TEST_ASSERT(r.phi > 0.0, "Methane phi is strictly positive");
        TEST_ASSERT(r.fugacity > 0.0, "Methane fugacity is strictly positive");

        // Audited benchmark comparison
        TEST_ASSERT(approx_rel_equal(r.Z, 0.901844, 1.0e-5), "Methane BM-01: Z ~ 0.901844");
        TEST_ASSERT(approx_equal(r.H_res, -902.3382, 0.05), "Methane BM-01: H_res ~ -902.34 J/mol");
        TEST_ASSERT(approx_equal(r.S_res, -2.144565, 0.0005), "Methane BM-01: S_res ~ -2.1446 J/(mol*K)");
        TEST_ASSERT(approx_rel_equal(r.phi, 0.901385, 1.0e-5), "Methane BM-01: phi ~ 0.901385");
        TEST_ASSERT(approx_rel_equal(r.fugacity, 4.506926e6, 1.0e-5), "Methane BM-01: f ~ 4.506926 MPa");

        // Rigorous Maxwell identity: S_res = H_res / T - R * ln(phi)
        const double maxwell_rhs = r.H_res / T - thermo::constants::R * r.ln_phi;
        const double identity_err = std::abs(r.S_res - maxwell_rhs);
        TEST_ASSERT(identity_err <= 1.0e-12, "Methane BM-01: S_res == H_res/T - R*ln(phi) within 1e-12");

        // 2. CO2 BM-02: 300 K, 4 MPa
        const auto comp_co2 = db.get_component("CO2");
        const thermo::PengRobinsonEOS eos_co2(comp_co2);
        const auto pvt_co2_4 = eos_co2.calculate_pvt(300.0, 4.0e6);
        TEST_ASSERT(pvt_co2_4.roots.size() == 1, "CO2 300K, 4MPa has exactly 1 root");
        const auto& r_co2_4 = pvt_co2_4.roots[0];

        TEST_ASSERT(approx_rel_equal(r_co2_4.Z, 0.752334, 1.0e-5), "CO2 BM-02: Z ~ 0.752334");
        TEST_ASSERT(approx_equal(r_co2_4.H_res, -1990.9877, 0.05), "CO2 BM-02: H_res ~ -1990.99 J/mol");
        TEST_ASSERT(approx_equal(r_co2_4.S_res, -4.734694, 0.0005), "CO2 BM-02: S_res ~ -4.7347 J/(mol*K)");
        TEST_ASSERT(approx_rel_equal(r_co2_4.phi, 0.795528, 1.0e-5), "CO2 BM-02: phi ~ 0.795528");
        TEST_ASSERT(approx_rel_equal(r_co2_4.fugacity, 3.182110e6, 1.0e-5), "CO2 BM-02: f ~ 3.182110 MPa");

        const double maxwell_rhs_co2_4 = r_co2_4.H_res / 300.0 - thermo::constants::R * r_co2_4.ln_phi;
        TEST_ASSERT(std::abs(r_co2_4.S_res - maxwell_rhs_co2_4) <= 1.0e-12,
                    "CO2 BM-02: S_res == H_res/T - R*ln(phi) within 1e-12");

        // 3. CO2 BM-03: 300 K, 7 MPa (single liquid root)
        const auto pvt_co2_7 = eos_co2.calculate_pvt(300.0, 7.0e6);
        TEST_ASSERT(pvt_co2_7.roots.size() == 1, "CO2 300K, 7MPa has exactly 1 root");
        const auto& r_co2_7 = pvt_co2_7.roots[0];

        TEST_ASSERT(approx_rel_equal(r_co2_7.Z, 0.198057, 1.0e-5), "CO2 BM-03: Z ~ 0.198057");
        TEST_ASSERT(approx_equal(r_co2_7.H_res, -9374.3721, 0.05), "CO2 BM-03: H_res ~ -9374.37 J/mol");
        TEST_ASSERT(approx_equal(r_co2_7.S_res, -27.514960, 0.0005), "CO2 BM-03: S_res ~ -27.5150 J/(mol*K)");
        TEST_ASSERT(approx_rel_equal(r_co2_7.phi, 0.638285, 1.0e-5), "CO2 BM-03: phi ~ 0.638285");
        TEST_ASSERT(approx_rel_equal(r_co2_7.fugacity, 4.467995e6, 1.0e-5), "CO2 BM-03: f ~ 4.467995 MPa");

        const double maxwell_rhs_co2_7 = r_co2_7.H_res / 300.0 - thermo::constants::R * r_co2_7.ln_phi;
        TEST_ASSERT(std::abs(r_co2_7.S_res - maxwell_rhs_co2_7) <= 1.0e-12,
                    "CO2 BM-03: S_res == H_res/T - R*ln(phi) within 1e-12");
    }

    // ========================================================================
    // Part B: Ideal-Gas & Low-Pressure Asymptotics (P -> 0)
    // ========================================================================
    std::cout << "\n--- Part B: Ideal-Gas & Low-Pressure Asymptotics ---\n";
    {
        const auto comp_ch4 = db.get_component("Methane");
        const thermo::PengRobinsonEOS eos_ch4(comp_ch4);
        const double T = 300.0;

        // Sequence of decreasing pressures: 100 kPa down to 10^-5 Pa
        const std::vector<double> pressures = {100000.0, 1000.0, 10.0, 0.1, 1.0e-3, 1.0e-5};

        for (double P_test : pressures) {
            const auto pvt = eos_ch4.calculate_pvt(T, P_test);
            TEST_ASSERT(pvt.roots.size() == 1, "P = " + std::to_string(P_test) + " Pa: exactly 1 root");
            const auto& r = pvt.roots[0];

            TEST_ASSERT(std::isfinite(r.H_res) && std::isfinite(r.S_res),
                        "P = " + std::to_string(P_test) + " Pa: H_res and S_res are finite");
            TEST_ASSERT(std::isfinite(r.ln_phi) && r.phi > 0.0,
                        "P = " + std::to_string(P_test) + " Pa: ln_phi finite and phi > 0");

            // Reciprocity f / P == phi
            TEST_ASSERT(approx_rel_equal(r.fugacity / P_test, r.phi, 1.0e-14),
                        "P = " + std::to_string(P_test) + " Pa: f / P == phi");

            // Maxwell identity at every pressure step
            const double maxwell_rhs = r.H_res / T - thermo::constants::R * r.ln_phi;
            TEST_ASSERT(std::abs(r.S_res - maxwell_rhs) <= 1.0e-12,
                        "P = " + std::to_string(P_test) + " Pa: Maxwell identity holds");
        }

        // Specifically check the ultra-low pressure state (10^-5 Pa)
        const auto pvt_ultra = eos_ch4.calculate_pvt(T, 1.0e-5);
        const auto& r_ultra = pvt_ultra.roots[0];
        TEST_ASSERT(approx_equal(r_ultra.Z, 1.0, 1.0e-9), "P = 1e-5 Pa: Z -> 1.0 exactly");
        TEST_ASSERT(approx_equal(r_ultra.phi, 1.0, 1.0e-9), "P = 1e-5 Pa: phi -> 1.0 exactly");
        TEST_ASSERT(approx_equal(r_ultra.H_res, 0.0, 1.0e-5), "P = 1e-5 Pa: H_res -> 0.0 J/mol");
        TEST_ASSERT(approx_equal(r_ultra.S_res, 0.0, 1.0e-5), "P = 1e-5 Pa: S_res -> 0.0 J/(mol*K)");

        // Continuity across the small-B Taylor regularization threshold (u = 1.0e-5)
        // Let's test calculate_I_L directly across the threshold
        const double Z_test = 1.0;
        const double B_below = 0.999e-5 / (2.0 * thermo::constants::SQRT2);
        const double B_above = 1.001e-5 / (2.0 * thermo::constants::SQRT2);

        const double IL_below = eos_ch4.calculate_I_L(Z_test, B_below);
        const double IL_above = eos_ch4.calculate_I_L(Z_test, B_above);

        TEST_ASSERT(std::isfinite(IL_below) && std::isfinite(IL_above), "I_L finite across threshold");
        TEST_ASSERT(approx_rel_equal(IL_below, IL_above, 1.0e-4),
                    "I_L is smooth and continuous across small-B Taylor boundary");
    }

    // ========================================================================
    // Part C: Pure-Component Fugacity and Three-Root Multiplicity
    // ========================================================================
    std::cout << "\n--- Part C: Pure-Component Fugacity & Three-Root Evaluation ---\n";
    {
        // Propane BM-04: 320 K, 1 MPa (Three real roots)
        const auto comp_c3h8 = db.get_component("Propane");
        const thermo::PengRobinsonEOS eos_c3h8(comp_c3h8);
        const double T = 320.0;
        const double P = 1.0e6;
        const auto pvt = eos_c3h8.calculate_pvt(T, P);

        TEST_ASSERT(pvt.roots.size() == 3, "Propane BM-04: exactly 3 roots found");

        const auto* liq = pvt.liquid_root();
        const auto* mid = pvt.intermediate_root();
        const auto* vap = pvt.vapor_root();

        TEST_ASSERT(liq != nullptr, "Liquid root exists");
        TEST_ASSERT(mid != nullptr, "Intermediate root exists");
        TEST_ASSERT(vap != nullptr, "Vapor root exists");

        // Verify root 0 (Liquid)
        TEST_ASSERT(approx_rel_equal(liq->Z, 0.036453, 1.0e-4), "Propane Root 0 (Liq): Z ~ 0.036453");
        TEST_ASSERT(approx_equal(liq->H_res, -14854.34, 0.5), "Propane Root 0 (Liq): H_res ~ -14854.34 J/mol");
        TEST_ASSERT(approx_equal(liq->S_res, -48.2035, 0.01), "Propane Root 0 (Liq): S_res ~ -48.2035 J/(mol*K)");
        TEST_ASSERT(approx_rel_equal(liq->phi, 1.239280, 1.0e-4), "Propane Root 0 (Liq): phi ~ 1.239280");
        TEST_ASSERT(approx_rel_equal(liq->fugacity, 1.239280e6, 1.0e-4), "Propane Root 0 (Liq): f ~ 1.239280 MPa");

        // Verify root 1 (Intermediate)
        TEST_ASSERT(approx_rel_equal(mid->Z, 0.091826, 1.0e-4), "Propane Root 1 (Mid): Z ~ 0.091826");
        TEST_ASSERT(approx_equal(mid->H_res, -8210.60, 0.5), "Propane Root 1 (Mid): H_res ~ -8210.60 J/mol");
        TEST_ASSERT(approx_equal(mid->S_res, -28.3698, 0.01), "Propane Root 1 (Mid): S_res ~ -28.3698 J/(mol*K)");
        TEST_ASSERT(approx_rel_equal(mid->phi, 1.385615, 1.0e-4), "Propane Root 1 (Mid): phi ~ 1.385615");
        TEST_ASSERT(mid->stability == thermo::MechanicalStability::UNSTABLE, "Propane Root 1 (Mid) is UNSTABLE");

        // Verify root 2 (Vapor)
        TEST_ASSERT(approx_rel_equal(vap->Z, 0.850555, 1.0e-4), "Propane Root 2 (Vap): Z ~ 0.850555");
        TEST_ASSERT(approx_equal(vap->H_res, -1131.04, 0.5), "Propane Root 2 (Vap): H_res ~ -1131.04 J/mol");
        TEST_ASSERT(approx_equal(vap->S_res, -2.3575, 0.01), "Propane Root 2 (Vap): S_res ~ -2.3575 J/(mol*K)");
        TEST_ASSERT(approx_rel_equal(vap->phi, 0.868006, 1.0e-4), "Propane Root 2 (Vap): phi ~ 0.868006");
        TEST_ASSERT(approx_rel_equal(vap->fugacity, 868006.2, 1.0e-4), "Propane Root 2 (Vap): f ~ 0.868006 MPa");

        // Maxwell identity holds on all 3 roots
        for (size_t i = 0; i < pvt.roots.size(); ++i) {
            const auto& r = pvt.roots[i];
            const double rhs = r.H_res / T - thermo::constants::R * r.ln_phi;
            TEST_ASSERT(std::abs(r.S_res - rhs) <= 1.0e-12,
                        "Propane Root " + std::to_string(i) + ": S_res == H_res/T - R*ln(phi) within 1e-12");
        }
    }

    // ========================================================================
    // Part D: Total Enthalpy & Entropy (Reference State Convention)
    // ========================================================================
    std::cout << "\n--- Part D: Total H & S (Reference State Convention) ---\n";
    {
        const auto comp_ch4 = db.get_component("Methane");
        const thermo::PengRobinsonEOS eos_ch4(comp_ch4);

        // 1. Exactly at Reference State (298.15 K, 100,000 Pa):
        // Total H and S must be zero within numerical tolerance
        const auto pvt_ref = eos_ch4.calculate_pvt(thermo::constants::T_REF, thermo::constants::P_REF);
        TEST_ASSERT(pvt_ref.roots.size() == 1, "Methane at ref state has 1 root");
        const auto& r_ref = pvt_ref.roots[0];
        TEST_ASSERT(r_ref.H_total.has_value(), "Ref state root has H_total value");
        TEST_ASSERT(r_ref.S_total.has_value(), "Ref state root has S_total value");
        TEST_ASSERT(std::abs(*r_ref.H_total) <= 1.0e-10, "Methane at T_ref, P_ref: H_total == 0.0 J/mol");
        TEST_ASSERT(std::abs(*r_ref.S_total) <= 1.0e-10, "Methane at T_ref, P_ref: S_total == 0.0 J/(mol*K)");

        // 2. Methane BM-01: 300 K, 5 MPa
        const auto pvt_bm1 = eos_ch4.calculate_pvt(300.0, 5.0e6);
        const auto& r_bm1 = pvt_bm1.roots[0];
        TEST_ASSERT(r_bm1.H_total.has_value() && r_bm1.S_total.has_value(), "Methane BM-01 has total H/S");
        TEST_ASSERT(approx_equal(*r_bm1.H_total, -818.1829, 0.1), "Methane BM-01: Total H ~ -818.18 J/mol");
        TEST_ASSERT(approx_equal(*r_bm1.S_total, -34.4078, 0.005), "Methane BM-01: Total S ~ -34.408 J/(mol*K)");

        // 3. CO2 BM-02: 300 K, 4 MPa
        const auto comp_co2 = db.get_component("CO2");
        const thermo::PengRobinsonEOS eos_co2(comp_co2);
        const auto pvt_co2_ref = eos_co2.calculate_pvt(thermo::constants::T_REF, thermo::constants::P_REF);
        TEST_ASSERT(std::abs(*pvt_co2_ref.roots[0].H_total) <= 1.0e-10, "CO2 at T_ref, P_ref: H_total == 0.0 J/mol");
        TEST_ASSERT(std::abs(*pvt_co2_ref.roots[0].S_total) <= 1.0e-10, "CO2 at T_ref, P_ref: S_total == 0.0 J/(mol*K)");

        const auto pvt_co2_4 = eos_co2.calculate_pvt(300.0, 4.0e6);
        TEST_ASSERT(approx_equal(*pvt_co2_4.roots[0].H_total, -1881.0417, 0.1), "CO2 BM-02: Total H ~ -1881.04 J/mol");
        TEST_ASSERT(approx_equal(*pvt_co2_4.roots[0].S_total, -35.0832, 0.005), "CO2 BM-02: Total S ~ -35.083 J/(mol*K)");

        // 4. Propane BM-04: 320 K, 1 MPa
        const auto comp_c3h8 = db.get_component("Propane");
        const thermo::PengRobinsonEOS eos_c3h8(comp_c3h8);
        const auto pvt_c3h8 = eos_c3h8.calculate_pvt(320.0, 1.0e6);
        const auto* vap = pvt_c3h8.vapor_root();
        const auto* liq = pvt_c3h8.liquid_root();

        TEST_ASSERT(vap->H_total.has_value() && vap->S_total.has_value(), "Propane Vapor root has total H/S");
        TEST_ASSERT(approx_equal(*vap->H_total, 829.0913, 0.1), "Propane BM-04 Vapor: Total H ~ +829.09 J/mol");
        TEST_ASSERT(approx_equal(*vap->S_total, -15.2835, 0.005), "Propane BM-04 Vapor: Total S ~ -15.283 J/(mol*K)");

        TEST_ASSERT(liq->H_total.has_value() && liq->S_total.has_value(), "Propane Liquid root has total H/S");
        TEST_ASSERT(approx_equal(*liq->H_total, -12894.2087, 0.2), "Propane BM-04 Liquid: Total H ~ -12894.21 J/mol");
        TEST_ASSERT(approx_equal(*liq->S_total, -61.1295, 0.01), "Propane BM-04 Liquid: Total S ~ -61.13 J/(mol*K)");

        // 5. Directional sanity checks:
        // Increasing T at constant P increases total H (positive heat capacity)
        const auto pvt_ch4_warm = eos_ch4.calculate_pvt(350.0, 100000.0);
        const auto pvt_ch4_cool = eos_ch4.calculate_pvt(300.0, 100000.0);
        TEST_ASSERT(*pvt_ch4_warm.roots[0].H_total > *pvt_ch4_cool.roots[0].H_total,
                    "Total H increases with temperature at constant pressure");

        // Increasing P at constant T decreases entropy (ideal gas pressure contribution -R*ln(P/P_ref))
        const auto pvt_ch4_high_p = eos_ch4.calculate_pvt(300.0, 200000.0);
        TEST_ASSERT(*pvt_ch4_high_p.roots[0].S_total < *pvt_ch4_cool.roots[0].S_total,
                    "Total S decreases with pressure at constant temperature");

        // 6. Water Vapor limitation: gas-phase Shomate is [500, 1700] K
        // T_ref = 298.15 K is outside range; total H/S must throw std::out_of_range
        const auto comp_h2o = db.get_component("Water");
        const thermo::PengRobinsonEOS eos_h2o(comp_h2o);
        bool threw_h2o = false;
        try {
            static_cast<void>(eos_h2o.calculate_H_total(600.0, 100000.0, 0.95));
        } catch (const std::out_of_range& e) {
            threw_h2o = true;
        }
        TEST_ASSERT(threw_h2o, "Water Vapor calculate_H_total throws out_of_range because T_ref < 500 K");

        // However, forward PVT on Water at 600 K, 1 bar works, with H_total left as nullopt
        const auto pvt_h2o = eos_h2o.calculate_pvt(600.0, 100000.0);
        TEST_ASSERT(pvt_h2o.roots.size() == 1, "Water at 600 K, 1 bar has 1 physical root");
        TEST_ASSERT(!pvt_h2o.roots[0].H_total.has_value(), "Water root has nullopt for H_total");
        TEST_ASSERT(std::isfinite(pvt_h2o.roots[0].H_res), "Water root still evaluates H_res");
    }

    // ========================================================================
    // Part E: Thermodynamic Phase Preference & Saturation Equilibrium
    // ========================================================================
    std::cout << "\n--- Part E: Phase Preference & Equilibrium Analysis ---\n";
    {
        const auto comp_c3h8 = db.get_component("Propane");
        const thermo::PengRobinsonEOS eos_c3h8(comp_c3h8);

        // Case 1: Propane at 320 K, 1.0 MPa (P < P_sat)
        // Expected: Vapor phase preferred (phi_vap < phi_liq, Delta G_res > 0)
        const auto pvt_1bar = eos_c3h8.calculate_pvt(320.0, 1.0e6);
        TEST_ASSERT(pvt_1bar.roots.size() == 3, "Propane at 1.0 MPa has 3 roots");
        TEST_ASSERT(pvt_1bar.equilibrium.preference == thermo::PhasePreference::VAPOR_PREFERRED,
                    "Propane at 1.0 MPa: VAPOR_PREFERRED");
        TEST_ASSERT(pvt_1bar.equilibrium.delta_G_res > 0.0,
                    "Propane at 1.0 MPa: Delta G_res = RT*ln(phi_l / phi_v) > 0");
        TEST_ASSERT(approx_equal(pvt_1bar.equilibrium.delta_G_res, 947.01, 1.0),
                    "Propane at 1.0 MPa: Delta G_res ~ +947.0 J/mol");
        TEST_ASSERT(pvt_1bar.preferred_root() == pvt_1bar.vapor_root(),
                    "preferred_root() points to vapor_root()");

        // Case 2: Propane at 320 K, 1.8 MPa (P > P_sat)
        // Expected: Liquid phase preferred (phi_liq < phi_vap, Delta G_res < 0)
        const auto pvt_18bar = eos_c3h8.calculate_pvt(320.0, 1.8e6);
        TEST_ASSERT(pvt_18bar.roots.size() == 3, "Propane at 1.8 MPa has 3 roots");
        TEST_ASSERT(pvt_18bar.equilibrium.preference == thermo::PhasePreference::LIQUID_PREFERRED,
                    "Propane at 1.8 MPa: LIQUID_PREFERRED");
        TEST_ASSERT(pvt_18bar.equilibrium.delta_G_res < 0.0,
                    "Propane at 1.8 MPa: Delta G_res < 0");
        TEST_ASSERT(approx_equal(pvt_18bar.equilibrium.delta_G_res, -200.91, 1.0),
                    "Propane at 1.8 MPa: Delta G_res ~ -200.9 J/mol");
        TEST_ASSERT(pvt_18bar.preferred_root() == pvt_18bar.liquid_root(),
                    "preferred_root() points to liquid_root()");

        // Case 3: Propane at 320 K, 1.603614 MPa (P ~= P_sat)
        // Expected: Phase equilibrium (phi_liq ~= phi_vap within tolerance)
        const double P_sat = 1603614.23;
        const auto pvt_sat = eos_c3h8.calculate_pvt(320.0, P_sat);
        TEST_ASSERT(pvt_sat.roots.size() == 3, "Propane at P_sat has 3 roots");
        TEST_ASSERT(pvt_sat.equilibrium.preference == thermo::PhasePreference::EQUILIBRIUM,
                    "Propane at P_sat: EQUILIBRIUM preference");
        TEST_ASSERT(pvt_sat.equilibrium.is_equilibrium, "is_equilibrium flag is true");
        TEST_ASSERT(approx_rel_equal(pvt_sat.equilibrium.phi_liq, pvt_sat.equilibrium.phi_vap, 1.0e-5),
                    "phi_liq == phi_vap within 1e-5 relative tolerance");
        TEST_ASSERT(std::abs(pvt_sat.equilibrium.delta_G_res) < 0.1,
                    "Delta G_res is effectively 0 at saturation");

        // Case 4: Single physical root states (Methane 300K, 5MPa and CO2 300K, 7MPa)
        const auto comp_ch4 = db.get_component("Methane");
        const thermo::PengRobinsonEOS eos_ch4(comp_ch4);
        const auto pvt_single = eos_ch4.calculate_pvt(300.0, 5.0e6);
        TEST_ASSERT(pvt_single.roots.size() == 1, "Methane at 300K, 5MPa is single root");
        TEST_ASSERT(pvt_single.equilibrium.preference == thermo::PhasePreference::SINGLE_PHASE,
                    "Single root reports SINGLE_PHASE preference (no false coexistence)");
        TEST_ASSERT(!pvt_single.equilibrium.is_equilibrium, "Single root is not reported as equilibrium");
    }

    // ========================================================================
    // Part F: Robustness & Domain Validation
    // ========================================================================
    std::cout << "\n--- Part F: Robustness & Domain Guardrails ---\n";
    {
        const auto comp_ch4 = db.get_component("Methane");
        const thermo::PengRobinsonEOS eos_ch4(comp_ch4);

        // 1. Z <= B domain rejection in individual property functions
        const double T = 300.0;
        const double P = 1.0e5;
        const double B = eos_ch4.B_parameter(T, P);

        bool threw_h = false;
        try {
            static_cast<void>(eos_ch4.calculate_H_res(T, P, B)); // Z == B
        } catch (const std::invalid_argument& e) {
            threw_h = true;
        }
        TEST_ASSERT(threw_h, "calculate_H_res rejects Z <= B");

        bool threw_s = false;
        try {
            static_cast<void>(eos_ch4.calculate_S_res(T, P, B * 0.5)); // Z < B
        } catch (const std::invalid_argument& e) {
            threw_s = true;
        }
        TEST_ASSERT(threw_s, "calculate_S_res rejects Z < B");

        bool threw_phi = false;
        try {
            static_cast<void>(eos_ch4.calculate_ln_phi(T, P, B));
        } catch (const std::invalid_argument& e) {
            threw_phi = true;
        }
        TEST_ASSERT(threw_phi, "calculate_ln_phi rejects Z <= B");

        // 2. Non-positive temperature or pressure
        bool threw_t = false;
        try {
            static_cast<void>(eos_ch4.calculate_H_res(-10.0, P, 1.0));
        } catch (const std::invalid_argument& e) {
            threw_t = true;
        }
        TEST_ASSERT(threw_t, "calculate_H_res rejects T <= 0");

        bool threw_p = false;
        try {
            static_cast<void>(eos_ch4.calculate_ln_phi(T, -100.0, 1.0));
        } catch (const std::invalid_argument& e) {
            threw_p = true;
        }
        TEST_ASSERT(threw_p, "calculate_ln_phi rejects P <= 0");

        // 3. NaN inputs
        bool threw_nan = false;
        try {
            static_cast<void>(eos_ch4.calculate_H_res(std::numeric_limits<double>::quiet_NaN(), P, 1.0));
        } catch (const std::invalid_argument& e) {
            threw_nan = true;
        }
        TEST_ASSERT(threw_nan, "calculate_H_res rejects NaN temperature");

        // 4. Shomate temperature bounds out of range
        bool threw_shomate_hi = false;
        try {
            static_cast<void>(eos_ch4.calculate_H_total(2000.0, P, 1.0)); // CH4 upper bound is 1300 K
        } catch (const std::out_of_range& e) {
            threw_shomate_hi = true;
        }
        TEST_ASSERT(threw_shomate_hi, "calculate_H_total rejects T > T_max (1300 K for Methane)");
    }

    std::cout << "\n========================================================\n";
    std::cout << "SUMMARY: Passed " << passed_assertions << " / " << total_assertions << " assertions.\n";
    std::cout << "========================================================\n";

    return (passed_assertions == total_assertions) ? 0 : 1;
}
