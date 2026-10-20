#include "thermo/Component.hpp"
#include "thermo/ComponentDatabase.hpp"
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

bool approx_equal(double a, double b, double tol = 1.0e-3) {
    return std::abs(a - b) <= tol;
}

} // namespace

int main() {
    std::cout << "==================================================\n";
    std::cout << "RUNNING COMPONENT & SHOMATE UNIT TESTS (Session 2)\n";
    std::cout << "==================================================\n";

    // ------------------------------------------------------------------------
    // Test 1: Component validation checks
    // ------------------------------------------------------------------------
    {
        thermo::PureComponent valid_comp{
            "TestGas", "TG", 0.028, 150.0, 3.5e6, 0.05, std::nullopt, "Test citation"
        };
        report_test("Test 1: Valid component passes is_valid()", valid_comp.is_valid());

        bool threw = false;
        try {
            valid_comp.validate();
        } catch (...) {
            threw = true;
        }
        report_test("Test 1: Valid component does not throw", !threw);

        // Invalid empty name
        thermo::PureComponent bad_name = valid_comp;
        bad_name.name = "";
        threw = false;
        try { bad_name.validate(); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test 1: Empty name throws invalid_argument", threw);

        // Invalid non-positive Tc
        thermo::PureComponent bad_tc = valid_comp;
        bad_tc.Tc = -10.0;
        threw = false;
        try { bad_tc.validate(); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test 1: Negative Tc throws invalid_argument", threw);

        // Invalid non-positive Pc
        thermo::PureComponent bad_pc = valid_comp;
        bad_pc.Pc = 0.0;
        threw = false;
        try { bad_pc.validate(); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test 1: Zero Pc throws invalid_argument", threw);

        // Invalid acentric factor <= -1.0
        thermo::PureComponent bad_omega = valid_comp;
        bad_omega.omega = -1.5;
        threw = false;
        try { bad_omega.validate(); } catch (const std::invalid_argument&) { threw = true; }
        report_test("Test 1: Acentric factor < -1.0 throws invalid_argument", threw);
    }

    // ------------------------------------------------------------------------
    // Test 2: Curated Component Database inspection
    // ------------------------------------------------------------------------
    thermo::ComponentDatabase db;
    {
        report_test("Test 2: Database contains exactly 8 components", db.size() == 8);
        const auto names = db.list_component_names();
        report_test("Test 2: List names returns 8 entries", names.size() == 8);

        // Verify all 8 gases exist
        report_test("Test 2: Has Methane", db.has_component("Methane") && db.has_component("CH4"));
        report_test("Test 2: Has Ethane", db.has_component("Ethane") && db.has_component("C2H6"));
        report_test("Test 2: Has Propane", db.has_component("Propane") && db.has_component("C3H8"));
        report_test("Test 2: Has n-Butane", db.has_component("n-Butane") && db.has_component("C4H10"));
        report_test("Test 2: Has Nitrogen", db.has_component("Nitrogen") && db.has_component("N2"));
        report_test("Test 2: Has Oxygen", db.has_component("Oxygen") && db.has_component("O2"));
        report_test("Test 2: Has Carbon Dioxide", db.has_component("Carbon Dioxide") && db.has_component("CO2"));
        report_test("Test 2: Has Water", db.has_component("Water") && db.has_component("H2O"));

        // Case-insensitive lookups
        report_test("Test 2: Case-insensitive 'ch4'", db.has_component("ch4"));
        report_test("Test 2: Case-insensitive 'co2'", db.has_component("co2"));
        report_test("Test 2: Case-insensitive 'carbon dioxide'", db.has_component("carbon dioxide"));

        // Non-existent component
        report_test("Test 2: Non-existent gas returns false", !db.has_component("Krypton"));
        bool threw = false;
        try {
            (void)db.get_component("Krypton");
        } catch (const std::runtime_error&) {
            threw = true;
        }
        report_test("Test 2: Non-existent gas throws runtime_error", threw);
    }

    // ------------------------------------------------------------------------
    // Test 3: Critical Properties accuracy against NIST SRD 69
    // ------------------------------------------------------------------------
    {
        const auto& ch4 = db.get_component("CH4");
        report_test("Test 3: Methane Tc == 190.56 K", approx_equal(ch4.Tc, 190.56));
        report_test("Test 3: Methane Pc == 4.599 MPa", approx_equal(ch4.Pc, 4.599e6));
        report_test("Test 3: Methane omega == 0.0115", approx_equal(ch4.omega, 0.0115));

        const auto& co2 = db.get_component("CO2");
        report_test("Test 3: CO2 Tc == 304.13 K", approx_equal(co2.Tc, 304.13));
        report_test("Test 3: CO2 Pc == 7.377 MPa", approx_equal(co2.Pc, 7.377e6));
        report_test("Test 3: CO2 omega == 0.2239", approx_equal(co2.omega, 0.2239));
    }

    // ------------------------------------------------------------------------
    // Test 4: NIST Shomate Cp evaluations at 298.15 K
    // ------------------------------------------------------------------------
    {
        const auto& ch4 = db.get_component("CH4");
        assert(ch4.shomate.has_value());
        const double cp_ch4 = ch4.shomate->evaluate_cp(298.15);
        report_test("Test 4: Methane Cp(298.15) == 35.65 J/(mol*K) [NIST: 35.69]",
                    approx_equal(cp_ch4, 35.648, 0.05));

        const auto& co2 = db.get_component("CO2");
        assert(co2.shomate.has_value());
        const double cp_co2 = co2.shomate->evaluate_cp(298.15);
        report_test("Test 4: CO2 Cp(298.15) == 37.13 J/(mol*K) [NIST: 37.13]",
                    approx_equal(cp_co2, 37.130, 0.05));

        const auto& n2 = db.get_component("N2");
        assert(n2.shomate.has_value());
        const double cp_n2 = n2.shomate->evaluate_cp(298.15);
        report_test("Test 4: N2 Cp(298.15) == 29.12 J/(mol*K) [NIST: 29.12]",
                    approx_equal(cp_n2, 29.124, 0.05));
    }

    // ------------------------------------------------------------------------
    // Test 5: Analytical definite integrals of NIST Shomate
    // ------------------------------------------------------------------------
    {
        const auto& co2 = db.get_component("CO2");
        const auto& sh = co2.shomate.value();

        // At T == T_ref, integral must be identically zero
        const double delta_H_0 = sh.integrate_cp(298.15, 298.15);
        const double delta_S_0 = sh.integrate_cp_over_T(298.15, 298.15);
        report_test("Test 5: delta_H_ig(T_ref -> T_ref) == 0.0", std::abs(delta_H_0) < 1.0e-9);
        report_test("Test 5: delta_S_ig(T_ref -> T_ref) == 0.0", std::abs(delta_S_0) < 1.0e-9);

        // Integrals to 400.0 K (audited values: delta_H = 4005.43 J/mol, delta_S = 11.5267 J/(mol*K))
        const double delta_H_400 = sh.integrate_cp(298.15, 400.0);
        const double delta_S_400 = sh.integrate_cp_over_T(298.15, 400.0);
        report_test("Test 5: CO2 delta_H_ig(298.15 -> 400 K) == 4005.43 J/mol",
                    approx_equal(delta_H_400, 4005.43, 0.5));
        report_test("Test 5: CO2 delta_S_ig(298.15 -> 400 K) == 11.5267 J/(mol*K)",
                    approx_equal(delta_S_400, 11.5267, 0.01));
    }

    // ------------------------------------------------------------------------
    // Test 6: Temperature bounds enforcement & rejection of extrapolation
    // ------------------------------------------------------------------------
    {
        const auto& co2 = db.get_component("CO2");
        const auto& sh = co2.shomate.value(); // Range [298.0, 1200.0 K]

        // Below T_min
        bool threw = false;
        try {
            (void)sh.evaluate_cp(200.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }
        report_test("Test 6: CO2 Cp at 200 K (< T_min) throws out_of_range", threw);

        // Above T_max
        threw = false;
        try {
            (void)sh.evaluate_cp(1300.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }
        report_test("Test 6: CO2 Cp at 1300 K (> T_max) throws out_of_range", threw);

        // Water vapor range [500.0, 1700.0 K] -> 298.15 K is out of range
        const auto& h2o = db.get_component("H2O");
        threw = false;
        try {
            (void)h2o.shomate->evaluate_cp(298.15);
        } catch (const std::out_of_range&) {
            threw = true;
        }
        report_test("Test 6: Water vapor at 298.15 K (< 500 K) throws out_of_range", threw);
    }

    // ------------------------------------------------------------------------
    // Test 7: Custom component registration
    // ------------------------------------------------------------------------
    {
        thermo::PureComponent argon{
            "Argon", "Ar", 0.039948, 150.86, 4.898e6, -0.002, std::nullopt, "NIST Chemistry WebBook"
        };
        db.register_custom_component(argon);
        report_test("Test 7: Database size increases to 9 after custom registration", db.size() == 9);
        report_test("Test 7: Can query custom component 'Argon'", db.has_component("Argon"));
        report_test("Test 7: Can query custom component 'Ar'", db.has_component("Ar"));

        const auto& retrieved = db.get_component("Ar");
        report_test("Test 7: Retrieved custom component Tc == 150.86 K", approx_equal(retrieved.Tc, 150.86));

        // Attempting to register duplicate throws invalid_argument
        bool threw = false;
        try {
            db.register_custom_component(argon);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        report_test("Test 7: Registering duplicate component throws invalid_argument", threw);
    }

    std::cout << "==================================================\n";
    std::cout << "SUMMARY: Passed " << g_pass_count << ", Failed " << g_fail_count << "\n";
    std::cout << "==================================================\n";

    return (g_fail_count == 0) ? 0 : 1;
}
