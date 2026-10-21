#include "cli/CLIController.hpp"
#include <iostream>
#include <sstream>
#include <cassert>
#include <string>
#include <cstdlib>

// Test macro with clear console feedback
#define CLI_ASSERT(cond, msg) \
    do { \
        if (!(cond)) { \
            std::cerr << "[FAIL] Line " << __LINE__ << ": " << msg << "\n"; \
            std::exit(1); \
        } else { \
            std::cout << "[PASS] " << msg << "\n"; \
        } \
    } while (0)

int main() {
    std::cout << "========================================================\n";
    std::cout << "RUNNING CLI AND USER EXPERIENCE TESTS (SESSION 6)\n";
    std::cout << "========================================================\n\n";

    // -------------------------------------------------------------------------
    // Test 1: Valid Component and PVT Calculation via Stream
    // -------------------------------------------------------------------------
    std::cout << "--- Test 1: Valid Component & PVT Calculation ---\n";
    {
        // 1: PVT
        // 1: Methane
        // 300: 300
        // 1: Kelvin
        // 5.0: 5.0
        // 3: MPa
        // 5: Exit
        std::istringstream input("1\n1\n300\n1\n5.0\n3\n5\n");
        std::ostringstream output;
        thermo::cli::CLIController controller(input, output);

        int exit_code = controller.run();
        CLI_ASSERT(exit_code == 0, "CLI returns 0 on clean menu exit");

        std::string text = output.str();
        CLI_ASSERT(text.find("Methane") != std::string::npos, "Output contains component name Methane");
        CLI_ASSERT(text.find("0.901844") != std::string::npos, "Output contains correct Z = 0.901844");
        CLI_ASSERT(text.find("4.4990e-04 m^3/mol") != std::string::npos, "Output contains molar volume");
        CLI_ASSERT(text.find("STABLE") != std::string::npos, "Output reports STABLE mechanical stability");
        CLI_ASSERT(text.find("SINGLE_PHASE") != std::string::npos, "Output reports SINGLE_PHASE fluid");
        CLI_ASSERT(text.find("Thank you for using CHL2002") != std::string::npos, "Output contains exit message");
    }

    // -------------------------------------------------------------------------
    // Test 2: Invalid Component Selection and Recovery
    // -------------------------------------------------------------------------
    std::cout << "\n--- Test 2: Invalid Component Selection & Recovery ---\n";
    {
        // 1: PVT
        // 99: Out of range
        // InvalidGas: Unknown name
        // 1: Methane (valid)
        // 300: 300 K
        // 1: Kelvin
        // 5: 5 MPa
        // 3: MPa
        // 5: Exit
        std::istringstream input("1\n99\nInvalidGas\n1\n300\n1\n5\n3\n5\n");
        std::ostringstream output;
        thermo::cli::CLIController controller(input, output);

        int exit_code = controller.run();
        CLI_ASSERT(exit_code == 0, "CLI recovers from invalid component selection and exits cleanly");

        std::string text = output.str();
        CLI_ASSERT(text.find("Error: Number 99 is out of range") != std::string::npos,
                   "Reports error for out-of-range integer component choice");
        CLI_ASSERT(text.find("Error: Component 'InvalidGas' not recognized") != std::string::npos,
                   "Reports error for unrecognized string component name");
        CLI_ASSERT(text.find("0.901844") != std::string::npos,
                   "Successfully completed PVT calculation after input recovery");
    }

    // -------------------------------------------------------------------------
    // Test 3: Invalid and Non-Finite Numeric Input Handling
    // -------------------------------------------------------------------------
    std::cout << "\n--- Test 3: Invalid & Non-Finite Numeric Input Handling ---\n";
    {
        // 1: PVT
        // 1: Methane
        // bad_num: string instead of number for T
        // -50: negative T
        // 300: valid T
        // 1: Kelvin
        // 0: zero P
        // 5.0: valid P
        // 3: MPa
        // 5: Exit
        std::istringstream input("1\n1\nbad_num\n-50\n300\n1\n0\n5.0\n3\n5\n");
        std::ostringstream output;
        thermo::cli::CLIController controller(input, output);

        int exit_code = controller.run();
        CLI_ASSERT(exit_code == 0, "CLI recovers from invalid numeric values without throwing or crashing");

        std::string text = output.str();
        CLI_ASSERT(text.find("Error: Invalid numeric input 'bad_num'") != std::string::npos,
                   "Rejects non-numeric characters for double input");
        CLI_ASSERT(text.find("Error: Value must be strictly positive") != std::string::npos,
                   "Rejects negative and zero values for physical temperature and pressure");
    }

    // -------------------------------------------------------------------------
    // Test 4: Graceful Handling of Unsupported Shomate Range (Water Vapor)
    // -------------------------------------------------------------------------
    std::cout << "\n--- Test 4: Water Vapor Reference-State Limitation Handling ---\n";
    {
        // 2: Thermodynamic properties (H & S)
        // 8: Water
        // 600: 600 K (within 500-1700 K, but T_ref 298.15 K is not)
        // 1: Kelvin
        // 1.0: 1.0 bar
        // 2: bar
        // 5: Exit
        std::istringstream input("2\n8\n600\n1\n1.0\n2\n5\n");
        std::ostringstream output;
        thermo::cli::CLIController controller(input, output);

        int exit_code = controller.run();
        CLI_ASSERT(exit_code == 0, "Handles Water Vapor total H/S limitation without crashing");

        std::string text = output.str();
        CLI_ASSERT(text.find("Residual Enthalpy:") != std::string::npos,
                   "Residual enthalpy is evaluated for Water Vapor");
        CLI_ASSERT(text.find("Total Enthalpy (H):  N/A (Reference state 298.15 K is below gas-phase Shomate range") != std::string::npos,
                   "Explicitly explains why total H is N/A instead of printing dummy zero");
        CLI_ASSERT(text.find("Total Entropy (S):   N/A (Reference state 298.15 K is below gas-phase Shomate range") != std::string::npos,
                   "Explicitly explains why total S is N/A instead of printing dummy zero");
    }

    // -------------------------------------------------------------------------
    // Test 5: Three-Root State and Phase Preference (Propane 320 K, 1 MPa)
    // -------------------------------------------------------------------------
    std::cout << "\n--- Test 5: Three-Root State & Phase Preference (Propane 320 K, 1 MPa) ---\n";
    {
        // 4: Comprehensive All Properties
        // 3: Propane
        // 320: 320 K
        // 1: Kelvin
        // 1.0: 1.0 MPa
        // 3: MPa
        // 5: Exit
        std::istringstream input("4\n3\n320\n1\n1.0\n3\n5\n");
        std::ostringstream output;
        thermo::cli::CLIController controller(input, output);

        int exit_code = controller.run();
        CLI_ASSERT(exit_code == 0, "Evaluates 3-root state cleanly");

        std::string text = output.str();
        CLI_ASSERT(text.find("Physical Roots Discovered: 3 physical root(s)") != std::string::npos,
                   "Identifies all 3 physical roots for Propane at 320 K, 1 MPa");
        CLI_ASSERT(text.find("LIQUID-LIKE") != std::string::npos, "Identifies liquid-like root");
        CLI_ASSERT(text.find("INTERMEDIATE (UNPHYSICAL)") != std::string::npos,
                   "Identifies intermediate root");
        CLI_ASSERT(text.find("UNSTABLE: (dP/dv)_T >= 0") != std::string::npos,
                   "Flags intermediate root as mechanically unstable");
        CLI_ASSERT(text.find("VAPOR-LIKE") != std::string::npos, "Identifies vapor-like root");
        CLI_ASSERT(text.find("VAPOR_PREFERRED") != std::string::npos,
                   "Correctly determines vapor phase preference");
        CLI_ASSERT(text.find("Delta G_res = G_liq - G_vap = +947") != std::string::npos,
                   "Calculates positive Delta G_res confirming vapor stability");
        CLI_ASSERT(text.find("Thermodynamically Selected Root: Z = 0.850555") != std::string::npos,
                   "Reports Root 3 (vapor) as thermodynamically selected root");
    }

    // -------------------------------------------------------------------------
    // Test 6: Pure Fugacity Calculation Menu Option (Option 3)
    // -------------------------------------------------------------------------
    std::cout << "\n--- Test 6: Pure Fugacity Calculation Menu Option ---\n";
    {
        // 3: Fugacity & Fugacity Coefficient
        // 7: Carbon Dioxide (CO2)
        // 300: 300 K
        // 1: Kelvin
        // 4.0: 4.0 MPa (BM-02)
        // 3: MPa
        // 5: Exit
        std::istringstream input("3\n7\n300\n1\n4.0\n3\n5\n");
        std::ostringstream output;
        thermo::cli::CLIController controller(input, output);

        int exit_code = controller.run();
        CLI_ASSERT(exit_code == 0, "Option 3 executes cleanly");

        std::string text = output.str();
        CLI_ASSERT(text.find("PURE-COMPONENT FUGACITY RESULTS") != std::string::npos,
                   "Displays pure fugacity header");
        CLI_ASSERT(text.find("0.795528") != std::string::npos,
                   "Displays correct fugacity coefficient phi = 0.795528 for CO2 at 4 MPa");
        CLI_ASSERT(text.find("3,182,110") != std::string::npos,
                   "Displays fugacity in Pa with comma formatting");
        CLI_ASSERT(text.find("3.182110 MPa") != std::string::npos,
                   "Displays fugacity in readable MPa");
    }

    // -------------------------------------------------------------------------
    // Test 7: Clean Exit and Immediate EOF Handling
    // -------------------------------------------------------------------------
    std::cout << "\n--- Test 7: Clean Exit & Immediate EOF Handling ---\n";
    {
        std::istringstream empty_input("");
        std::ostringstream output;
        thermo::cli::CLIController controller(empty_input, output);

        int exit_code = controller.run();
        CLI_ASSERT(exit_code == 0, "CLI handles immediate EOF gracefully without looping or crashing");

        std::string text = output.str();
        CLI_ASSERT(text.find("End of input reached. Exiting gracefully.") != std::string::npos,
                   "Prints graceful exit message on EOF");
    }

    // -------------------------------------------------------------------------
    // Test 8: Predefined Demonstrations Direct Execution
    // -------------------------------------------------------------------------
    std::cout << "\n--- Test 8: Predefined Demonstrations Direct Execution ---\n";
    {
        std::istringstream dummy_in("");
        std::ostringstream output;
        thermo::cli::CLIController controller(dummy_in, output);

        controller.run_demonstration(1); // Methane
        controller.run_demonstration(2); // CO2
        controller.run_demonstration(3); // Propane

        std::string text = output.str();
        CLI_ASSERT(text.find("DEMONSTRATION CASE 1: METHANE") != std::string::npos,
                   "Demonstration 1 executes successfully");
        CLI_ASSERT(text.find("DEMONSTRATION CASE 2: CARBON DIOXIDE") != std::string::npos,
                   "Demonstration 2 executes successfully");
        CLI_ASSERT(text.find("DEMONSTRATION CASE 3: PROPANE") != std::string::npos,
                   "Demonstration 3 executes successfully");
    }

    std::cout << "\n========================================================\n";
    std::cout << "ALL CLI TESTS COMPLETED SUCCESSFULLY!\n";
    std::cout << "========================================================\n";
    return 0;
}
