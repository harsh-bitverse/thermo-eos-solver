#include "cli/CLIController.hpp"
#include "thermo/Constants.hpp"
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cctype>
#include <algorithm>

namespace thermo::cli {

namespace {

// Trims leading and trailing whitespace from a string
std::string trim(const std::string& str) {
    const auto first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

constexpr double CELSIUS_TO_KELVIN = 273.15;
constexpr double BAR_TO_PA = 100000.0;


// Helper to format large integers with commas (e.g. 1,000,000)
std::string format_with_commas(double val, int precision = 0) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << val;
    std::string s = ss.str();
    size_t dot_pos = s.find('.');
    std::string int_part = (dot_pos == std::string::npos) ? s : s.substr(0, dot_pos);
    std::string frac_part = (dot_pos == std::string::npos) ? "" : s.substr(dot_pos);

    bool negative = false;
    if (!int_part.empty() && int_part[0] == '-') {
        negative = true;
        int_part = int_part.substr(1);
    }

    int n = static_cast<int>(int_part.length()) - 3;
    while (n > 0) {
        int_part.insert(static_cast<size_t>(n), ",");
        n -= 3;
    }
    return (negative ? "-" : "") + int_part + frac_part;
}

const char* phase_tag_str(PhaseTag tag) noexcept {
    switch (tag) {
        case PhaseTag::SINGLE_PHASE: return "SINGLE-PHASE FLUID";
        case PhaseTag::VAPOR:        return "VAPOR-LIKE";
        case PhaseTag::LIQUID:       return "LIQUID-LIKE";
        case PhaseTag::INTERMEDIATE: return "INTERMEDIATE (UNPHYSICAL)";
    }
    return "UNKNOWN";
}

const char* stability_str(MechanicalStability stab) noexcept {
    switch (stab) {
        case MechanicalStability::STABLE:   return "STABLE";
        case MechanicalStability::UNSTABLE: return "UNSTABLE: (dP/dv)_T >= 0";
    }
    return "UNKNOWN";
}

} // anonymous namespace

CLIController::CLIController(std::istream& in, std::ostream& out)
    : in_(in), out_(out), db_() {}

bool CLIController::read_line(std::string& line) {
    if (!std::getline(in_, line)) {
        return false;
    }
    line = trim(line);
    return true;
}

std::optional<int> CLIController::prompt_int(const std::string& prompt_text,
                                             int min_val, int max_val,
                                             std::optional<int> default_val) {
    while (true) {
        out_ << prompt_text;
        if (default_val.has_value()) {
            out_ << " [" << *default_val << "]: ";
        } else {
            out_ << ": ";
        }

        std::string line;
        if (!read_line(line)) {
            return std::nullopt; // EOF
        }

        if (line.empty()) {
            if (default_val.has_value()) {
                return default_val;
            }
            out_ << "Error: Input cannot be empty. Please enter an integer between "
                 << min_val << " and " << max_val << ".\n";
            continue;
        }

        try {
            size_t idx = 0;
            int val = std::stoi(line, &idx);
            if (idx != line.length()) {
                out_ << "Error: Unexpected trailing characters '" << line.substr(idx)
                     << "'. Please enter a valid integer.\n";
                continue;
            }
            if (val < min_val || val > max_val) {
                out_ << "Error: Value " << val << " is out of allowed range ["
                     << min_val << ", " << max_val << "]. Please try again.\n";
                continue;
            }
            return val;
        } catch (const std::exception&) {
            out_ << "Error: Invalid input '" << line << "'. Please enter an integer between "
                 << min_val << " and " << max_val << ".\n";
        }
    }
}

std::optional<double> CLIController::prompt_positive_double(const std::string& prompt_text,
                                                            std::optional<double> default_val) {
    while (true) {
        out_ << prompt_text;
        if (default_val.has_value()) {
            out_ << " [" << *default_val << "]: ";
        } else {
            out_ << ": ";
        }

        std::string line;
        if (!read_line(line)) {
            return std::nullopt; // EOF
        }

        if (line.empty()) {
            if (default_val.has_value()) {
                return default_val;
            }
            out_ << "Error: Input cannot be empty. Please enter a positive number.\n";
            continue;
        }

        try {
            size_t idx = 0;
            double val = std::stod(line, &idx);
            if (idx != line.length()) {
                out_ << "Error: Unexpected trailing characters '" << line.substr(idx)
                     << "'. Please enter a clean numeric value.\n";
                continue;
            }
            if (std::isnan(val) || std::isinf(val)) {
                out_ << "Error: Value must be finite (not NaN or Infinity). Please try again.\n";
                continue;
            }
            if (val <= 0.0) {
                out_ << "Error: Value must be strictly positive (> 0). Please try again.\n";
                continue;
            }
            return val;
        } catch (const std::exception&) {
            out_ << "Error: Invalid numeric input '" << line << "'. Please enter a positive number.\n";
        }
    }
}

std::optional<PureComponent> CLIController::prompt_component() {
    auto names = db_.list_component_names();

    std::vector<std::string> clean_names;
    clean_names.reserve(names.size());
    for (const auto& item : names) {
        size_t p = item.find(" (");
        clean_names.push_back((p != std::string::npos) ? item.substr(0, p) : item);
    }

    while (true) {
        out_ << "\nAvailable Curated Pure Components:\n";
        for (size_t i = 0; i < clean_names.size(); ++i) {
            const auto& comp = db_.get_component(clean_names[i]);
            out_ << "  [" << (i + 1) << "] " << std::left << std::setw(16) << comp.name
                 << " (" << std::left << std::setw(6) << comp.formula << ") "
                 << "Tc = " << std::setw(7) << std::fixed << std::setprecision(2) << comp.Tc << " K, "
                 << "Pc = " << std::setw(6) << std::fixed << std::setprecision(3) << (comp.Pc / 1.0e6) << " MPa\n";
        }
        out_ << "Select component by number [1-" << clean_names.size()
             << "] or name/formula (or 'q' to cancel): ";

        std::string line;
        if (!read_line(line)) {
            return std::nullopt;
        }

        if (line.empty()) {
            out_ << "Error: Selection cannot be empty. Please enter 1-" << clean_names.size() << " or name.\n";
            continue;
        }

        if (line == "q" || line == "Q" || line == "cancel") {
            out_ << "Operation cancelled by user.\n";
            return std::nullopt;
        }

        // Check if numeric choice [1..N]
        bool is_num = true;
        for (char c : line) {
            if (!std::isdigit(static_cast<unsigned char>(c))) {
                is_num = false;
                break;
            }
        }

        if (is_num) {
            try {
                int choice = std::stoi(line);
                if (choice >= 1 && choice <= static_cast<int>(clean_names.size())) {
                    return db_.get_component(clean_names[static_cast<size_t>(choice - 1)]);
                }
                out_ << "Error: Number " << choice << " is out of range [1-" << clean_names.size() << "].\n";
                continue;
            } catch (...) {
                // fall through to text match
            }
        }

        // Match by name or formula
        if (db_.has_component(line)) {
            return db_.get_component(line);
        }

        out_ << "Error: Component '" << line << "' not recognized. Please choose from the list.\n";
    }
}

std::optional<double> CLIController::prompt_temperature() {
    out_ << "\n--- Temperature Input ---\n";
    auto val_opt = prompt_positive_double("Enter temperature value (e.g. 300)");
    if (!val_opt.has_value()) return std::nullopt;
    double raw_val = *val_opt;

    out_ << "Select temperature unit:\n"
         << "  [1] Kelvin [K] (SI default)\n"
         << "  [2] Celsius [°C]\n";
    auto unit_opt = prompt_int("Choice", 1, 2, 1);
    if (!unit_opt.has_value()) return std::nullopt;

    double T_K = raw_val;
    if (*unit_opt == 2) {
        T_K = raw_val + CELSIUS_TO_KELVIN;
        if (T_K <= 0.0) {
            out_ << "Error: Temperature " << raw_val << " °C (" << T_K
                 << " K) is at or below absolute zero (0 K). Please try again.\n";
            return prompt_temperature();
        }
        out_ << "[Unit Conversion] " << std::fixed << std::setprecision(2) << raw_val
             << " °C = " << std::setprecision(4) << T_K << " K\n";
    } else {
        out_ << "[Temperature] " << std::fixed << std::setprecision(4) << T_K << " K\n";
    }
    return T_K;
}

std::optional<double> CLIController::prompt_pressure() {
    out_ << "\n--- Pressure Input ---\n";
    auto val_opt = prompt_positive_double("Enter pressure value (e.g. 5.0 or 100000)");
    if (!val_opt.has_value()) return std::nullopt;
    double raw_val = *val_opt;

    out_ << "Select pressure unit:\n"
         << "  [1] Pascal [Pa] (SI default)\n"
         << "  [2] Bar [bar] (1 bar = 100,000 Pa)\n"
         << "  [3] Megapascal [MPa] (1 MPa = 1,000,000 Pa)\n";
    auto unit_opt = prompt_int("Choice", 1, 3, 1);
    if (!unit_opt.has_value()) return std::nullopt;

    double P_Pa = raw_val;
    if (*unit_opt == 2) {
        P_Pa = raw_val * BAR_TO_PA;
        out_ << "[Unit Conversion] " << std::fixed << std::setprecision(4) << raw_val
             << " bar = " << format_with_commas(P_Pa, 2) << " Pa ("
             << std::setprecision(4) << (P_Pa / 1.0e6) << " MPa)\n";
    } else if (*unit_opt == 3) {
        P_Pa = raw_val * 1.0e6;
        out_ << "[Unit Conversion] " << std::fixed << std::setprecision(4) << raw_val
             << " MPa = " << format_with_commas(P_Pa, 2) << " Pa ("
             << std::setprecision(4) << (P_Pa / BAR_TO_PA) << " bar)\n";
    } else {
        out_ << "[Pressure] " << format_with_commas(P_Pa, 2) << " Pa ("
             << std::fixed << std::setprecision(4) << (P_Pa / BAR_TO_PA) << " bar = "
             << (P_Pa / 1.0e6) << " MPa)\n";
    }
    return P_Pa;
}

void CLIController::print_header(const std::string& title) const {
    out_ << "\n" << std::string(80, '=') << "\n";
    out_ << "  " << title << "\n";
    out_ << std::string(80, '=') << "\n";
}

void CLIController::print_component_summary(const PureComponent& comp) const {
    out_ << "Component:             " << comp.name << " (" << comp.formula << ")\n"
         << "Molar Mass:            " << std::fixed << std::setprecision(6) << comp.molar_mass
         << " kg/mol (" << std::setprecision(3) << (comp.molar_mass * 1000.0) << " g/mol)\n"
         << "Critical Temperature:  " << std::setprecision(2) << comp.Tc << " K\n"
         << "Critical Pressure:     " << std::setprecision(3) << (comp.Pc / 1.0e6)
         << " MPa (" << format_with_commas(comp.Pc, 0) << " Pa)\n"
         << "Acentric Factor (w):   " << std::setprecision(4) << comp.omega << "\n";
    if (comp.shomate.has_value() && comp.shomate->is_defined()) {
        out_ << "NIST Shomate Range:    [" << std::setprecision(1) << comp.shomate->T_min
             << ", " << comp.shomate->T_max << "] K\n";
    } else {
        out_ << "NIST Shomate Range:    Not defined for this component\n";
    }
}

void CLIController::print_state_conditions(double T, double P, const PureComponent& comp) const {
    const double Tr = T / comp.Tc;
    const double Pr = P / comp.Pc;
    out_ << "--------------------------------------------------------------------------------\n"
         << "State Conditions:\n"
         << "  Temperature (T):     " << std::fixed << std::setprecision(4) << T << " K  ("
         << std::setprecision(2) << (T - CELSIUS_TO_KELVIN) << " °C)  [Tr = "
         << std::setprecision(4) << Tr << "]\n"
         << "  Pressure (P):        " << format_with_commas(P, 2) << " Pa  ("
         << std::setprecision(4) << (P / BAR_TO_PA) << " bar = "
         << (P / 1.0e6) << " MPa)  [Pr = " << std::setprecision(4) << Pr << "]\n";
}

void CLIController::print_eos_parameters(const EOSParameters& params, double A, double B, double discr) const {
    out_ << "--------------------------------------------------------------------------------\n"
         << "Peng-Robinson EOS Parameters:\n"
         << "  Critical Attract. ac:" << std::setw(14) << std::fixed << std::setprecision(6) << params.ac << " Pa*m^6/mol^2\n"
         << "  Attractive a(T):     " << std::setw(14) << std::fixed << std::setprecision(6) << params.a << " Pa*m^6/mol^2\n"
         << "  Co-volume b:         " << std::setw(14) << std::scientific << std::setprecision(5) << params.b << " m^3/mol\n"
         << "  Temperature Deriv:   " << std::setw(14) << std::fixed << std::setprecision(6) << params.da_dT << " Pa*m^6/(mol^2*K)\n"
         << "  Dimensionless A:     " << std::setw(14) << std::fixed << std::setprecision(6) << A << "\n"
         << "  Dimensionless B:     " << std::setw(14) << std::fixed << std::setprecision(6) << B << "\n"
         << "  Cubic Discriminant:  " << std::setw(14) << std::scientific << std::setprecision(5) << discr;
    if (discr < 0.0) {
        out_ << " (< 0: 3 Real Roots)\n";
    } else if (discr > 0.0) {
        out_ << " (> 0: 1 Real Root)\n";
    } else {
        out_ << " (== 0: Real Roots with Coalescence)\n";
    }
}

void CLIController::print_root_pvt(size_t index, const PVTRoot& root, double B) const {
    out_ << "\n[Root " << (index + 1) << " - " << phase_tag_str(root.phase) << "]\n"
         << "  Compressibility (Z): " << std::fixed << std::setprecision(6) << root.Z
         << "  (Physical: Z > B = " << std::setprecision(6) << B << ")\n"
         << "  Molar Volume (v):    " << std::scientific << std::setprecision(4) << root.molar_volume
         << " m^3/mol  (" << std::fixed << std::setprecision(6) << (root.molar_volume * 1000.0) << " L/mol)\n"
         << "  Polynomial Residual: " << std::scientific << std::setprecision(2) << root.poly_residual << "\n"
         << "  Isothermal (dP/dv)_T:" << std::setw(14) << std::scientific << std::setprecision(4) << root.dP_dv
         << " Pa/(m^3/mol)  [" << stability_str(root.stability) << "]\n";
    if (root.stability == MechanicalStability::UNSTABLE) {
        out_ << "  Notice:              Mechanically unstable root within spinodal region; unobservable.\n";
    }
}

void CLIController::print_root_thermo(size_t index, const PVTRoot& root, const PureComponent& comp, double T) const {
    out_ << "\n[Root " << (index + 1) << " - " << phase_tag_str(root.phase) << "]\n"
         << "  Residual Enthalpy:   " << std::setw(12) << std::fixed << std::setprecision(2) << root.H_res << " J/mol\n"
         << "  Residual Entropy:    " << std::setw(12) << std::fixed << std::setprecision(4) << root.S_res << " J/(mol*K)\n";

    if (root.H_total.has_value()) {
        out_ << "  Total Enthalpy (H):  " << std::setw(12) << std::fixed << std::setprecision(2) << *root.H_total
             << " J/mol  [Ref: 298.15 K, 1.0 bar ideal gas]\n";
    } else {
        if (comp.formula == "H2O") {
            out_ << "  Total Enthalpy (H):  N/A (Reference state 298.15 K is below gas-phase Shomate range [500.0, 1700.0] K)\n";
        } else if (!comp.shomate.has_value() || !comp.shomate->is_defined()) {
            out_ << "  Total Enthalpy (H):  N/A (No ideal-gas Shomate heat capacity defined)\n";
        } else {
            out_ << "  Total Enthalpy (H):  N/A (Temperature " << std::fixed << std::setprecision(2) << T
                 << " K is outside valid Shomate range [" << comp.shomate->T_min << ", " << comp.shomate->T_max << "] K)\n";
        }
    }

    if (root.S_total.has_value()) {
        out_ << "  Total Entropy (S):   " << std::setw(12) << std::fixed << std::setprecision(4) << *root.S_total
             << " J/(mol*K) [Ref: 298.15 K, 1.0 bar ideal gas]\n";
    } else {
        if (comp.formula == "H2O") {
            out_ << "  Total Entropy (S):   N/A (Reference state 298.15 K is below gas-phase Shomate range [500.0, 1700.0] K)\n";
        } else if (!comp.shomate.has_value() || !comp.shomate->is_defined()) {
            out_ << "  Total Entropy (S):   N/A (No ideal-gas Shomate heat capacity defined)\n";
        } else {
            out_ << "  Total Entropy (S):   N/A (Temperature " << std::fixed << std::setprecision(2) << T
                 << " K is outside valid Shomate range [" << comp.shomate->T_min << ", " << comp.shomate->T_max << "] K)\n";
        }
    }
}

void CLIController::print_root_fugacity(size_t index, const PVTRoot& root) const {
    out_ << "\n[Root " << (index + 1) << " - " << phase_tag_str(root.phase) << "]\n"
         << "  Fugacity Coeff (phi):" << std::setw(12) << std::fixed << std::setprecision(6) << root.phi
         << "  (ln(phi) = " << std::setprecision(6) << root.ln_phi << ")\n"
         << "  Fugacity (f):        " << format_with_commas(root.fugacity, 2) << " Pa  ("
         << std::fixed << std::setprecision(6) << (root.fugacity / 1.0e6) << " MPa)\n";
}

void CLIController::print_root_all(size_t index, const PVTRoot& root, const PureComponent& comp, double T, double B) const {
    out_ << "\n[Root " << (index + 1) << " - " << phase_tag_str(root.phase) << "]\n"
         << "  Compressibility (Z): " << std::fixed << std::setprecision(6) << root.Z
         << "  (Physical: Z > B = " << std::setprecision(6) << B << ")\n"
         << "  Molar Volume (v):    " << std::scientific << std::setprecision(4) << root.molar_volume
         << " m^3/mol  (" << std::fixed << std::setprecision(6) << (root.molar_volume * 1000.0) << " L/mol)\n"
         << "  Polynomial Residual: " << std::scientific << std::setprecision(2) << root.poly_residual << "\n"
         << "  Isothermal (dP/dv)_T:" << std::setw(14) << std::scientific << std::setprecision(4) << root.dP_dv
         << " Pa/(m^3/mol)  [" << stability_str(root.stability) << "]\n"
         << "  Residual Enthalpy:   " << std::setw(12) << std::fixed << std::setprecision(2) << root.H_res << " J/mol\n"
         << "  Residual Entropy:    " << std::setw(12) << std::fixed << std::setprecision(4) << root.S_res << " J/(mol*K)\n"
         << "  Fugacity Coeff (phi):" << std::setw(12) << std::fixed << std::setprecision(6) << root.phi
         << "  (ln(phi) = " << std::setprecision(6) << root.ln_phi << ")\n"
         << "  Fugacity (f):        " << format_with_commas(root.fugacity, 2) << " Pa  ("
         << std::fixed << std::setprecision(6) << (root.fugacity / 1.0e6) << " MPa)\n";

    if (root.H_total.has_value()) {
        out_ << "  Total Enthalpy (H):  " << std::setw(12) << std::fixed << std::setprecision(2) << *root.H_total
             << " J/mol  [Ref: 298.15 K, 1.0 bar ideal gas]\n";
    } else {
        if (comp.formula == "H2O") {
            out_ << "  Total Enthalpy (H):  N/A (Reference state 298.15 K is below gas-phase Shomate range [500.0, 1700.0] K)\n";
        } else if (!comp.shomate.has_value() || !comp.shomate->is_defined()) {
            out_ << "  Total Enthalpy (H):  N/A (No ideal-gas Shomate heat capacity defined)\n";
        } else {
            out_ << "  Total Enthalpy (H):  N/A (Temperature " << std::fixed << std::setprecision(2) << T
                 << " K outside Shomate range [" << comp.shomate->T_min << ", " << comp.shomate->T_max << "] K)\n";
        }
    }

    if (root.S_total.has_value()) {
        out_ << "  Total Entropy (S):   " << std::setw(12) << std::fixed << std::setprecision(4) << *root.S_total
             << " J/(mol*K) [Ref: 298.15 K, 1.0 bar ideal gas]\n";
    } else {
        if (comp.formula == "H2O") {
            out_ << "  Total Entropy (S):   N/A (Reference state 298.15 K is below gas-phase Shomate range [500.0, 1700.0] K)\n";
        } else if (!comp.shomate.has_value() || !comp.shomate->is_defined()) {
            out_ << "  Total Entropy (S):   N/A (No ideal-gas Shomate heat capacity defined)\n";
        } else {
            out_ << "  Total Entropy (S):   N/A (Temperature " << std::fixed << std::setprecision(2) << T
                 << " K outside Shomate range [" << comp.shomate->T_min << ", " << comp.shomate->T_max << "] K)\n";
        }
    }

    if (root.stability == MechanicalStability::UNSTABLE) {
        out_ << "  Notice:              Mechanically unstable root within spinodal region; unobservable.\n";
    }
}

void CLIController::print_equilibrium_summary(const PVTResult& pvt) const {
    out_ << "--------------------------------------------------------------------------------\n"
         << "Phase Preference & Equilibrium Analysis:\n";

    if (pvt.roots.size() <= 1) {
        out_ << "  Thermodynamic Status: SINGLE_PHASE (Homogeneous single-phase fluid)\n"
             << "  Phase Coexistence:    None (Exactly 1 physical root exists at this state)\n";
        return;
    }

    const PVTRoot* pref = pvt.preferred_root();
    const PVTRoot* liq = pvt.liquid_root();
    const PVTRoot* vap = pvt.vapor_root();

    if (pvt.equilibrium.preference == PhasePreference::EQUILIBRIUM) {
        out_ << "  Thermodynamic Status: EQUILIBRIUM (Two-Phase Saturated Coexistence)\n"
             << "  Fugacity Matching:    phi_liq (" << std::fixed << std::setprecision(6) << pvt.equilibrium.phi_liq
             << ") == phi_vap (" << pvt.equilibrium.phi_vap << ") within tolerance\n"
             << "  Residual Gibbs Diff:  Delta G_res = " << std::fixed << std::setprecision(3)
             << pvt.equilibrium.delta_G_res << " J/mol (~ 0)\n"
             << "  Conclusion:           Both liquid and vapor phases are simultaneously stable in thermodynamic equilibrium.\n";
    } else if (pvt.equilibrium.preference == PhasePreference::VAPOR_PREFERRED) {
        out_ << "  Thermodynamic Status: VAPOR_PREFERRED (Superheated Vapor)\n"
             << "  Preferred Root:       Vapor-like root (Z = " << std::fixed << std::setprecision(6)
             << (vap ? vap->Z : 0.0) << ")\n"
             << "  Selection Criterion:  phi_vap (" << std::fixed << std::setprecision(6) << pvt.equilibrium.phi_vap
             << ") < phi_liq (" << pvt.equilibrium.phi_liq << ")\n"
             << "  Residual Gibbs Diff:  Delta G_res = G_liq - G_vap = +" << std::fixed << std::setprecision(2)
             << pvt.equilibrium.delta_G_res << " J/mol (> 0)\n"
             << "  Conclusion:           Vapor state has lower chemical potential and represents the stable phase.\n";
    } else if (pvt.equilibrium.preference == PhasePreference::LIQUID_PREFERRED) {
        out_ << "  Thermodynamic Status: LIQUID_PREFERRED (Subcooled Compressed Liquid)\n"
             << "  Preferred Root:       Liquid-like root (Z = " << std::fixed << std::setprecision(6)
             << (liq ? liq->Z : 0.0) << ")\n"
             << "  Selection Criterion:  phi_liq (" << std::fixed << std::setprecision(6) << pvt.equilibrium.phi_liq
             << ") < phi_vap (" << pvt.equilibrium.phi_vap << ")\n"
             << "  Residual Gibbs Diff:  Delta G_res = G_liq - G_vap = " << std::fixed << std::setprecision(2)
             << pvt.equilibrium.delta_G_res << " J/mol (< 0)\n"
             << "  Conclusion:           Liquid state has lower chemical potential and represents the stable phase.\n";
    }

    if (pref) {
        out_ << "  Thermodynamically Selected Root: Z = " << std::fixed << std::setprecision(6) << pref->Z
             << " (" << phase_tag_str(pref->phase) << ")\n";
    }
}

void CLIController::handle_pvt() {
    print_header("1. PENG-ROBINSON PVT BEHAVIOR CALCULATION");

    auto comp_opt = prompt_component();
    if (!comp_opt.has_value()) return;
    const auto& comp = *comp_opt;

    auto T_opt = prompt_temperature();
    if (!T_opt.has_value()) return;
    double T = *T_opt;

    auto P_opt = prompt_pressure();
    if (!P_opt.has_value()) return;
    double P = *P_opt;

    try {
        PengRobinsonEOS eos(comp);
        PVTResult pvt = eos.calculate_pvt(T, P);
        EOSParameters params = eos.parameters(T);

        print_header("PVT CALCULATION RESULTS");
        print_component_summary(comp);
        print_state_conditions(T, P, comp);
        print_eos_parameters(params, pvt.A, pvt.B, pvt.discriminant);

        out_ << "--------------------------------------------------------------------------------\n"
             << "Physical Roots Discovered: " << pvt.roots.size() << " physical root(s) (Z > B = "
             << std::fixed << std::setprecision(6) << pvt.B << ")\n";

        for (size_t i = 0; i < pvt.roots.size(); ++i) {
            print_root_pvt(i, pvt.roots[i], pvt.B);
        }

        print_equilibrium_summary(pvt);
        out_ << std::string(80, '=') << "\n";
    } catch (const std::exception& ex) {
        out_ << "\nCalculation Error: " << ex.what() << "\n";
    }
}

void CLIController::handle_thermo() {
    print_header("2. ENTHALPY (H) AND ENTROPY (S) CALCULATION");

    auto comp_opt = prompt_component();
    if (!comp_opt.has_value()) return;
    const auto& comp = *comp_opt;

    auto T_opt = prompt_temperature();
    if (!T_opt.has_value()) return;
    double T = *T_opt;

    auto P_opt = prompt_pressure();
    if (!P_opt.has_value()) return;
    double P = *P_opt;

    try {
        PengRobinsonEOS eos(comp);
        PVTResult pvt = eos.calculate_pvt(T, P);

        print_header("THERMODYNAMIC PROPERTIES (H & S)");
        print_component_summary(comp);
        print_state_conditions(T, P, comp);

        out_ << "--------------------------------------------------------------------------------\n"
             << "Reference State Convention:\n"
             << "  Ideal gas at T_ref = 298.15 K, P_ref = 100,000 Pa (1.0 bar)\n"
             << "  H_total(T_ref, P_ref) = 0.0 J/mol, S_total(T_ref, P_ref) = 0.0 J/(mol*K)\n"
             << "  Departure Definitions: H^res = H - H^ig(T, P),  S^res = S - S^ig(T, P)\n"
             << "--------------------------------------------------------------------------------\n"
             << "Physical Roots Discovered: " << pvt.roots.size() << " root(s)\n";

        for (size_t i = 0; i < pvt.roots.size(); ++i) {
            print_root_thermo(i, pvt.roots[i], comp, T);
        }

        if (pvt.roots.size() >= 2) {
            const PVTRoot* liq = pvt.liquid_root();
            const PVTRoot* vap = pvt.vapor_root();
            if (liq && vap) {
                out_ << "\nPhase Transition Properties (Vaporization):\n"
                     << "  Delta H_res (Vap - Liq): " << std::fixed << std::setprecision(2)
                     << (vap->H_res - liq->H_res) << " J/mol\n"
                     << "  Delta S_res (Vap - Liq): " << std::fixed << std::setprecision(4)
                     << (vap->S_res - liq->S_res) << " J/(mol*K)\n";
                if (vap->H_total.has_value() && liq->H_total.has_value()) {
                    out_ << "  Enthalpy of Vaporization (Delta H_vap): " << std::fixed << std::setprecision(2)
                         << (*vap->H_total - *liq->H_total) << " J/mol\n";
                }
            }
        }

        print_equilibrium_summary(pvt);
        out_ << std::string(80, '=') << "\n";
    } catch (const std::exception& ex) {
        out_ << "\nCalculation Error: " << ex.what() << "\n";
    }
}

void CLIController::handle_fugacity() {
    print_header("3. FUGACITY & FUGACITY COEFFICIENT CALCULATION");

    auto comp_opt = prompt_component();
    if (!comp_opt.has_value()) return;
    const auto& comp = *comp_opt;

    auto T_opt = prompt_temperature();
    if (!T_opt.has_value()) return;
    double T = *T_opt;

    auto P_opt = prompt_pressure();
    if (!P_opt.has_value()) return;
    double P = *P_opt;

    try {
        PengRobinsonEOS eos(comp);
        PVTResult pvt = eos.calculate_pvt(T, P);

        print_header("PURE-COMPONENT FUGACITY RESULTS");
        print_component_summary(comp);
        print_state_conditions(T, P, comp);

        out_ << "--------------------------------------------------------------------------------\n"
             << "Governing Fugacity Identity:\n"
             << "  ln(phi) = (Z - 1) - ln(Z - B) - [A / (2*sqrt(2)*B)] * ln[(Z + (1+sqrt2)B) / (Z + (1-sqrt2)B)]\n"
             << "  f = phi * P\n"
             << "--------------------------------------------------------------------------------\n"
             << "Physical Roots Discovered: " << pvt.roots.size() << " root(s)\n";

        for (size_t i = 0; i < pvt.roots.size(); ++i) {
            print_root_fugacity(i, pvt.roots[i]);
        }

        print_equilibrium_summary(pvt);
        out_ << std::string(80, '=') << "\n";
    } catch (const std::exception& ex) {
        out_ << "\nCalculation Error: " << ex.what() << "\n";
    }
}

void CLIController::handle_all_properties() {
    print_header("4. COMPREHENSIVE PENG-ROBINSON STATE ANALYSIS");

    auto comp_opt = prompt_component();
    if (!comp_opt.has_value()) return;
    const auto& comp = *comp_opt;

    auto T_opt = prompt_temperature();
    if (!T_opt.has_value()) return;
    double T = *T_opt;

    auto P_opt = prompt_pressure();
    if (!P_opt.has_value()) return;
    double P = *P_opt;

    try {
        PengRobinsonEOS eos(comp);
        PVTResult pvt = eos.calculate_pvt(T, P);
        EOSParameters params = eos.parameters(T);

        print_header("COMPREHENSIVE STATE REPORT");
        print_component_summary(comp);
        print_state_conditions(T, P, comp);
        print_eos_parameters(params, pvt.A, pvt.B, pvt.discriminant);

        out_ << "--------------------------------------------------------------------------------\n"
             << "Physical Roots Discovered: " << pvt.roots.size() << " physical root(s) (Z > B = "
             << std::fixed << std::setprecision(6) << pvt.B << ")\n";

        for (size_t i = 0; i < pvt.roots.size(); ++i) {
            print_root_all(i, pvt.roots[i], comp, T, pvt.B);
        }

        print_equilibrium_summary(pvt);
        out_ << std::string(80, '=') << "\n";
    } catch (const std::exception& ex) {
        out_ << "\nCalculation Error: " << ex.what() << "\n";
    }
}

void CLIController::run_demonstration(int case_num) {
    if (case_num == 1 || case_num == 0) {
        print_header("DEMONSTRATION CASE 1: METHANE AT 300 K, 5.0 MPa (SUPERCRITICAL / GAS)");
        out_ << "Objective: Demonstrate single-phase gas PVT, fugacity, and residual/total H and S.\n"
             << "Expected (PR-1976): Z ~ 0.901844, phi ~ 0.901385, H_res ~ -902.34 J/mol, S_res ~ -2.1446 J/(mol*K)\n";
        const auto& comp = db_.get_component("Methane");
        const double T = 300.0;
        const double P = 5.0e6;
        PengRobinsonEOS eos(comp);
        PVTResult pvt = eos.calculate_pvt(T, P);
        EOSParameters params = eos.parameters(T);

        print_component_summary(comp);
        print_state_conditions(T, P, comp);
        print_eos_parameters(params, pvt.A, pvt.B, pvt.discriminant);
        out_ << "--------------------------------------------------------------------------------\n"
             << "Roots Discovered: " << pvt.roots.size() << " physical root(s)\n";
        for (size_t i = 0; i < pvt.roots.size(); ++i) {
            print_root_all(i, pvt.roots[i], comp, T, pvt.B);
        }
        print_equilibrium_summary(pvt);
        out_ << std::string(80, '=') << "\n\n";
    }

    if (case_num == 2 || case_num == 0) {
        print_header("DEMONSTRATION CASE 2: CARBON DIOXIDE AT 300 K, 7.0 MPa (DENSE COMPRESSED FLUID)");
        out_ << "Objective: Demonstrate dense near-critical fluid with 1 physical liquid-like root.\n"
             << "Expected (PR-1976): Z ~ 0.198057, phi ~ 0.638285, H_res ~ -9374.37 J/mol, S_res ~ -27.5150 J/(mol*K)\n";
        const auto& comp = db_.get_component("CO2");
        const double T = 300.0;
        const double P = 7.0e6;
        PengRobinsonEOS eos(comp);
        PVTResult pvt = eos.calculate_pvt(T, P);
        EOSParameters params = eos.parameters(T);

        print_component_summary(comp);
        print_state_conditions(T, P, comp);
        print_eos_parameters(params, pvt.A, pvt.B, pvt.discriminant);
        out_ << "--------------------------------------------------------------------------------\n"
             << "Roots Discovered: " << pvt.roots.size() << " physical root(s)\n";
        for (size_t i = 0; i < pvt.roots.size(); ++i) {
            print_root_all(i, pvt.roots[i], comp, T, pvt.B);
        }
        print_equilibrium_summary(pvt);
        out_ << std::string(80, '=') << "\n\n";
    }

    if (case_num == 3 || case_num == 0) {
        print_header("DEMONSTRATION CASE 3: PROPANE AT 320 K, 1.0 MPa (THREE PHYSICAL ROOTS & PHASE PREFERENCE)");
        out_ << "Objective: Demonstrate 3 physical roots, unstable intermediate root, and vapor phase preference.\n"
             << "Expected (PR-1976):\n"
             << "  Liquid root:       Z ~ 0.036453, phi ~ 1.239280, H_res ~ -14854.34 J/mol [Stable]\n"
             << "  Intermediate root: Z ~ 0.091826, phi ~ 1.385615, H_res ~ -8210.60 J/mol  [Unstable]\n"
             << "  Vapor root:        Z ~ 0.850555, phi ~ 0.868006, H_res ~ -1131.04 J/mol  [Stable, Preferred]\n"
             << "  Delta G_res:       G_liq - G_vap = +947.42 J/mol (> 0 -> Vapor preferred)\n";
        const auto& comp = db_.get_component("Propane");
        const double T = 320.0;
        const double P = 1.0e6;
        PengRobinsonEOS eos(comp);
        PVTResult pvt = eos.calculate_pvt(T, P);
        EOSParameters params = eos.parameters(T);

        print_component_summary(comp);
        print_state_conditions(T, P, comp);
        print_eos_parameters(params, pvt.A, pvt.B, pvt.discriminant);
        out_ << "--------------------------------------------------------------------------------\n"
             << "Roots Discovered: " << pvt.roots.size() << " physical root(s)\n";
        for (size_t i = 0; i < pvt.roots.size(); ++i) {
            print_root_all(i, pvt.roots[i], comp, T, pvt.B);
        }
        print_equilibrium_summary(pvt);
        out_ << std::string(80, '=') << "\n\n";
    }
}

int CLIController::run() {
    print_header("CHL2002 THERMODYNAMICS SOFTWARE - PENG-ROBINSON EOS SOLVER");
    out_ << "Pure-Component Equation of State & Departure Function Engineering Suite\n"
         << "C++17 Implementation | Thermodynamic Foundations & Rigorous Verification\n";

    while (true) {
        out_ << "\n" << std::string(60, '-') << "\n"
             << "MAIN MENU:\n"
             << "  [1] Calculate PVT properties for a pure component\n"
             << "  [2] Calculate thermodynamic properties H and S (Residual & Total)\n"
             << "  [3] Calculate pure-component fugacity and fugacity coefficient\n"
             << "  [4] Display all available properties for one state (Combined)\n"
             << "  [5] Exit\n"
             << std::string(60, '-') << "\n"
             << "Enter choice [1-5, or 'd' for demonstration cases]: ";

        std::string line;
        if (!read_line(line)) {
            // End of input stream (EOF)
            out_ << "\nEnd of input reached. Exiting gracefully.\n";
            break;
        }

        if (line.empty()) {
            continue;
        }

        if (line == "d" || line == "D" || line == "demo" || line == "demos") {
            out_ << "\nSelect demonstration case:\n"
                 << "  [1] Methane at 300 K, 5 MPa (Supercritical single phase)\n"
                 << "  [2] CO2 at 300 K, 7 MPa (Dense fluid near critical)\n"
                 << "  [3] Propane at 320 K, 1 MPa (Three roots & phase preference)\n"
                 << "  [4] Run all demonstration cases in sequence\n";
            auto demo_choice = prompt_int("Choice", 1, 4, 4);
            if (demo_choice.has_value()) {
                if (*demo_choice == 4) {
                    run_demonstration(0);
                } else {
                    run_demonstration(*demo_choice);
                }
            }
            continue;
        }

        int choice = 0;
        try {
            choice = std::stoi(line);
        } catch (...) {
            out_ << "Error: Unrecognized option '" << line << "'. Please enter 1-5 or 'd'.\n";
            continue;
        }

        switch (choice) {
            case 1:
                handle_pvt();
                break;
            case 2:
                handle_thermo();
                break;
            case 3:
                handle_fugacity();
                break;
            case 4:
                handle_all_properties();
                break;
            case 5:
                out_ << "\nThank you for using CHL2002 Thermodynamics Software. Exiting.\n";
                return 0;
            default:
                out_ << "Error: Choice " << choice << " is not an option. Please select 1 through 5.\n";
                break;
        }
    }

    return 0;
}

} // namespace thermo::cli
