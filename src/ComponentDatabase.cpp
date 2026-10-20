#include "thermo/ComponentDatabase.hpp"
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <sstream>

namespace thermo {

std::string ComponentDatabase::normalize(const std::string& str) {
    std::string result;
    result.reserve(str.size());
    for (char ch : str) {
        if (!std::isspace(static_cast<unsigned char>(ch)) && ch != '-' && ch != '_') {
            result.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
        }
    }
    return result;
}

ComponentDatabase::ComponentDatabase() {
    // 1. Methane (CH4)
    add_builtin(PureComponent{
        "Methane",
        "CH4",
        0.016043,     // kg/mol
        190.56,       // K
        4.599e6,      // Pa (45.99 bar)
        0.0115,       // acentric factor
        ShomateParameters{
            -0.703029, 108.4773, -42.52157, 5.862788, 0.678565,
            298.0, 1300.0 // K
        },
        "NIST Chemistry WebBook (SRD 69); Poling et al. (5th Ed.) Table A-1"
    });

    // 2. Ethane (C2H6)
    add_builtin(PureComponent{
        "Ethane",
        "C2H6",
        0.030070,     // kg/mol
        305.32,       // K
        4.872e6,      // Pa (48.72 bar)
        0.0995,       // acentric factor
        ShomateParameters{
            -3.864770, 227.1738, -94.75782, 14.50904, 0.222340,
            298.0, 1200.0 // K
        },
        "NIST Chemistry WebBook (SRD 69); Poling et al. (5th Ed.) Table A-1"
    });

    // 3. Propane (C3H8)
    add_builtin(PureComponent{
        "Propane",
        "C3H8",
        0.044096,     // kg/mol
        369.83,       // K
        4.248e6,      // Pa (42.48 bar)
        0.1523,       // acentric factor
        ShomateParameters{
            -23.16245, 361.6441, -157.9157, 24.63945, 0.985935,
            298.0, 1200.0 // K
        },
        "NIST Chemistry WebBook (SRD 69); Poling et al. (5th Ed.) Table A-1"
    });

    // 4. n-Butane (C4H10)
    add_builtin(PureComponent{
        "n-Butane",
        "C4H10",
        0.058122,     // kg/mol
        425.12,       // K
        3.796e6,      // Pa (37.96 bar)
        0.2002,       // acentric factor
        ShomateParameters{
            -25.29177, 478.1883, -217.4729, 35.10976, 1.050679,
            298.0, 1200.0 // K
        },
        "NIST Chemistry WebBook (SRD 69); Poling et al. (5th Ed.) Table A-1"
    });

    // 5. Nitrogen (N2)
    add_builtin(PureComponent{
        "Nitrogen",
        "N2",
        0.028013,     // kg/mol
        126.20,       // K
        3.398e6,      // Pa (33.98 bar)
        0.0372,       // acentric factor
        ShomateParameters{
            28.98641, 1.853978, -9.647459, 16.63537, 0.000117,
            100.0, 500.0  // K
        },
        "NIST Chemistry WebBook (SRD 69); Poling et al. (5th Ed.) Table A-1"
    });

    // 6. Oxygen (O2)
    add_builtin(PureComponent{
        "Oxygen",
        "O2",
        0.031999,     // kg/mol
        154.58,       // K
        5.043e6,      // Pa (50.43 bar)
        0.0222,       // acentric factor
        ShomateParameters{
            31.32234, -20.23531, 57.86644, -36.50624, -0.007374,
            100.0, 700.0  // K
        },
        "NIST Chemistry WebBook (SRD 69); Poling et al. (5th Ed.) Table A-1"
    });

    // 7. Carbon Dioxide (CO2)
    add_builtin(PureComponent{
        "Carbon Dioxide",
        "CO2",
        0.044010,     // kg/mol
        304.13,       // K
        7.377e6,      // Pa (73.77 bar)
        0.2239,       // acentric factor
        ShomateParameters{
            24.99735, 55.18696, -33.69137, 7.948387, -0.136638,
            298.0, 1200.0 // K
        },
        "NIST Chemistry WebBook (SRD 69); Poling et al. (5th Ed.) Table A-1"
    });

    // 8. Water Vapor (H2O)
    add_builtin(PureComponent{
        "Water",
        "H2O",
        0.018015,     // kg/mol
        647.10,       // K
        22.064e6,     // Pa (220.64 bar)
        0.3443,       // acentric factor
        ShomateParameters{
            30.09200, 6.832514, 6.793435, -2.534480, 0.082139,
            500.0, 1700.0 // K (gas phase)
        },
        "NIST Chemistry WebBook (SRD 69); Poling et al. (5th Ed.) Table A-1"
    });
}

void ComponentDatabase::add_builtin(PureComponent comp) {
    comp.validate();
    const size_t index = components_.size();
    lookup_map_[normalize(comp.name)] = index;
    lookup_map_[normalize(comp.formula)] = index;
    components_.push_back(std::move(comp));
}

bool ComponentDatabase::has_component(const std::string& query) const {
    return lookup_map_.find(normalize(query)) != lookup_map_.end();
}

const PureComponent& ComponentDatabase::get_component(const std::string& query) const {
    const auto it = lookup_map_.find(normalize(query));
    if (it == lookup_map_.end()) {
        std::ostringstream oss;
        oss << "Component '" << query << "' was not found in the component database.";
        throw std::runtime_error(oss.str());
    }
    return components_[it->second];
}

std::vector<std::string> ComponentDatabase::list_component_names() const {
    std::vector<std::string> names;
    names.reserve(components_.size());
    for (const auto& comp : components_) {
        names.push_back(comp.name + " (" + comp.formula + ")");
    }
    return names;
}

void ComponentDatabase::register_custom_component(const PureComponent& component) {
    component.validate();
    const std::string norm_name = normalize(component.name);
    const std::string norm_form = normalize(component.formula);

    if (lookup_map_.find(norm_name) != lookup_map_.end() || lookup_map_.find(norm_form) != lookup_map_.end()) {
        throw std::invalid_argument("A component with name '" + component.name +
                                   "' or formula '" + component.formula + "' already exists.");
    }

    const size_t index = components_.size();
    lookup_map_[norm_name] = index;
    if (!norm_form.empty()) {
        lookup_map_[norm_form] = index;
    }
    components_.push_back(component);
}

size_t ComponentDatabase::size() const noexcept {
    return components_.size();
}

} // namespace thermo
