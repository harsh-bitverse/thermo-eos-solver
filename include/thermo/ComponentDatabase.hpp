#ifndef THERMO_COMPONENT_DATABASE_HPP
#define THERMO_COMPONENT_DATABASE_HPP

#include "thermo/Component.hpp"
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>

namespace thermo {

/**
 * @brief Catalog of pure fluid thermodynamic properties.
 * Provides curated datasets for 8 industrial gases and allows user-registered custom fluids.
 */
class ComponentDatabase {
public:
    ComponentDatabase();

    /**
     * @brief Checks if a component exists by name or chemical formula (case-insensitive).
     */
    [[nodiscard]] bool has_component(const std::string& query) const;

    /**
     * @brief Retrieves a component by name or chemical formula (case-insensitive).
     * @throws std::runtime_error if component is not found.
     */
    [[nodiscard]] const PureComponent& get_component(const std::string& query) const;

    /**
     * @brief Returns a list of all component display names in the database.
     */
    [[nodiscard]] std::vector<std::string> list_component_names() const;

    /**
     * @brief Registers a custom user-defined component into the database.
     * Validates physical properties before registration.
     * @throws std::invalid_argument if the component data is invalid.
     */
    void register_custom_component(const PureComponent& component);

    /**
     * @brief Returns total count of registered components.
     */
    [[nodiscard]] size_t size() const noexcept;

private:
    std::vector<PureComponent> components_;
    std::unordered_map<std::string, size_t> lookup_map_; // Normalized string -> index in components_

    void add_builtin(PureComponent comp);
    static std::string normalize(const std::string& str);
};

} // namespace thermo

#endif // THERMO_COMPONENT_DATABASE_HPP
