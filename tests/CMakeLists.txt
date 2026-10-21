# Tests for cubic polynomial solver
add_executable(test_cubic_solver test_cubic_solver.cpp)
target_link_libraries(test_cubic_solver PRIVATE thermo_core)
add_test(NAME test_cubic_solver COMMAND test_cubic_solver)

# Tests for components and Shomate heat capacity
add_executable(test_components test_components.cpp)
target_link_libraries(test_components PRIVATE thermo_core)
add_test(NAME test_components COMMAND test_components)

# Tests for PR EOS and PVT calculations
add_executable(test_pvt_eos test_pvt_eos.cpp)
target_link_libraries(test_pvt_eos PRIVATE thermo_core)
add_test(NAME test_pvt_eos COMMAND test_pvt_eos)

# Tests for residual properties, fugacity, total H/S, and phase equilibrium (Session 4)
add_executable(test_thermo_properties test_thermo_properties.cpp)
target_link_libraries(test_thermo_properties PRIVATE thermo_core)
add_test(NAME test_thermo_properties COMMAND test_thermo_properties)

