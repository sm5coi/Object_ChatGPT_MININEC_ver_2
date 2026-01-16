cmake_minimum_required(VERSION 3.16)

project(Object_ChatGPT_MININEC_ver_2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Object_ChatGPT_MININEC_ver_2 main.cpp
    mininec3_compat.hpp
    mininec3_compat.cpp
    vec3.hpp
    mininec3_loader.hpp
    mininec3_loader.cpp
    mininec3_state.hpp
    frequencyInput.hpp
    frequencyInput.cpp)

include(GNUInstallDirs)
install(TARGETS Object_ChatGPT_MININEC_ver_2
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
