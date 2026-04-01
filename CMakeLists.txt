# ------------------------------------------------------------
# AMORE2Onl top-level CMakeLists.txt (Linux Optimized & Relocatable)
# ------------------------------------------------------------
cmake_minimum_required ( VERSION 3.16 )
project ( AMORE2Onl VERSION 1.0.0 LANGUAGES CXX )

# ---- C++ Standard ----
set ( CMAKE_CXX_STANDARD 17 )
set ( CMAKE_CXX_STANDARD_REQUIRED ON )
set ( CMAKE_CXX_EXTENSIONS OFF )

# Use standard GNU installation directories
include ( GNUInstallDirs )

# ---- Relocatable RPATH Settings (Linux) ----
set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
set ( CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}" )

# ---- ROOT ----
find_package ( ROOT REQUIRED COMPONENTS Core RIO )
include ( ${ROOT_USE_FILE} )

# Prevent absolute paths from being hardcoded in ROOT dictionaries
set ( CMAKE_ROOT_DICT_ARGS --noinclude )
set ( ROOT_GENERATE_DICTIONARY_RELATIVE ON )

# ---- libusb (via pkg-config) ----
find_package ( PkgConfig REQUIRED )
pkg_check_modules ( LIBUSB REQUIRED IMPORTED_TARGET libusb-1.0 )

# ---- yaml-cpp for configuration ----
find_package ( yaml-cpp REQUIRED )

# ---- zeromq (cppzmq) for message server ----
find_package ( cppzmq REQUIRED )

# ---- HDF5 (Required) ----
# Find HDF5 and halt the build process with a fatal error if not found.
find_package ( HDF5 REQUIRED COMPONENTS C CXX HL )
message ( STATUS "HDF5 found (${HDF5_VERSION}): enabling HDF5-based utilities" )
add_compile_definitions ( ENABLE_HDF5 )

# ---- Submodule: CUPDAQ ----
# Force the CUPDAQ submodule to enable HDF5 features since it is now required.
set ( CUPDAQ_ENABLE_HDF5 ON CACHE BOOL "Force HDF5 in CUPDAQ" FORCE )

# Integrate the CUPDAQ core system
message ( STATUS "Adding CUPDAQ submodule into AMORE2Onl..." )
add_subdirectory ( CUPDAQ )

# ---- Global Compiler Settings ----
add_compile_options ( -O2 )

# ---- AMORE2Onl Specific Subdirectories ----
# Add your specific modules here when ready
add_subdirectory ( AMORE )
add_subdirectory ( AMOREAlgs )
add_subdirectory ( AMORESystem )
add_subdirectory ( AMOREHDF5 )
add_subdirectory ( AMOREDAQ )

# ---- Environment Setup Script Generation (setup_amore2onl.sh) ----
# Create a unified setup script for both AMORE2Onl and CUPDAQ
file ( WRITE ${CMAKE_BINARY_DIR}/setup_amore2onl.sh
  "#!/bin/bash
# AMORE2Onl & CUPDAQ Unified Environment Setup Script
# Generated automatically by CMake.

export AMORE2ONL_DIR=\"\$(builtin cd \"\$(dirname \"\${BASH_SOURCE[0]}\")\" && pwd)\"

# Update system PATH
export PATH=\"\$AMORE2ONL_DIR/bin:\$PATH\"

# Update library search path
export LD_LIBRARY_PATH=\"\$AMORE2ONL_DIR/${CMAKE_INSTALL_LIBDIR}:\$LD_LIBRARY_PATH\"

# Update ROOT include path
export ROOT_INCLUDE_PATH=\"\$AMORE2ONL_DIR/include:\$ROOT_INCLUDE_PATH\"

echo \"[AMORE2Onl] Environment loaded from: \$AMORE2ONL_DIR\"
" )

install ( FILES ${CMAKE_BINARY_DIR}/setup_amore2onl.sh DESTINATION . )
