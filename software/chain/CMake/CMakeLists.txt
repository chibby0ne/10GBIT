project(LDPC_QUANT)
cmake_minimum_required(VERSION 2.6)

# Set Build type to release if none is given
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Compiler options for MSVC
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/O2 /arch:SSE2")
  set(CMAKE_CXX_FLAGS_RELEASE "/MT /O2 /Ob2 /D NDEBUG")
  add_definitions(XML_LIBRARY)
ENDIF(MSVC)

# Compiler and linker options for gcc
IF(CMAKE_COMPILER_IS_GNUCXX)
  # C code compiler
  set(CMAKE_C_FLAGS "-Wall")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -march=nocona -mmmx -msse -msse2 -msse3 -mfpmath=sse -fomit-frame-pointer -fstrict-aliasing -momit-leaf-frame-pointer -fno-tree-pre -falign-loops -ftree-vectorize")
  set(CMAKE_C_FLAGS_DEBUG "-ggdb")
  
  # C++ code compiler
  set(CMAKE_CXX_FLAGS "-Wall")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=nocona -mmmx -msse -msse2 -msse3 -mfpmath=sse -fomit-frame-pointer -fstrict-aliasing -momit-leaf-frame-pointer -fno-tree-pre -falign-loops -ftree-vectorize")
  set(CMAKE_CXX_FLAGS_DEBUG "-ggdb")
  set(CMAKE_CXX_FLAGS_PROF "-pg -O3 -march=nocona -mmmx -msse -msse2 -msse3 -mfpmath=sse -fstrict-aliasing -momit-leaf-frame-pointer -fno-tree-pre -falign-loops -ftree-vectorize")
  
  # Linker
  set(CMAKE_EXE_LINKER_FLAGS "")
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-static-libgcc")
  
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# cse library directories
set(CSE_LIBRARY_HOME_PATH "${LDPC_QUANT_SOURCE_DIR}/../../cse")
include_directories(${CSE_LIBRARY_HOME_PATH}/include)
link_directories (${CSE_LIBRARY_HOME_PATH}/lib)

#ems library directories
set(EMS_LIBRARY_HOME_PATH "${LDPC_QUANT_SOURCE_DIR}/../../ems")
include_directories(${EMS_LIBRARY_HOME_PATH}/include)
link_directories (${EMS_LIBRARY_HOME_PATH}/lib)

# itpp library directories
set(ITPP_LIBRARY_HOME_PATH "/usr/local")
include_directories(${ITPP_LIBRARY_HOME_PATH}/itpp/include)
link_directories (${ITPP_LIBRARY_HOME_PATH}/lib)

# additional source directories
set(SRC_TOP_DIR ${LDPC_QUANT_SOURCE_DIR}/..)
include_directories(${SRC_TOP_DIR}/ldpc_dec)
AUX_SOURCE_DIRECTORY(${SRC_TOP_DIR}/ldpc_dec LDPC_DEC_SOURCES)
include_directories(${SRC_TOP_DIR}/ldpc_enc)
AUX_SOURCE_DIRECTORY(${SRC_TOP_DIR}/ldpc_enc  LDPC_ENC_SOURCES)

# Make executable and target LDPC_QUANT
add_executable (LDPC_QUANT ${LDPC_QUANT_SOURCE_DIR}/../WPAN_chain ${LDPC_ENC_SOURCES} ${LDPC_DEC_SOURCES})

# Set directory for executable
set_target_properties(LDPC_QUANT PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${LDPC_QUANT_SOURCE_DIR}/../bin)

# Set the name of the executable
set_target_properties(LDPC_QUANT PROPERTIES RELEASE_OUTPUT_NAME "ieee802_3an_chain_release")
set_target_properties(LDPC_QUANT PROPERTIES DEBUG_OUTPUT_NAME   "ieee802_3an_chain_debug")
set_target_properties(LDPC_QUANT PROPERTIES PROF_OUTPUT_NAME    "ieee802_3an_chain_prof")


# target_link_libraries (LDPC_QUANT cse_static ems_static itpp)
target_link_libraries (LDPC_QUANT cse ems itpp)

execute_process(COMMAND ctags -R WORKING_DIRECTORY ${LDPC_QUANT_SOURCE_DIR}/../../../.)

