cmake_minimum_required(VERSION 3.10)

if (NOT $ENV{CONDA_PREFIX} STREQUAL "")
    message(STATUS "CURRENT CONDA PREFIX: $ENV{CONDA_PREFIX}")
    include_directories($ENV{CONDA_PREFIX}/include)
    link_directories($ENV{CONDA_PREFIX}/lib)
endif ()

set(CMAKE_CXX_COMPILE g++)
set(CMAKE_C_COMPILE gcc)

project(ExtractHAIRs)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXX_FLAGS} -O0 -Wall -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 ")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_subdirectory(lib/bamtools-2.3.0)
#include_directories (lib/bamtools-2.3.0/src)
#include_directories(lib/tclap-1.2.1/include)
aux_source_directory(. EXTRACT_SRC)
aux_source_directory(./htslib HTS_SRC)
aux_source_directory(./samtools SAM_SRC)

#set(LIBSAM libbam.a)
set(LIBHTS libhts.a)

if (NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
    set(CONDA_HTS_LIBRARY "$ENV{CONDA_PREFIX}/lib/libhts.so")
    set(CONDA_ZLIB_LIBRARY "$ENV{CONDA_PREFIX}/lib/libz.so")
endif ()

add_executable(ExtractHAIRs ${EXTRACT_SRC})
if (EXISTS "${CONDA_HTS_LIBRARY}")
    target_link_libraries(ExtractHAIRs "${CONDA_HTS_LIBRARY}")
else ()
    target_link_libraries(ExtractHAIRs -lhts)
endif ()

if (EXISTS "${CONDA_ZLIB_LIBRARY}")
    target_link_libraries(ExtractHAIRs "${CONDA_ZLIB_LIBRARY}")
else ()
    target_link_libraries(ExtractHAIRs z)
endif ()
