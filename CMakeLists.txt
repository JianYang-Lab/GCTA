CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

project(GCTA2 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
    message("Setting to release build")
endif()

# Fetch 3rd parties
include(FetchContent)

set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
set(FETCHCONTENT_QUIET OFF)

set(BOOST_ENABLE_CMAKE ON)
set(BOOST_INCLUDE_LIBRARIES algorithm math crc)
FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.84.0/boost-1.84.0.tar.xz
    URL_MD5 893b5203b862eb9bbd08553e24ff146a
    DOWNLOAD_EXTRACT_TIMESTAMP ON)

FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE)

FetchContent_Declare(
    spectra
    GIT_REPOSITORY https://github.com/yixuan/spectra.git
    GIT_TAG v1.1.0
    GIT_SHALLOW TRUE)

set(ZSTD_BUILD_SHARED ON)
FetchContent_Declare(
    zstd
    GIT_REPOSITORY https://github.com/facebook/zstd.git
    GIT_TAG v1.5.6
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR build/cmake
)

FetchContent_MakeAvailable(Boost Eigen spectra zstd)

find_package(Threads REQUIRED)

# NOTE: Install gsl sqlite3 with your system package manager
find_package(GSL REQUIRED)
find_package(SQLite3 REQUIRED)

# Find OpenMP
if(APPLE)
  # Use manual flags for OpenMP with Apple Clang
  set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_C_INCLUDE_DIR "/opt/homebrew/opt/libomp/include")
  set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
  set(OpenMP_CXX_LIB_NAMES "omp")
  set(OpenMP_CXX_INCLUDE_DIR "/opt/homebrew/opt/libomp/include")
  set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif()

if(LINUX)
    if(NOT OMP_LIBRARY)
        set(OMP_LIBRARY "/usr/lib/libgomp.so")
    endif()
    if(NOT MKL_DIR)
        set(MKL_DIR "/opt/intel/mkl/lib/intel64/cmake/mkl")
    endif()
    set(MKL_INTERFACE "lp64")
    set(MKL_THREADING "gnu_thread")
    find_package(MKL REQUIRED)
endif()

find_package(OpenMP REQUIRED)

# get system architecture
EXECUTE_PROCESS(COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE)
set(CMAKE_LIBRARY_ARCHITECTURE ${ARCHITECTURE})

# fix win stack size issue
if(WIN32)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--stack,10000000")
    MESSAGE(FATAL_ERROR "Use the Visual Studio project file which located in gcta_win64 directory")
endif()

if("${CMAKE_SYSTEM_NAME}" MATCHES "Linux" OR "${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    MESSAGE("Linux and Mac are supported.")
else()
    MESSAGE("Your system: ${CMAKE_SYSTEM_NAME}")
    MESSAGE(FATAL_ERROR "Only Linux and Mac are currently supported.")
endif()

if (CMAKE_COMPILER_IS_GNUCC)
    MESSAGE("Compiling by GCC${CMAKE_CXX_COMPILER_VERSION}")
    if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0))
        MESSAGE(FATAL_ERROR "GCTA can't build without GCC version >= 6.1")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    MESSAGE("Compiling by Clang${CMAKE_CXX_COMPILER_VERSION}")
    if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.0))
        MESSAGE(FATAL_ERROR "GCTA can't build without Clang version >= 8.0")
    endif()
else()
    MESSAGE(FATAL_ERROR "GCTA can only compile with GCC or Clang")
endif()

message(STATUS "arch:${CMAKE_LIBRARY_ARCHITECTURE}, compiler_id:${CMAKE_CXX_COMPILER_ID}")
if(CMAKE_LIBRARY_ARCHITECTURE STREQUAL "aarch64")
    SET(OPENBLAS "$ENV{OPENBLAS}")
    IF(NOT OPENBLAS)
        MESSAGE(FATAL_ERROR "Specify environment variable OPENBLAS to the location of OPENBLAS package")
    ENDIF()
    INCLUDE_DIRECTORIES("${OPENBLAS}/include")

    find_library(lib_openblas NAMES openblas PATHS "${OPENBLAS}" "${OPENBLAS}/lib")
    set(BLAS_LIB ${lib_openblas})

    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -pthread")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
        add_compile_options(-openmp)
        set(CMAKE_EXE_LINKER_FLAGS "-Wl,-no_implicit_dylibs")
        if(NOT APPLE)
            # don't delete this it will result in _implibs
            set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        endif()
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -pthread -fvisibility=hidden -fvisibility-inlines-hidden")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG") #-flto
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fno-inline -fno-implicit-inline-templates -g3")
        add_compile_options(-fopenmp)
        set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -pthread")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ipo -DNDEBUG")
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
    endif()

else()
    IF(NOT SSE)
        set(SSE "-msse2")
        message("SSE build")
    ENDIF()

    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 ${SSE} -DMKL_LP64 -pthread")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
        add_compile_options(-openmp)
        set(CMAKE_EXE_LINKER_FLAGS "-Wl,-no_implicit_dylibs")
        if(NOT APPLE)
            # don't delete this it will result in _implibs
            set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        endif()
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 ${SSE} -pthread -DMKL_LP64 -fvisibility=hidden -fvisibility-inlines-hidden")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG") #-flto
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fno-inline -fno-implicit-inline-templates -g3")
        add_compile_options(-fopenmp)
        # set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -axAVX2 -pthread")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ipo -DNDEBUG")
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        message(FALTAL_ERROR "GCTA can not compile by Microsft C++ compiler in CMake currently, please use gcta_win64 folder to build")
    endif()
endif()

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
set(MAIN_SOURCE "${PROJECT_SOURCE_DIR}/src/main.cpp")

file(GLOB SRCS "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRCS "${MAIN_SOURCE}")
list(APPEND SRCS "${PROJECT_SOURCE_DIR}/src/acat.c")

set(GCTA_LIBS "")
foreach(lib_src ${SRCS})
    get_filename_component(file_name "${lib_src}" NAME_WE)
    string(TOLOWER ${file_name} lib_name)
    add_library(${lib_name} "${lib_src}")

    target_include_directories(${lib_name} PUBLIC
        ${COMMON_INCLUDES}
        ${MKL_INCLUDE}
        # PLINK2.0 pgen library;
        "${CMAKE_SOURCE_DIR}/submods/plink-ng/2.0/include"
        "${CMAKE_SOURCE_DIR}/submods/plink-ng/2.0/simde")

    target_link_libraries(${lib_name}
        Boost::algorithm Boost::math Boost::crc
        Eigen3::Eigen
        Spectra)

    list(APPEND GCTA_LIBS ${lib_name})
endforeach(lib_src)

add_executable(gcta64 ${MAIN_SOURCE})

add_subdirectory(main)
target_include_directories(mainV1 PUBLIC ${COMMON_INCLUDES} ${MKL_INCLUDE})
target_link_libraries(mainV1 Eigen3::Eigen)

set(GCTA_TARGETS gcta64 mainV1 ${GCTA_LIBS})
function(targets_includes targets include_dir)
    foreach(target IN LISTS ${targets})
        target_include_directories(${target} PUBLIC ${include_dir})
    endforeach()
endfunction()

if(APPLE)
    if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "arm64")
        # For Eigen backend, OpenBLAS, LAPACK
        # Get macOS SDK path
        execute_process(
            COMMAND xcrun --sdk macosx --show-sdk-path
            OUTPUT_VARIABLE MACOSX_SDK_PATH
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(VECLIB_INCLUDE_DIR "${MACOSX_SDK_PATH}/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers")
        set(OPENBLAS_INCLUDE_DIR "/opt/homebrew/opt/openblas/include")
        set(GSL_INCLUDE_DIR "/opt/homebrew/opt/gsl/include")

        targets_includes(GCTA_TARGETS ${VECLIB_INCLUDE_DIR})
        targets_includes(GCTA_TARGETS ${OPENBLAS_INCLUDE_DIR})
        targets_includes(GCTA_TARGETS ${GSL_INCLUDE_DIR})
        targets_includes(GCTA_TARGETS ${zstd_SOURCE_DIR}/lib)

        target_link_libraries(gcta64 OpenMP::OpenMP_CXX)
        target_link_libraries(mainV1 OpenMP::OpenMP_CXX)
        foreach(gcta_lib ${GCTA_LIBS})
            target_link_libraries(${gcta_lib} OpenMP::OpenMP_CXX)
        endforeach()

        find_library(ACCELERATE_LIBRARY Accelerate)
        target_link_libraries(gcta64 ${ACCELERATE_LIBRARY})
        target_link_libraries(mainV1 ${ACCELERATE_LIBRARY})
        foreach(gcta_lib ${GCTA_LIBS})
            target_link_libraries(${gcta_lib} ${ACCELERATE_LIBRARY})
        endforeach()

    elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        # Do nothing
    endif()
endif()

add_subdirectory(submods)

if(LINUX AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_directories(gcta64 ${MKLROOT}/../../../lib/intel64)
endif()

target_include_directories(gcta64 PUBLIC
    ${COMMON_INCLUDES}
    ${MKL_INCLUDE}
    # PLINK2.0 pgen library;
    "${CMAKE_SOURCE_DIR}/submods/plink-ng/2.0/include"
    "${CMAKE_SOURCE_DIR}/submods/plink-ng/2.0/simde")

if(APPLE)
    # iomp5
    target_link_libraries(gcta64 mainV1 ${GCTA_LIBS} Pgenlib z m GSL::gsl GSL::gslcblas SQLite::SQLite3 libzstd_shared Threads::Threads)
else()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # error keeps only
        #target_link_libraries(gcta64 ${GCTA_LIBS} mainV1 z ${MKL_LIB} iomp5 pthread m dl )
        #target_link_libraries(gcta64 mainV1 ${GCTA_LIBS} z ${MKL_LIB} pthread m dl)
        # Dynamic link works
        #target_link_libraries(gcta64 mainV1 ${GCTA_LIBS} -Wl,-Bstatic z -Wl,--start-group ${MKL_LIB} -Wl,--end-group iomp5 -Wl,-Bdynamic pthread m dl)
        target_link_libraries(gcta64 mainV1 ${GCTA_LIBS} Pgenlib z sqlite3 zstd -Wl,--start-group ${BLAS_LIB} -Wl,--end-group iomp5 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive m dl)
    else()
        #target_link_libraries(gcta64 mainV1 ${GCTA_LIBS} -Wl,-Bstatic z -Wl,--start-group ${MKL_LIB} -Wl,--end-group gomp -Wl,-Bdynamic pthread m dl)
        #target_link_libraries(gcta64 mainV1 ${GCTA_LIBS} Pgenlib -static z sqlite3 zstd -Wl,--start-group ${MKL_LIB} -Wl,--end-group gomp -Wl,--whole-archive -lpthread -Wl,--no-whole-archive m dl)
        target_link_libraries(gcta64 mainV1 ${GCTA_LIBS} Pgenlib z m GSL::gsl GSL::gslcblas sqlite3 libzstd_shared MKL::MKL Threads::Threads)
    endif()
endif()

# Testing has some problems currently
#enable_testing()
#ADD_SUBDIRECTORY(test)

IF(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
    )
ENDIF()

# List all targets
get_property(imported_targets DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY IMPORTED_TARGETS)
foreach(target IN LISTS imported_targets)
    message(STATUS "imported target: ${target}")
endforeach()

get_property(local_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
foreach(t IN LISTS local_targets)
    message(STATUS "built target: ${t}")
endforeach()

# Installation
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)

# Targets
install(TARGETS gcta64 mainV1 ${GCTA_LIBS}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

# Resource
install(FILES ${CMAKE_SOURCE_DIR}/resources/gcta64.desktop
    DESTINATION share/applications)
install(FILES ${CMAKE_SOURCE_DIR}/resources/gcta64.png
    DESTINATION share/icons/hicolor/256x256/apps)
install(FILES ${CMAKE_SOURCE_DIR}/resources/gcta64.svg
    DESTINATION share/icons/hicolor/scalable/apps)
