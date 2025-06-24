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

FetchContent_MakeAvailable(Boost Eigen spectra)

find_package(Threads REQUIRED)

# NOTE: Install gsl sqlite3 zstd with your system package manager
find_package(GSL REQUIRED)
find_package(SQLite3 REQUIRED)
find_package(zstd REQUIRED CONFIG)

if(LINUX)
    if(NOT OMP_LIBRARY)
        set(OMP_LIBRARY "/usr/lib/libgomp.so")
    endif()
    if(NOT MKL_DIR)
        set(MKL_DIR "/opt/intel/mkl/lib/intel64/cmake/mkl")
    endif()
    set(MKL_THREADING "gnu_thread")
endif()
find_package(OpenMP REQUIRED)
find_package(MKL REQUIRED)

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

message(STATUS "ARCHITECTURE:${CMAKE_LIBRARY_ARCHITECTURE}")
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
        "${CMAKE_SOURCE_DIR}/submods/plink-ng/2.0"
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

add_subdirectory(submods)

if(APPLE)
    target_link_directories(gcta64 ${MKLROOT}/lib)
    target_link_directories(gcta64 /usr/local/lib)
elseif(LINUX AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_directories(gcta64 ${MKLROOT}/../../../lib/intel64)
endif()

target_include_directories(gcta64 PUBLIC
    ${COMMON_INCLUDES}
    ${MKL_INCLUDE}
    # PLINK2.0 pgen library;
    "${CMAKE_SOURCE_DIR}/submods/plink-ng/2.0"
    "${CMAKE_SOURCE_DIR}/submods/plink-ng/2.0/simde")

if(APPLE)
    # iomp5
    target_link_libraries(gcta64 ${GCTA_LIBS} Pgenlib mainV1 z ${BLAS_LIB} sqlite3 zstd omp pthread m dl gsl gslcblas)
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
        target_link_libraries(gcta64 mainV1 ${GCTA_LIBS} Pgenlib z m GSL::gsl GSL::gslcblas sqlite3 zstd MKL::MKL Threads::Threads)
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

