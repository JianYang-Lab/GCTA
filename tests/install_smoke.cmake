if(NOT DEFINED GCTA_INSTALL_TEST_BINARY_DIR)
    message(FATAL_ERROR "GCTA_INSTALL_TEST_BINARY_DIR is required")
endif()
if(NOT DEFINED GCTA_INSTALL_TEST_PREFIX)
    message(FATAL_ERROR "GCTA_INSTALL_TEST_PREFIX is required")
endif()
if(NOT DEFINED GCTA_INSTALL_TEST_WORK)
    message(FATAL_ERROR "GCTA_INSTALL_TEST_WORK is required")
endif()
if(NOT DEFINED GCTA_INSTALL_TEST_DATA)
    message(FATAL_ERROR "GCTA_INSTALL_TEST_DATA is required")
endif()

set(ENV{DESTDIR} "${GCTA_INSTALL_TEST_WORK}/root")
execute_process(
    COMMAND "${CMAKE_COMMAND}" --install "${GCTA_INSTALL_TEST_BINARY_DIR}" --config Release
    RESULT_VARIABLE install_result
    OUTPUT_VARIABLE install_output
    ERROR_VARIABLE install_error
)
unset(ENV{DESTDIR})
if(NOT install_result EQUAL 0)
    message(FATAL_ERROR "Installation failed:\n${install_output}\n${install_error}")
endif()

set(installed_gcta "${GCTA_INSTALL_TEST_WORK}/root${GCTA_INSTALL_TEST_PREFIX}/bin/gcta64")
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env --unset=LD_LIBRARY_PATH
        "${installed_gcta}"
        --bfile "${GCTA_INSTALL_TEST_DATA}/gcta_legacy/toy"
        --make-grm
        --out "${GCTA_INSTALL_TEST_WORK}/smoke"
        --thread-num 1
    RESULT_VARIABLE smoke_result
    OUTPUT_VARIABLE smoke_output
    ERROR_VARIABLE smoke_error
)
if(NOT smoke_result EQUAL 0)
    message(FATAL_ERROR "Installed GCTA failed:\n${smoke_output}\n${smoke_error}")
endif()
foreach(suffix IN ITEMS grm.bin grm.N.bin grm.id)
    if(NOT EXISTS "${GCTA_INSTALL_TEST_WORK}/smoke.${suffix}")
        message(FATAL_ERROR "Missing installed-run output: smoke.${suffix}")
    endif()
endforeach()
message(STATUS "Installed GCTA ran without LD_LIBRARY_PATH")
