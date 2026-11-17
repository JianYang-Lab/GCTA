set(GCTA_TEST_DATA "${CMAKE_CURRENT_SOURCE_DIR}/data")
set(GCTA_TEST_REGRESSION "${CMAKE_CURRENT_SOURCE_DIR}/regression")

add_test(
    NAME data.integrity
    COMMAND ${Python3_EXECUTABLE} "${GCTA_TEST_DATA}/download.py" --check
)
set_tests_properties(data.integrity PROPERTIES LABELS "data" TIMEOUT 30)

function(gcta_add_grm_test scenario label)
	set(test_name "grm.${scenario}")
	add_test(
		NAME ${test_name}
		COMMAND ${Python3_EXECUTABLE} "${GCTA_TEST_REGRESSION}/test_grm.py"
			--gcta $<TARGET_FILE:gcta64>
			--data "${GCTA_TEST_DATA}"
			--work "${CMAKE_CURRENT_BINARY_DIR}/output/${test_name}"
			--scenario ${scenario}
	)
	set_tests_properties(${test_name} PROPERTIES
		LABELS "${label}"
		TIMEOUT 120
		ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1"
	)
endfunction()

gcta_add_grm_test(standard unit)
gcta_add_grm_test(weighted unit)
gcta_add_grm_test(threads integration)
gcta_add_grm_test(partitions integration)
