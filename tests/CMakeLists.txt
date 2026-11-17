set(GCTA_TEST_DATA "${CMAKE_CURRENT_SOURCE_DIR}/data")
set(GCTA_TEST_REGRESSION "${CMAKE_CURRENT_SOURCE_DIR}/regression")

add_test(
    NAME data.integrity
    COMMAND ${Python3_EXECUTABLE} "${GCTA_TEST_DATA}/download.py" --check
)
set_tests_properties(data.integrity PROPERTIES LABELS "data" TIMEOUT 30)

function(gcta_add_regression_test test_name script scenario label)
	add_test(
		NAME ${test_name}
		COMMAND ${Python3_EXECUTABLE} "${GCTA_TEST_REGRESSION}/${script}"
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

gcta_add_regression_test(grm.standard test_grm.py standard unit)
gcta_add_regression_test(grm.weighted test_grm.py weighted unit)
gcta_add_regression_test(grm.missing test_grm.py missing unit)
gcta_add_regression_test(grm.threads test_grm.py threads integration)
gcta_add_regression_test(grm.partitions test_grm.py partitions integration)
gcta_add_regression_test(geno.decoding test_geno.py decoding integration)
