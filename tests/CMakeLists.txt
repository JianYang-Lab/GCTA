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
gcta_add_regression_test(grm.sparse test_grm.py sparse integration)
gcta_add_regression_test(grm.snp_merge test_grm.py snp_merge integration)
gcta_add_regression_test(grm.partitions test_grm.py partitions integration)
gcta_add_regression_test(geno.decoding test_geno.py decoding integration)
gcta_add_regression_test(workflow.greml test_workflows.py greml integration)
gcta_add_regression_test(workflow.pca test_workflows.py pca integration)
gcta_add_regression_test(workflow.fastgwa_lr test_workflows.py fastgwa_lr integration)
gcta_add_regression_test(workflow.fastgwa_lr_nofilter test_workflows.py fastgwa_lr_nofilter integration)
gcta_add_regression_test(workflow.mlma test_workflows.py mlma integration)
gcta_add_regression_test(workflow.fastgwa_exact test_workflows.py fastgwa_exact integration)
gcta_add_regression_test(workflow.fastgwa_mlm test_workflows.py fastgwa_mlm integration)
gcta_add_regression_test(workflow.fastgwa_glmm test_workflows.py fastgwa_glmm integration)
gcta_add_regression_test(workflow.mlma_loco test_workflows.py mlma_loco integration)
gcta_add_regression_test(binary.fastgwa_lr test_binary.py fastgwa_lr_savebin integration)
gcta_add_regression_test(binary.fastgwa_glmm test_binary.py fastgwa_glmm_savebin integration)
gcta_add_regression_test(error.missing_input test_errors.py missing_input integration)
gcta_add_regression_test(error.bad_magic test_errors.py bad_magic integration)
gcta_add_regression_test(error.truncated_bed test_errors.py truncated_bed integration)
gcta_add_regression_test(error.unwritable_output test_errors.py unwritable_output integration)
gcta_add_regression_test(error.grm_short_write test_errors.py grm_short_write integration)
set_tests_properties(error.grm_short_write PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:")
gcta_add_regression_test(error.legacy_invalid_option test_errors.py legacy_invalid_option integration)
gcta_add_regression_test(error.ld_invalid_window test_errors.py ld_invalid_window integration)
gcta_add_regression_test(analysis.ld test_ld.py matrix integration)
gcta_add_regression_test(analysis.acat test_acat.py output integration)
gcta_add_regression_test(error.acat_missing_input test_acat.py missing_input integration)
gcta_add_regression_test(model.lr_roundtrip test_model.py model_lr integration)
gcta_add_regression_test(model.mlm_roundtrip test_model.py model_mlm integration)
gcta_add_regression_test(model.exact_roundtrip test_model.py model_exact integration)
gcta_add_regression_test(model.glmm_roundtrip test_model.py model_glmm integration)
gcta_add_regression_test(format.bgen_subset test_formats.py bgen_subset integration)

add_executable(async_buffer_test async_buffer_test.cpp)
target_include_directories(async_buffer_test PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(async_buffer_test Threads::Threads)
add_test(NAME async.buffer COMMAND async_buffer_test)
set_tests_properties(async.buffer PROPERTIES LABELS "unit" TIMEOUT 120)
gcta_add_regression_test(linkage.numeric_symbols test_linkage.py numeric_symbols integration)
set_tests_properties(linkage.numeric_symbols PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:")
gcta_add_regression_test(backend.grm test_backends.py grm integration)
set_tests_properties(backend.grm PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:")
