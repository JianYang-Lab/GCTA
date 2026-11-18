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
gcta_add_regression_test(grm.xchr test_grm_x.py x integration)
gcta_add_regression_test(grm.dominance test_grm_x.py dominance integration)
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
gcta_add_regression_test(workflow.bivar_reml test_bivar.py reml integration)
gcta_add_regression_test(workflow.hereg test_he.py output integration)
gcta_add_regression_test(error.bivar_priors test_bivar.py bad_priors integration)
gcta_add_regression_test(workflow.cojo_joint test_summary.py cojo_joint integration)
gcta_add_regression_test(workflow.fastbat test_summary.py fastbat integration)
gcta_add_regression_test(workflow.mbat test_summary.py mbat integration)
gcta_add_regression_test(workflow.gsmr test_gsmr.py output integration)
gcta_add_regression_test(workflow.gsmr_compressed_flip test_gsmr.py compressed_flip integration)
gcta_add_regression_test(workflow.gsmr_multi_trait test_gsmr.py multi_trait integration)
gcta_add_regression_test(workflow.gsmr_correlated test_gsmr.py correlated integration)
gcta_add_regression_test(workflow.gsmr_heidi test_gsmr.py heidi integration)
gcta_add_regression_test(workflow.gsmr_overlap_ldsc test_gsmr.py overlap_ldsc integration)
gcta_add_regression_test(workflow.gsmr_overlap_correlation test_gsmr.py overlap_correlation integration)
gcta_add_regression_test(error.gsmr_overlap_invalid test_gsmr.py invalid_overlap integration)
gcta_add_regression_test(workflow.mtcojo test_mtcojo.py output integration)
gcta_add_regression_test(workflow.mtcojo_compressed_flip test_mtcojo.py compressed_flip integration)
gcta_add_regression_test(workflow.mtcojo_internal_bxy test_mtcojo.py internal_bxy integration)
gcta_add_regression_test(error.mtcojo_missing_ld test_mtcojo.py missing_ld integration)
gcta_add_regression_test(workflow.erm_std test_expression.py erm_std integration)
gcta_add_regression_test(workflow.erm_center test_expression.py erm_center integration)
gcta_add_regression_test(workflow.erm_individual test_expression.py erm_individual integration)
gcta_add_regression_test(workflow.erm_gz test_expression.py erm_gz integration)
gcta_add_regression_test(error.expression_duplicate_probe test_expression.py duplicate_probe integration)
gcta_add_regression_test(workflow.ecojo_select test_expression.py ecojo_select integration)
gcta_add_regression_test(workflow.ecojo_blup test_expression.py ecojo_blup integration)
gcta_add_regression_test(workflow.ecojo_no_signal test_expression.py ecojo_no_signal integration)
gcta_add_regression_test(workflow.ecojo_collinear test_expression.py ecojo_collinear integration)
gcta_add_regression_test(format.raw_modern test_raw.py modern integration)
gcta_add_regression_test(format.raw_old test_raw.py old integration)
gcta_add_regression_test(error.raw_snp_mismatch test_raw.py mismatch integration)
gcta_add_regression_test(workflow.simu_qt test_simu.py qt integration)
gcta_add_regression_test(workflow.simu_cc test_simu.py cc integration)
gcta_add_regression_test(error.simu_seed test_simu.py bad_seed integration)
gcta_add_regression_test(workflow.fst test_fst.py output integration)
gcta_add_regression_test(binary.fastgwa_lr test_binary.py fastgwa_lr_savebin integration)
gcta_add_regression_test(binary.fastgwa_glmm test_binary.py fastgwa_glmm_savebin integration)
gcta_add_regression_test(error.missing_input test_errors.py missing_input integration)
gcta_add_regression_test(error.bad_magic test_errors.py bad_magic integration)
gcta_add_regression_test(error.truncated_bed test_errors.py truncated_bed integration)
gcta_add_regression_test(error.unwritable_output test_errors.py unwritable_output integration)
gcta_add_regression_test(error.grm_short_write test_errors.py grm_short_write integration)
set_tests_properties(error.grm_short_write PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:")
gcta_add_regression_test(error.legacy_invalid_option test_errors.py legacy_invalid_option integration)
gcta_add_regression_test(legacy.grm_v1 test_legacy.py grm integration)
gcta_add_regression_test(legacy.freq_v1 test_legacy.py frequency integration)
gcta_add_regression_test(cli.threads_alias test_legacy.py threads_alias integration)
gcta_add_regression_test(error.thread_conflict test_legacy.py thread_conflict integration)
gcta_add_regression_test(error.ld_invalid_window test_errors.py ld_invalid_window integration)
gcta_add_regression_test(analysis.ld test_ld.py matrix integration)
gcta_add_regression_test(legacy.ld_score test_ld_legacy.py score integration)
gcta_add_regression_test(legacy.ld_pruning test_ld_legacy.py pruning integration)
gcta_add_regression_test(analysis.acat test_acat.py output integration)
gcta_add_regression_test(error.acat_missing_input test_acat.py missing_input integration)
gcta_add_regression_test(model.lr_roundtrip test_model.py model_lr integration)
gcta_add_regression_test(model.mlm_roundtrip test_model.py model_mlm integration)
gcta_add_regression_test(model.exact_roundtrip test_model.py model_exact integration)
gcta_add_regression_test(model.glmm_roundtrip test_model.py model_glmm integration)
gcta_add_regression_test(format.bgen_subset test_formats.py bgen_subset integration)
gcta_add_regression_test(format.pgen_hardcall test_pgen.py hardcall integration)
gcta_add_regression_test(format.pgen_dosage test_pgen.py dosage integration)
gcta_add_regression_test(format.mach_dosage test_dosage.py mach integration)
gcta_add_regression_test(format.beagle_dosage test_dosage.py beagle integration)
gcta_add_regression_test(format.mbfile test_multifile.py mbfile integration)
gcta_add_regression_test(format.mbgen test_multifile.py mbgen integration)
gcta_add_regression_test(format.mpfile test_multifile.py mpfile integration)
gcta_add_regression_test(geno.freqx test_freqx.py output integration)
gcta_add_regression_test(error.freqx_no_x test_freqx.py no_x integration)

add_executable(async_buffer_test async_buffer_test.cpp)
target_include_directories(async_buffer_test PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(async_buffer_test Threads::Threads)
add_test(NAME async.buffer COMMAND async_buffer_test)
set_tests_properties(async.buffer PROPERTIES LABELS "unit" TIMEOUT 120)
add_executable(thread_pool_test thread_pool_test.cpp)
target_include_directories(thread_pool_test PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(thread_pool_test threadpool Threads::Threads)
add_test(NAME thread.pool COMMAND thread_pool_test)
set_tests_properties(thread.pool PROPERTIES LABELS "unit" TIMEOUT 120)
gcta_add_regression_test(linkage.numeric_symbols test_linkage.py numeric_symbols integration)
set_tests_properties(linkage.numeric_symbols PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:")
gcta_add_regression_test(backend.grm test_backends.py grm integration)
set_tests_properties(backend.grm PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:")

if(LINUX)
    add_test(
        NAME install.runtime
        COMMAND ${CMAKE_COMMAND}
            -DGCTA_INSTALL_TEST_BINARY_DIR=${CMAKE_BINARY_DIR}
            -DGCTA_INSTALL_TEST_PREFIX=${CMAKE_INSTALL_PREFIX}
            -DGCTA_INSTALL_TEST_WORK=${CMAKE_CURRENT_BINARY_DIR}/output/install.runtime
            -DGCTA_INSTALL_TEST_DATA=${GCTA_TEST_DATA}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/install_smoke.cmake
    )
    set_tests_properties(install.runtime PROPERTIES LABELS "integration" TIMEOUT 120)
endif()
