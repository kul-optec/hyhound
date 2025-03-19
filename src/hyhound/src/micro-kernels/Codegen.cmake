function(codegen_hyhound_microkernels tgt)

    set(OUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/${tgt}/codegen")
    set(TPP_DIR "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../../include/hyhound/micro-kernels")
    foreach(alg "householder-updowndate")
        target_precompile_headers(${tgt} PRIVATE "${TPP_DIR}/../${alg}-serial.tpp")
        foreach (op "diag" "full" "tail")
            target_precompile_headers(${tgt} PRIVATE "${TPP_DIR}/${alg}-${op}.tpp")
        endforeach()
        foreach(R 32 16 12 8 4 2 1)
            if (R GREATER HYHOUND_MAX_HYH_KERNEL_WIDTH)
                continue()
            endif()
            configure_file("${CMAKE_CURRENT_FUNCTION_LIST_DIR}/${alg}-diag.cpp.in" "${OUT_DIR}/${alg}-diag-${R}.cpp" @ONLY)
            target_sources(${tgt} PRIVATE "${OUT_DIR}/${alg}-diag-${R}.cpp")
            foreach(S 32 24 16 12 8 4 2 1)
                if (S GREATER HYHOUND_MAX_HYH_KERNEL_HEIGHT)
                    continue()
                endif()
                configure_file("${CMAKE_CURRENT_FUNCTION_LIST_DIR}/${alg}-tail.cpp.in" "${OUT_DIR}/${alg}-tail-${R}-${S}.cpp" @ONLY)
                target_sources(${tgt} PRIVATE "${OUT_DIR}/${alg}-tail-${R}-${S}.cpp")
                configure_file("${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../${alg}-serial.cpp.in" "${OUT_DIR}/${alg}-serial-${R}-${S}.cpp" @ONLY)
                target_sources(${tgt} PRIVATE "${OUT_DIR}/${alg}-serial-${R}-${S}.cpp")
            endforeach()
        endforeach()
        foreach(R RANGE 1 32)
            if (R GREATER HYHOUND_MAX_HYH_KERNEL_WIDTH)
                continue()
            endif()
            configure_file("${CMAKE_CURRENT_FUNCTION_LIST_DIR}/${alg}-full.cpp.in" "${OUT_DIR}/${alg}-full-${R}.cpp" @ONLY)
            target_sources(${tgt} PRIVATE "${OUT_DIR}/${alg}-full-${R}.cpp")
        endforeach()
        foreach(R 15 14 13 12 11 10 9 7 6 5 3)
            if (R GREATER HYHOUND_MAX_HYH_KERNEL_WIDTH)
                continue()
            endif()
            configure_file("${CMAKE_CURRENT_FUNCTION_LIST_DIR}/${alg}-diag.cpp.in" "${OUT_DIR}/${alg}-diag-${R}.cpp" @ONLY)
            target_sources(${tgt} PRIVATE "${OUT_DIR}/${alg}-diag-${R}.cpp")
            foreach(S 32 24 16 12 8 4 2 1)
                if (S GREATER HYHOUND_MAX_HYH_KERNEL_HEIGHT)
                    continue()
                endif()
                configure_file("${CMAKE_CURRENT_FUNCTION_LIST_DIR}/${alg}-tail.cpp.in" "${OUT_DIR}/${alg}-tail-${R}-${S}.cpp" @ONLY)
                target_sources(${tgt} PRIVATE "${OUT_DIR}/${alg}-tail-${R}-${S}.cpp")
            endforeach()
        endforeach()
    endforeach()

endfunction()
