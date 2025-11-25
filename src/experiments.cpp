#include "basics.hpp"
#include "predefined_automata.hpp"
#include "quantum_program.hpp"
#include "swta.hpp"
#include "swta_builders.hpp"

#include <chrono>

/*
 * RUN_EXPERIMENT 巨集
 * 用途：
 *   執行驗證實驗，自動計時並輸出結果摘要。
 * 怎麼算：
 *   1. 輸出實驗名稱和開始訊息
 *   2. 記錄開始時間
 *   3. 執行實驗函數
 *   4. 記錄結束時間
 *   5. 輸出驗證狀態和執行時間（毫秒）
 */
#define RUN_EXPERIMENT(name, experiment_function) \
{ \
    std::cout << "Running experiment \"" << name << "\"..."; \
    auto begin = std::chrono::steady_clock::now(); \
    auto status = experiment_function(); \
    std::cout << " Done. Summary:\n"; \
    auto end = std::chrono::steady_clock::now(); \
    std::cout << " > Status : " << get_word_for_status(status) << "\n"; \
    std::cout << " > Runtime: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << "\n"; \
}

/*
 * Verification_Status 列舉
 * 用途：
 *   表示驗證實驗的結果狀態。
 * 怎麼算：
 *   FAILURE = 0：驗證失敗
 *   SUCCESS = 1：驗證成功
 */
enum Verification_Status : u32 {
    FAILURE = 0,
    SUCCESS = 1,
};

/*
 * get_word_for_status
 * 用途：
 *   根據驗證狀態返回對應的字串表示。
 * 怎麼算：
 *   若狀態為 FAILURE 則返回 "FAIL"，否則返回 "SUCCESS"。
 */
const char* get_word_for_status(Verification_Status status) {
    return status == Verification_Status::FAILURE ? "FAIL" : "SUCCESS";
}

/*
 * verify_bv
 * 用途：
 *   驗證 Bernstein-Vazirani (BV) 量子算法的正確性。該算法在預設條件上應用 Hadamard 和 CNOT 閘，
 *   然後檢查結果是否與後置條件等價。
 * 怎麼算：
 *   1. 加載 BV 預條件 SWTA
 *   2. 加載 Hadamard 和 Parity CNOT 轉導器
 *   3. 依序應用轉導器到 SWTA
 *   4. 加載 BV 後置條件
 *   5. 使用顏色等價檢查驗證結果與後置條件的等價性
 */
Verification_Status verify_bv() {
    /*
     * 步驟 1 - 對前置條件應用 Hadamard 閘
     */
    SWTA precondition = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_PRE);
    auto metadata = precondition.get_metadata();

    WTT hadamard  = get_predefined_wtt(Predefined_WTT_Names::HADAMARD, metadata);
    WTT parity_cnot = get_predefined_wtt(Predefined_WTT_Names::PARITY_CNOT, metadata);
    std::vector<WTT> required_transducers = {
        hadamard,
        parity_cnot
    };

    SWTA_Program program = {
        .initial_swta = precondition,
        .transducers = required_transducers,
        .applications = {
            Transducer_Application(0),
            Transducer_Application(1),
            Transducer_Application(0),
        }
    };

    SWTA result = run_swta_program(program);
    SWTA postcondition = get_predefined_swta(Predefined_SWTA_Names::BV_EXAMPLE_10STAR_POST);

    bool are_equivalent = are_two_swtas_color_equivalent(result, postcondition);

    return are_equivalent ? Verification_Status::SUCCESS : Verification_Status::FAILURE;
}


/*
 * verify_grover
 * 用途：
 *   驗證 Grover 量子搜尋算法（振幅放大）的正確性。
 *   比較理想實現版本與硬體相容版本是否在語義上等價。
 * 怎麼算：
 *   1. 定義 SWTA 元數據（工作量子位和輔助位）
 *   2. 加載預定義的轉導器（X、H、多控制 Z 閘等）
 *   3. 構建兩個程式：理想版本和硬體版本（使用不同的 CnZ 實現）
 *   4. 依序執行兩個程式
 *   5. 比較結果的顏色等價性
 */
Verification_Status verify_grover() {
    SWTA::Metadata swta_metadata = {
        /*
         * 工作量子位和輔助位
         */
        .number_of_internal_symbols = 3,
        .number_of_colors = 1
    };

    SWTA all_basis = get_predefined_swta(Predefined_SWTA_Names::GROVER_ALL_BASIS);

    std::vector<WTT> needed_transducers {
        /* 0 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_X, swta_metadata),
        /* 1 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_H, swta_metadata),
        /* 2 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_FIRST_MULTI_Z, swta_metadata),
        /* 3 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_FIRST_MULTI_Z_USING_CCX, swta_metadata),
        /* 4 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_SECOND_MULTI_Z, swta_metadata),
        /* 5 */ get_predefined_wtt(Predefined_WTT_Names::GROVER_SECOND_MULTI_Z_USING_CCX, swta_metadata),
    };

    SWTA_Program first_program = {
        .initial_swta = all_basis,
        .transducers = needed_transducers,
        .applications = {
            Transducer_Application(0), // X
            Transducer_Application(2), // CnZ to last ancilla
            Transducer_Application(0), // X
            Transducer_Application(1), // H
            Transducer_Application(0), // X
            Transducer_Application(4), // CnZ to last working qubit
            Transducer_Application(0), // X
            Transducer_Application(1), // H
        }
    };

    SWTA_Program program_for_hw = {
        .initial_swta = all_basis,
        .transducers = needed_transducers,
        .applications = {
            Transducer_Application(0), // X
            Transducer_Application(3), // CnZ to last ancilla
            Transducer_Application(0), // X
            Transducer_Application(1), // H
            Transducer_Application(0), // X
            Transducer_Application(5), // CnZ to last working qubit
            Transducer_Application(0), // X
            Transducer_Application(1), // H
        }
    };

    SWTA result_idealistic = run_swta_program(first_program);
    SWTA result_for_hw     = run_swta_program(program_for_hw);

    bool are_equivalent = are_two_swtas_color_equivalent(
        result_idealistic,
        result_for_hw
    );

    return are_equivalent ? Verification_Status::SUCCESS : Verification_Status::FAILURE;
}

/*
 * verify_adder
 * 用途：
 *   驗證量子加法器電路的正確性。使用梯形建構和順序組合 WTT 來實現加法器。
 * 怎麼算：
 *   1. 加載 MAJ（多數）和 UMA（反向多數）轉導器
 *   2. 順序組合 UMA 轉導器
 *   3. 為 MAJ 和 UMA 執行梯形建構（生成多層結構）
 *   4. 水平組合梯形結構與恆等轉導器
 *   5. 順序組合所有部分形成完整電路
 *   6. 驗證結果與後置條件的等價性
 */
Verification_Status verify_adder() {
    
    SWTA::Metadata metadata = {
        /*
         * 內部符號表示工作量子位
         */
        .number_of_internal_symbols = 1,
        .number_of_colors = 1
    };

    WTT maj = get_predefined_wtt(Predefined_WTT_Names::ADDER_MAJ_RESULT, metadata);

    WTT uma1 = get_predefined_wtt(Predefined_WTT_Names::ADDER_UMA1, metadata);
    WTT uma2 = get_predefined_wtt(Predefined_WTT_Names::ADDER_UMA2, metadata);
    WTT uma3 = get_predefined_wtt(Predefined_WTT_Names::ADDER_UMA3, metadata);

    WTT uma12  = compose_wtts_sequentially(uma1, uma2);
    WTT uma123 = compose_wtts_sequentially(uma12, uma3);

    std::vector<u64> box_inputs {0, 0, 0};
    u64 box_offset = 2;
    u64 new_symbol = 1;

    WTT maj_staircase = perform_staircase_construction(maj,    box_inputs, box_offset, new_symbol, Staircase_Direction::LEFT_RIGHT);
    WTT uma_staircase = perform_staircase_construction(uma123, box_inputs, box_offset, new_symbol, Staircase_Direction::RIGHT_LEFT);

    WTT id1 = get_predefined_wtt(Predefined_WTT_Names::TEST_FIXED_ID1, metadata);
    WTT extended_maj_staircase = compose_wtts_horizontally(maj_staircase, id1);
    WTT extended_uma_staircase = compose_wtts_horizontally(uma_staircase, id1);

    WTT middle_piece = get_predefined_wtt(Predefined_WTT_Names::ADDER_MIDDLE, metadata);

    WTT result12  = compose_wtts_sequentially(extended_maj_staircase, middle_piece);
    WTT adder_circuit = compose_wtts_sequentially(result12, extended_uma_staircase);

    SWTA precondition = get_predefined_swta(Predefined_SWTA_Names::ADDER_PRE);
    SWTA postcondition = get_predefined_swta(Predefined_SWTA_Names::ADDER_POST);

    SWTA result_swta = apply_wtt_to_swta(precondition, adder_circuit);

    auto precondition_dfa = build_frontier_automaton(precondition).determinize();
    precondition_dfa.complete();
    auto result_dfa       = build_frontier_automaton(postcondition).determinize();
    result_dfa.complete();

    bool are_equivalent = are_two_swtas_color_equivalent(result_swta, postcondition);
    return are_equivalent ? Verification_Status::SUCCESS : Verification_Status::FAILURE;
}

/*
 * verify_ecc
 * 用途：
 *   驗證量子糾錯碼 (ECC) 的正確性。測試編碼電路對預定義測試案例的驗證。
 * 怎麼算：
 *   1. 加載 ECC 前置條件和後置條件 SWTA
 *   2. 加載 ECC BOX1 和 BOX2 轉導器
 *   3. 順序組合兩個 BOX 轉導器
 *   4. 執行梯形建構以創建多層 ECC 結構
 *   5. 將梯形結構應用於前置條件
 *   6. 驗證結果與後置條件的等價性
 */
Verification_Status verify_ecc() {
    auto ecc_pre  = get_predefined_swta(Predefined_SWTA_Names::ECC_PRE);
    auto ecc_post  = get_predefined_swta(Predefined_SWTA_Names::ECC_POST);

    SWTA::Metadata metadata = {
        /*
         * 四種顏色表示不同的糾錯碼操作
         */
        .number_of_internal_symbols = 1,
        .number_of_colors = 4
    };

    auto box_part1 = get_predefined_wtt(Predefined_WTT_Names::ECC_BOX1, metadata);
    auto box_part2 = get_predefined_wtt(Predefined_WTT_Names::ECC_BOX2, metadata);
    auto box = compose_wtts_sequentially(box_part1, box_part2);

    u64 new_symbol = 1;
    u64 offset = 2;
    auto staircase = perform_staircase_construction(box, {0, 0, 0, 0}, offset, new_symbol);

    auto ecc_result = apply_wtt_to_swta(ecc_pre, staircase);

    bool are_equivalent = are_two_swtas_color_equivalent(ecc_result, ecc_post);
    return are_equivalent ? Verification_Status::SUCCESS : Verification_Status::FAILURE;
}

/*
 * verify_hamiltonian_simulation
 * 用途：
 *   驗證 Hamiltonian 模擬量子算法的優化正確性。比較天真實現與優化實現的等價性。
 * 怎麼算：
 *   1. 加載所有必要的轉導器（RZZ、RXX、RYY、UZZ、SQRT_X、S、H 等）
 *   2. 為每個轉導器執行梯形建構生成各階段
 *   3. 構建天真程式：依序應用所有轉導器
 *   4. 構建優化程式：僅應用 RXX、RYY、RZZ 的組合
 *   5. 執行兩個程式
 *   6. 驗證結果的顏色等價性
 */
Verification_Status verify_hamiltonian_simulation() {
    using ACN = Algebraic_Complex_Number;

    SWTA::Metadata metadata;

    /*
     * 加載各個 Hamiltonian 模擬階段的轉導器
     */
    WTT rzz_box   = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_RZZ, metadata);
    WTT rzz_stage = perform_staircase_construction(rzz_box, {0, 0}, 1, 1);

    WTT rxx_box   = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_RXX, metadata);
    WTT rxx_stage = perform_staircase_construction(rxx_box, {0, 0}, 1, 1);

    WTT ryy_box   = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_RYY, metadata);
    WTT ryy_stage = perform_staircase_construction(ryy_box, {0, 0}, 1, 1);

    WTT uzz_box   = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_UZZ, metadata);
    WTT uzz_stage = perform_staircase_construction(uzz_box, {0, 0}, 1, 1);

    WTT sqrt_x_stage = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_SQRT_X_STAGE, metadata);
    WTT s_stage = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_S_STAGE, metadata);
    WTT h_stage = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_H_STAGE, metadata);
    WTT last_x_stage = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_LAST_X_STAGE, metadata);

    std::vector<WTT> circuit_stages = {
        h_stage,       /*
                        * 0 - Hadamard 階段
                        */
        rzz_stage,     /*
                        * 1 - RZZ 旋轉階段
                        */
        sqrt_x_stage,  /*
                        * 2 - sqrt(X) 階段
                        */
        uzz_stage,     /*
                        * 3 - UZZ 階段
                        */
        last_x_stage,  /*
                        * 4 - 最後的 X 階段
                        */
        s_stage,       /*
                        * 5 - S 閘階段
                        */
        rxx_stage,     /*
                        * 6 - RXX 旋轉階段
                        */
        ryy_stage,     /*
                        * 7 - RYY 旋轉階段
                        */
    };

    SWTA precondition = get_predefined_swta(Predefined_SWTA_Names::HAMILTONIAN_ALL_BASIS);

    SWTA_Program naive_program = {
        .initial_swta = precondition,
        .transducers = circuit_stages,
        .applications = {
            /*
             * 天真實現：依序應用所有階段
             */
            Transducer_Application(0),
            Transducer_Application(1),
            Transducer_Application(2),
            Transducer_Application(3),
            Transducer_Application(4),
            Transducer_Application(0),
            Transducer_Application(5),
            Transducer_Application(1),
        },
    };

    SWTA_Program optimized_program = {
        .initial_swta = precondition,
        .transducers = circuit_stages,
        .applications = {
            /*
             * 優化實現：僅應用 RXX、RYY、RZZ 的組合
             */
            Transducer_Application(6),
            Transducer_Application(7),
            Transducer_Application(1),
        },
    };

    auto naive_program_result     = run_swta_program(naive_program);
    auto optimized_program_result = run_swta_program(optimized_program);

    bool are_equivalent = are_two_swtas_color_equivalent(naive_program_result, optimized_program_result);
    return are_equivalent ? Verification_Status::SUCCESS : Verification_Status::FAILURE;
}


/*
 * main
 * 用途：
 *   程式入口點。執行所有量子算法驗證實驗。
 * 怎麼算：
 *   1. 執行 Bernstein-Vazirani 驗證
 *   2. 執行 Grover 算法驗證
 *   3. 執行量子加法器驗證
 *   4. 執行糾錯碼驗證
 *   5. 執行 Hamiltonian 模擬驗證
 *   6. 返回 0 表示程式正常結束
 */
int main() {

    RUN_EXPERIMENT("Bernstein-Vazirani", verify_bv);
    RUN_EXPERIMENT("Grover (amplitude amplification)", verify_grover);
    RUN_EXPERIMENT("Adder", verify_adder);
    RUN_EXPERIMENT("Error correction code", verify_ecc);
    RUN_EXPERIMENT("Hamiltonian simulation", verify_hamiltonian_simulation);

    return 0;
}
