#pragma once

#include "arith.hpp"
#include "swta.hpp"


enum class Predefined_WTT_Names : u64 {
    HADAMARD    = 0,
    PARITY_CNOT = 1,

    TEST_STAIRCASE_IDENTITY3 = 2,

    GROVER_FIRST_MULTI_Z  =  8,
    GROVER_SECOND_MULTI_Z =  9,
    GROVER_X              = 10,
    GROVER_H              = 11,
    GROVER_FIRST_MULTI_Z_USING_CCX  = 12,
    GROVER_SECOND_MULTI_Z_USING_CCX = 13,

    ADDER_UMA1 = 20,
    ADDER_UMA2 = 21,
    ADDER_UMA3 = 22,
    ADDER_UMA_RESULT = 23, // 未定義
    ADDER_MAJ1 = 24,
    ADDER_MAJ2 = 21, // same as ADDER_UMA2
    ADDER_MAJ3 = 20, // same as ADDER_UMA1
    ADDER_MAJ_RESULT = 30,
    ADDER_MAJ_RESULT_12 = 31, // 第一階段與第二階段的組合
    ADDER_MIDDLE = 32,

    ECC_BOX1 = 40,
    ECC_BOX2 = 41,

    HAMILTONIAN_BC_CNOT = 50,
    HAMILTONIAN_BC_RZ   = 51,
    HAMILTONIAN_BC_H    = 52,
    HAMILTONIAN_BC_S    = 53,
    HAMILTONIAN_BC_X    = 54,

    HAMILTONIAN_RZZ = 60,
    HAMILTONIAN_RXX = 61,
    HAMILTONIAN_RYY = 62,
    HAMILTONIAN_UZZ = 63,

    HAMILTONIAN_H_STAGE      = 64,
    HAMILTONIAN_S_STAGE      = 65,
    HAMILTONIAN_SQRT_X_STAGE = 66,
    HAMILTONIAN_LAST_X_STAGE = 67,

    TEST_FIXED_ID1 = 100
};

enum class Predefined_SWTA_Names : u64 {
    BV_EXAMPLE_10STAR_PRE = 0,
    BV_EXAMPLE_10STAR_POST = 1,
    BV_EXAMPLE_10STAR_RESULT = 2,

    TRIVIAL_BOT = 3,   // 用於測試
    TRIVIAL_ONES = 4,  // 用於測試
    TEST_BV_EXAMPLE_AFTER_STEP1 = 5,
    TEST_BV_EXAMPLE_AFTER_STEP2 = 6,
    TEST_BV_EXAMPLE_AFTER_STEP3 = 7,


    GROVER_ALL_BASIS = 8,

    TEST_ADDER_ALL_3BASIS = 10,
    ADDER_PRE        = 11,
    ADDER_POST       = 12,

    ECC_PRE        = 20,
    ECC_POST       = 21,

    HAMILTONIAN_ALL_BASIS = 30,
};

struct Def_State;
struct Def_Coef;
struct Def_Linear_Form;

/**
 * 封裝 WTT/SWTA 的狀態以便在建構線性形式或轉移時直接引用，避免在不同物件之間傳遞裸 State。
 * 主要用途是讓組合線性形式時可讀性更佳，計算時僅攜帶原始狀態值，不影響狀態編號。
 */
struct Def_State {
    State state;

    // 封裝狀態編號，直接儲存給定的狀態值，方便在建構線性形式時引用。
    Def_State(State state) : state(state) {}
};

struct Def_Coef {
    using ACN = Algebraic_Complex_Number;

    ACN number;

    // 初始化複數係數，後續可與狀態相乘生成線性形式項。
    Def_Coef(const ACN& number) : number(number) {}

    // 將係數與狀態相乘，產生帶有該係數與狀態的線性形式元素。
    Def_Linear_Form operator*(const Def_State& other);
};

struct Def_Linear_Form {
    using ACN = Algebraic_Complex_Number;

    ACN         coef;
    State       state;
    Subtree_Tag tag;

    // 以係數與狀態建構線性形式，並預設不帶子樹標記。
    Def_Linear_Form(const Def_Coef& def_coef, const Def_State& def_state) : coef(def_coef.number), state(def_state.state), tag(Subtree_Tag::NONE) {};

    // 複製另一個線性形式，保留係數、狀態與子樹標記。
    Def_Linear_Form(const Def_Linear_Form& other) : coef(other.coef), state(other.state), tag(other.tag) {};

    // 為線性形式標記來源子樹方向（左或右），以便合成轉移時過濾。
    Def_Linear_Form operator*(const Subtree_Tag& tag) {
        this->tag = tag;
        return *this;
    }
};



// 依照預設名稱生成對應的 WTT：會根據傳入的名稱挑選對應的建構流程，
// 使用中繼資料中的內部符號與顏色數量建立 WTT_Builder，
// 並逐步新增線性形式描述的轉移、初始與終止狀態，最後回傳完成的樹轉換器。
WTT get_predefined_wtt(Predefined_WTT_Names name, const SWTA::Metadata& metadata);

// 依照預設名稱生成對應的 SWTA：依需求配置內部符號數、顏色數及轉移函式，
// 並標記初始與接受葉節點集合，使呼叫端能直接取得可執行的樹自動機。
SWTA get_predefined_swta(Predefined_SWTA_Names name);

// 根據標記好的線性形式集合，將左右子樹的分量收集起來並組成 SWTA 轉移：
// 會掃描左右集合中子樹標記為 NONE 的項目，將係數與狀態打包成 Linear_Form，
// 進而產生 SWTA 轉移的左右線性組合 (f_left, f_right)。
SWTA::Transition synthetize_swta_transition(const std::vector<Def_Linear_Form>& left_subtree, const std::vector<Def_Linear_Form>& right_subtree);

// 根據標記好的線性形式集合，收集四個子樹分量並建立 WTT 轉移（LL、LR、RL、RR）：
// 會將左右集合中標記 LEFT/RIGHT 的項目分別整理成四份 Linear_Form，
// 以符合 WTT 轉移矩陣對左／右子樹的映射要求。
WTT::Transition synthetize_wtt_transition(const std::vector<Def_Linear_Form>& left_subtree, const std::vector<Def_Linear_Form>& right_subtree);
