#include "predefined_automata.hpp"
#include "arith.hpp"
#include "basics.hpp"
#include "swta.hpp"

#include <stdexcept>
#include <string>

/**
 * 輔助結構：在多重預設自動機中追蹤兩種「身分」狀態，
 * 以便在配置轉移時快速辨識是否落在等價的 Id 分支。
 * 建構時計算量為 O(1)，僅儲存狀態與預先分好的線性形式列表，
 * 使後續判斷或收集轉移時不需要重新排列資料。
 */
struct Alternative_Id_Helper {
    State ida;
    State idb;
    std::vector<Def_Linear_Form> ida_branch;
    std::vector<Def_Linear_Form> idb_branch;
    std::vector<Def_Linear_Form> ida_fin_branch;
    std::vector<Def_Linear_Form> idb_fin_branch;

    // 檢查給定狀態是否為替代的身分狀態：直接比較 state 是否等於任一紀錄的 Id，時間複雜度 O(1)。
    bool is_id(State state) const {
        return state == ida || state == idb;
    }
};

// 將係數包裝成可與狀態相乘的線性形式：回傳一個帶有給定係數與狀態的 Def_Linear_Form，供後續標記子樹或累加。
Def_Linear_Form Def_Coef::operator*(const Def_State& other) {
    return Def_Linear_Form(*this, other);
}

// 收集指定子樹標記的線性形式分量：篩選 components 中 tag 符合要求的項目，
// 將其係數與狀態放入 Linear_Form::components 中；未匹配的項目會被忽略，
// 形成對應轉移矩陣的單一線性組合。
Linear_Form collect_tagged_subtrees_into_form(const std::vector<Def_Linear_Form>& components, Subtree_Tag tag) {
    Linear_Form form;
    for (auto& component : components) {
        if (component.tag == tag) {
            form.components.push_back({component.coef, component.state});
        }
    }
    return form;
}

// 從左右子樹的線性形式集合中分別收集 LL/LR 與 RL/RR 部分：
// 先依標記整理出左子樹的 LEFT/RIGHT 與右子樹的 LEFT/RIGHT 線性形式，
// 再以這四個 Linear_Form 初始化 WTT::Transition，對應 (LL, LR, RL, RR) 四個象限的振幅。
WTT::Transition synthetize_wtt_transition(const std::vector<Def_Linear_Form>& left_subtree, const std::vector<Def_Linear_Form>& right_subtree) {
    Linear_Form ll = collect_tagged_subtrees_into_form(left_subtree, Subtree_Tag::LEFT);
    Linear_Form lr = collect_tagged_subtrees_into_form(left_subtree, Subtree_Tag::RIGHT);
    Linear_Form rl = collect_tagged_subtrees_into_form(right_subtree, Subtree_Tag::LEFT);
    Linear_Form rr = collect_tagged_subtrees_into_form(right_subtree, Subtree_Tag::RIGHT);

    return WTT::Transition(ll, lr, rl, rr);
}

// 從左右子樹的線性形式集合中收集未標記的項目：
// 僅挑選 tag 為 NONE 的項目填入左、右 Linear_Form，
// 產生 SWTA::Transition 需要的 (f_left, f_right) 兩組幅度分佈。
SWTA::Transition synthetize_swta_transition(const std::vector<Def_Linear_Form>& left_subtree, const std::vector<Def_Linear_Form>& right_subtree) {
    Linear_Form ll = collect_tagged_subtrees_into_form(left_subtree,  Subtree_Tag::NONE);
    Linear_Form rr = collect_tagged_subtrees_into_form(right_subtree, Subtree_Tag::NONE);

    return SWTA::Transition(ll, rr);
}

// 依據給定的四個複數幅度構造哈密頓量階段的 WTT：
// 會建立三個狀態（應用、終止前、葉），並以傳入的 ll/lr/rl/rr 係數填入轉移的四個象限，
// 使工作量子位在掃描過程中套用對應幅度，最終在葉節點收束成給定的幅度矩陣。
WTT construct_stage_for_hamiltonian(const Algebraic_Complex_Number& ll, const Algebraic_Complex_Number& lr, const Algebraic_Complex_Number& rl, const Algebraic_Complex_Number& rr) {
    using DLF = std::vector<Def_Linear_Form>;

    State q_apply      = 0;
    State q_apply_last = 1;
    State q_leaf       = 2;

    u64 number_of_states = 3;

    Internal_Symbol working_qubit = 0;
    Internal_Symbol working_qubit_stop = 1;

    SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
    WTT_Builder builder (metadata);

    builder.mark_state_initial(q_apply);
    builder.mark_state_final(q_leaf);

    Def_Coef def_ll(ll);
    Def_Coef def_lr(lr);
    Def_Coef def_rl(rl);
    Def_Coef def_rr(rr);

    auto add_single_target_transition = [&builder, &def_ll, &def_lr, &def_rl, &def_rr](State source, State target, Internal_Symbol internal_sym) {
         DLF left_subtree  { def_ll * target * Subtree_Tag::LEFT };
         if (!def_lr.number.is_zero()) {
             left_subtree.push_back( def_lr * target * Subtree_Tag::RIGHT );
         }

         DLF right_subtree { def_rr * target * Subtree_Tag::RIGHT };
         if (!def_rl.number.is_zero()) {
             right_subtree.push_back( def_rl * target * Subtree_Tag::LEFT );
         }

         auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
         builder.add_transition(source, internal_sym, transition);
    };

    add_single_target_transition(q_apply, q_apply, working_qubit);
    add_single_target_transition(q_apply, q_apply_last, working_qubit_stop);
    add_single_target_transition(q_apply_last, q_leaf, working_qubit);

    WTT result = builder.build(number_of_states);

    if (DEBUG) {
        result.debug_data = new WTT::Debug_Data;
        auto& state_names = result.debug_data->state_names;

        state_names[q_apply]      = "apply";
        state_names[q_apply_last] = "apply.last";
        state_names[q_leaf]       = "leaf";
    }

    return result;
}

// 依照預設的名稱組裝對應的 WTT：透過名稱選擇對應的建構模板，
// 在每個模板中設定狀態數、內部符號、初末狀態，並以線性形式描述各內部符號的振幅流向，
// 將複數係數分配到左右子樹，使得 WTT 能模擬相應的量子門或控制流程。
WTT get_predefined_wtt(Predefined_WTT_Names name, const SWTA::Metadata& swta_metadata) {
    using DLF = std::vector<Def_Linear_Form>;
    using ACN = Algebraic_Complex_Number;

    if (name == Predefined_WTT_Names::HADAMARD) {
        // q0 -> LEFT{ 1/sqrt(2)(q0, L) + 1/sqrt(2)(q0, R) }, RIGHT{ 1/sqrt(2)(q0, L) - 1/sqrt(2)(q0, R) }
        // q0(left) -> (left)
        
        State q0 = 0;

        std::vector<std::vector<WTT::Transition>> transitions_by_symbol;
        transitions_by_symbol.resize(1);

        Def_Coef one_over_sqrt2       ( Algebraic_Complex_Number::ONE_OVER_SQRT2());
        Def_Coef minus_one_over_sqrt2 (-Algebraic_Complex_Number::ONE_OVER_SQRT2());

        for (Internal_Symbol symbol = 0; symbol < swta_metadata.number_of_internal_symbols; symbol++) {
            DLF left_subtree  {one_over_sqrt2 * q0 * Subtree_Tag::LEFT, one_over_sqrt2 * q0 * Subtree_Tag::RIGHT};
            DLF right_subtree {one_over_sqrt2 * q0 * Subtree_Tag::LEFT, minus_one_over_sqrt2 * q0 * Subtree_Tag::RIGHT};
            WTT::Transition transition = synthetize_wtt_transition(left_subtree, right_subtree);

            transitions_by_symbol[0].push_back(transition);
        }

        WTT transducer (transitions_by_symbol, {0}, {0});
        return transducer;
    }

    if (name == Predefined_WTT_Names::PARITY_CNOT) {
        // 依照 Bernstein-Vazirani 電路的秘密字串 (10)* 計算位元奇偶性。
        // 量子位/狀態以 A、B 標記。
        State a0 = 0; // 奇數位量子位，偶數奇偶性
        State b0 = 1; // 偶數位量子位，偶數奇偶性
        State a1 = 2; // 奇數位量子位，奇數奇偶性
        State b1 = 3; // 偶數位量子位，奇數奇偶性

        size_t state_cnt = 4;
        size_t internal_symbol_cnt = 2;

        Internal_Symbol work_qubit = 0, ancilla = 1;

        std::vector<std::vector<WTT::Transition>> transitions;
        transitions.resize(4);
        for (State state = 0; state < state_cnt; state++) {
            transitions[state].resize(internal_symbol_cnt);
        }

        // 移動：a0 -w-> b0(L), b1(R)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), b0);
            Linear_Form ll ({ll_component});

            Linear_Form::Component rr_component (Algebraic_Complex_Number::ONE(), b1);
            Linear_Form rr ({rr_component});

            WTT::Transition transition (ll, {}, {}, rr);

            transitions[a0][work_qubit] = transition;
        }

        // 移動：b0 -w-> a0(L), a0(R)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), a0);
            Linear_Form ll ({ll_component});

            WTT::Transition transition (ll, {}, {}, ll);

            transitions[b0][work_qubit] = transition;
        }

        // 移動：a1 -w-> b1(L), b0(R)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), b1);
            Linear_Form ll ({ll_component});

            Linear_Form::Component rr_component (Algebraic_Complex_Number::ONE(), b0);
            Linear_Form rr ({rr_component});

            WTT::Transition transition (ll, {}, {}, rr);

            transitions[a1][work_qubit] = transition;
        }

        // 移動：b1 -w-> a1(L), a1(R)
        {
            Linear_Form::Component ll_component (Algebraic_Complex_Number::ONE(), a1);
            Linear_Form ll ({ll_component});

            WTT::Transition transition (ll, {}, {}, ll);

            transitions[b1][work_qubit] = transition;
        }

        // 移動： a0 -a-> b0(L), b0(r)
        // 移動： b0 -a-> b0(L), b0(r)
        {
            std::vector<Def_Linear_Form> left_subtree  {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::LEFT};
            std::vector<Def_Linear_Form> right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::RIGHT};
            WTT::Transition transition = synthetize_wtt_transition(left_subtree, right_subtree);

            transitions[a0][ancilla] = transition;
            transitions[b0][ancilla] = transition;
        }

        // 移動： a1 -a-> b0(R), b0(L)
        // 移動： b1 -a-> b0(R), b0(L)
        {
            std::vector<Def_Linear_Form> left_subtree  {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::RIGHT};
            std::vector<Def_Linear_Form> right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * b0 * Subtree_Tag::LEFT};
            WTT::Transition transition = synthetize_wtt_transition(left_subtree, right_subtree);

            transitions[a1][ancilla] = transition;
            transitions[b1][ancilla] = transition;
        }

        WTT transducer (transitions, {b0}, {a0});
        return transducer;
    }

    if (name == Predefined_WTT_Names::GROVER_FIRST_MULTI_Z) {
        State q_init         = 0;
        State q_w            = 1;
        State q_a            = 2;
        State q_last_ancilla = 3;
        State q_id           = 4;
        State q_leaf         = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);
        builder.mark_state_final(q_id);

        Def_Coef one (ACN::ONE());
        { // q_init --w--> ( Id(L), q_w(R) )
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_w  * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_w --w--> ( Id(L), q_a(R) )
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_a  * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_w, working_qubit, transition);
        }
        { // q_w --W'--> ( Id(L), q_last_ancilla(R) )
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_last_ancilla  * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_w, last_working_qubit, transition);
        }

        { // q_a --a--> ( q_w(L), q_w(R) )
             DLF left_subtree  {one * q_w * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_w * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_a, ancilla, transition);
        }

        // q_id --w-->  ( q_id(L), q_id(R) )
        // q_id --W'--> ( q_id(L), q_id(R) )
        // q_id --a-->  ( q_id(L), q_id(R) )
        {
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_id * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);

             builder.add_transition(q_id, working_qubit, transition);
             builder.add_transition(q_id, last_working_qubit, transition);
             builder.add_transition(q_id, ancilla, transition);
        }

        // q_last_ancilla --a-->  (q_leaf(L), -1 q_leaf(R) )
        {
             DLF left_subtree  {one * q_leaf * Subtree_Tag::LEFT};
             DLF right_subtree {Def_Coef(-ACN::ONE()) * q_leaf * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);

             builder.add_transition(q_last_ancilla, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::GROVER_SECOND_MULTI_Z) {
        State q_init         = 0;
        State q_w            = 1;
        State q_a            = 2;
        State q_last_ancilla = 3;
        State q_id           = 4;
        State q_leaf         = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);
        builder.mark_state_final(q_id);

        Def_Coef one (ACN::ONE());
        { // q_init --w--> ( Id(L), q_w(R) )
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_w  * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_w --w--> ( Id(L), q_a(R) )
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_a  * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_w, working_qubit, transition);
        }
        { // q_w --W'--> ( 1 * q_last_ancilla(L), -1 q_last_ancilla(R) )
             DLF left_subtree  {                  one * q_last_ancilla * Subtree_Tag::LEFT};
             DLF right_subtree {Def_Coef(-ACN::ONE()) * q_last_ancilla * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_w, last_working_qubit, transition);
        }

        { // q_a --a--> ( q_w(L), q_w(R) )
             DLF left_subtree  {one * q_w * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_w * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_a, ancilla, transition);
        }

        // q_id --w-->  ( q_id(L), q_id(R) )
        // q_id --W'--> ( q_id(L), q_id(R) )
        // q_id --a-->  ( q_id(L), q_id(R) )
        {
             DLF left_subtree  {one * q_id * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_id * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);

             builder.add_transition(q_id, working_qubit, transition);
             builder.add_transition(q_id, last_working_qubit, transition);
             builder.add_transition(q_id, ancilla, transition);
        }

        // q_last_ancilla --a-->  (q_leaf(L), q_leaf(R) )
        {
             DLF left_subtree  {one * q_leaf * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_leaf * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);

             builder.add_transition(q_last_ancilla, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::GROVER_X) {
        State q_init         = 0;
        State q_x            = 1;
        State q_a            = 2;
        State q_last_ancilla = 3;
        State q_leaf         = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);

        Def_Coef one (ACN::ONE());
        { // q_init --w--> ( q_x(R), q_x(L) )
             DLF left_subtree  {one * q_x * Subtree_Tag::RIGHT };
             DLF right_subtree {one * q_x * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_x --w--> ( q_a(R), q_a(L) )
             DLF left_subtree  {one * q_a * Subtree_Tag::RIGHT };
             DLF right_subtree {one * q_a * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_x, working_qubit, transition);
        }

        { // q_x --W'--> ( q_last_ancilla(R), q_last_ancilla(L) )
             DLF left_subtree  {one * q_last_ancilla * Subtree_Tag::RIGHT };
             DLF right_subtree {one * q_last_ancilla * Subtree_Tag::LEFT  };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_x, last_working_qubit, transition);
        }

        { // q_a --a--> ( q_x(L), q_x(R) )
             DLF left_subtree  {one * q_x * Subtree_Tag::LEFT  };
             DLF right_subtree {one * q_x * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_a, ancilla, transition);
        }

        { // q_last_ancilla --a--> ( q_leaf(L), q_leafw(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_last_ancilla, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::GROVER_H) {
        State q_init         = 0;
        State q_h            = 1;
        State q_a            = 2;
        State q_last_ancilla = 3;
        State q_leaf         = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);

        Def_Coef omega       ( ACN::ONE_OVER_SQRT2());
        Def_Coef minus_omega (-ACN::ONE_OVER_SQRT2());
        Def_Coef one         (ACN::ONE());

        { // q_init --w--> ( omega*q_h(L) + omega*q_h(R), omega*q_h(L) - omega*q_h(R) )
             DLF left_subtree  { omega * q_h * Subtree_Tag::LEFT, omega * q_h * Subtree_Tag::RIGHT};
             DLF right_subtree { omega * q_h * Subtree_Tag::LEFT, minus_omega * q_h * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_h --w--> ( omega*q_a(L) + omega*q_a(R), omega*q_a(L) - omega*q_a(R) )
             DLF left_subtree  { omega * q_a * Subtree_Tag::LEFT, omega * q_a * Subtree_Tag::RIGHT};
             DLF right_subtree { omega * q_a * Subtree_Tag::LEFT, minus_omega * q_a * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_h, working_qubit, transition);
        }

        { // q_h --W'--> ( omega*q_last_ancilla(L) + omega*q_last_ancilla(R), omega*q_last_ancilla(L) - omega*q_last_ancilla(R) )
             DLF left_subtree  { omega * q_last_ancilla * Subtree_Tag::LEFT, omega * q_last_ancilla * Subtree_Tag::RIGHT};
             DLF right_subtree { omega * q_last_ancilla * Subtree_Tag::LEFT, minus_omega * q_last_ancilla * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_h, last_working_qubit, transition);
        }

        { // q_a --a--> ( q_h(L), q_h(R) )
             DLF left_subtree  { one * q_h * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_h * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_a, ancilla, transition);
        }

        { // q_last_ancilla --a--> ( q_leaf(L), q_leaf(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_last_ancilla, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::GROVER_FIRST_MULTI_Z_USING_CCX) {
        State q_init         = 0;
        State q_no_swap      = 1;
        State q_maybe_swap   = 2;
        State q_do_swap      = 3;
        State q_no_apply     = 4;
        State q_apply        = 5;
        State q_leaf         = 6;

        u64 number_of_states = 7;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);

        Def_Coef one       ( ACN::ONE());
        Def_Coef minus_one (-ACN::ONE());

        { // q_init --w--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  {one * q_no_swap    * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_maybe_swap * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_no_swap --w--> ( q_no_swap(L), q_no_swap(R) )
             DLF left_subtree  {one * q_no_swap * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_no_swap * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, working_qubit, transition);
        }

        { // q_no_swap --a--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  { one * q_no_swap    * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_maybe_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, ancilla, transition);
        }

        { // q_maybe_swap --w--> ( q_no_swap(L), q_swap(R) )
             DLF left_subtree  { one * q_no_swap * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_do_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap, working_qubit, transition);
        }

        { // q_swap --a--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  { one * q_no_swap    * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_maybe_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_do_swap, ancilla, transition);
        }

        { // q_no_swap --W'--> ( q_no_apply(L), q_no_apply(R) )
             DLF left_subtree  { one * q_no_apply * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_no_apply * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, last_working_qubit, transition);
        }

        { // q_maybe_swap --W'--> ( q_no_apply(L), q_apply(R) )
             DLF left_subtree  { one * q_no_apply * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_apply * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap, last_working_qubit, transition);
        }

        { // q_no_apply --a--> ( leaf(L), leaf(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_apply, ancilla, transition);
        }

        { // q_apply --a--> ( leaf(L), -leaf(R) )
             DLF left_subtree  {       one * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { minus_one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_apply, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::GROVER_SECOND_MULTI_Z_USING_CCX) {
        State q_init         = 0;
        State q_no_swap      = 1;
        State q_maybe_swap   = 2;
        State q_do_swap      = 3;
        State q_last_anc     = 4;
        State q_leaf         = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_init);
        builder.mark_state_final(q_leaf);

        Def_Coef one       ( ACN::ONE());
        Def_Coef minus_one (-ACN::ONE());

        { // q_init --w--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  {one * q_no_swap    * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_maybe_swap * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_init, working_qubit, transition);
        }

        { // q_no_swap --w--> ( q_no_swap(L), q_no_swap(R) )
             DLF left_subtree  {one * q_no_swap * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_no_swap * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, working_qubit, transition);
        }

        { // q_no_swap --a--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  { one * q_no_swap    * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_maybe_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, ancilla, transition);
        }

        { // q_maybe_swap --w--> ( q_no_swap(L), q_swap(R) )
             DLF left_subtree  { one * q_no_swap * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_do_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap, working_qubit, transition);
        }

        { // q_swap --a--> ( q_no_swap(L), q_maybe_swap(R) )
             DLF left_subtree  { one * q_no_swap    * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_maybe_swap * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_do_swap, ancilla, transition);
        }

        { // q_no_swap --W'--> ( q_last_anc(L), q_last_anc(R) )
             DLF left_subtree  { one * q_last_anc * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_last_anc * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, last_working_qubit, transition);
        }

        { // q_maybe_swap --W'--> ( q_last_anc(L), -q_last_anc(R) )
             DLF left_subtree  {       one * q_last_anc * Subtree_Tag::LEFT };
             DLF right_subtree { minus_one * q_last_anc * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap, last_working_qubit, transition);
        }

        { // q_last_anc --a--> ( leaf(L), leaf(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_last_anc, ancilla, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::TEST_STAIRCASE_IDENTITY3) {
        State q3     = 0;
        State q2     = 1;
        State q1     = 2;
        State q_leaf = 3;

        u64 number_of_states = 4;

        Internal_Symbol working_qubit = 0, ancilla = 1;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q3);
        builder.mark_state_final(q_leaf);

        Def_Coef one (  ACN::ONE());
        Def_Coef two (ACN(2, 0, 0, 0, 0));

        { // q3 --w--> ( q2(L), q2(R) )
             DLF left_subtree  {two * q2 * Subtree_Tag::LEFT};
             DLF right_subtree {two * q2 * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q3, working_qubit, transition);
        }

        { // q2 --a--> ( q1(L), q1(R) )
             DLF left_subtree  {two * q1 * Subtree_Tag::LEFT};
             DLF right_subtree {two * q1 * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q2, ancilla, transition);
        }

        { // q1 --w--> ( q_leaf(L), q_leaf(R) )
             DLF left_subtree  {one * q_leaf * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_leaf * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q1, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_UMA1 || name == Predefined_WTT_Names::ADDER_MAJ3) {
        State q_maybe_swap0 = 0;
        State q_maybe_swap1 = 1;
        State q_no_swap1    = 2;
        State q_no_swap2    = 3;
        State q_swap        = 4;
        State q_leaf        = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_maybe_swap0);
        builder.mark_state_final(q_leaf);

        Def_Coef one (  ACN::ONE());

        { // q_maybe_swap0 --w--> ( q_no_swap1(L), q_maybe_swap1(R) )
             DLF left_subtree  {one * q_no_swap1    * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_maybe_swap1 * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap0, working_qubit, transition);
        }

        { // q_no_swap1 --w--> ( q_no_swap2(L), q_no_swap2(R) )
             DLF left_subtree  {one * q_no_swap2 * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_no_swap2 * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap1, working_qubit, transition);
        }

        { // q_no_swap2 --w--> ( q_leaf(L), q_leaf(R) )
             DLF left_subtree  {one * q_leaf * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_leaf * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap2, working_qubit, transition);
        }

        { // q_maybe_swap1 --w--> ( q_no_swap2(L), q_swap(R) )
             DLF left_subtree  {one * q_no_swap2 * Subtree_Tag::LEFT};
             DLF right_subtree {one * q_swap     * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_maybe_swap1, working_qubit, transition);
        }

        { // q_swap --w--> ( q_leaf(R), q_leaf(L) )
             DLF left_subtree  {one * q_leaf * Subtree_Tag::RIGHT};
             DLF right_subtree {one * q_leaf * Subtree_Tag::LEFT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_swap, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_UMA2 || name == Predefined_WTT_Names::ADDER_MAJ2) {
        State root = 0;
        State p_star1  = 1; // 投影掉所有帶「1」的分支
        State p_1      = 2;
        State e_star1  = 3; // 萃取並保留所有帶「1」的分支值
        State e_1      = 4;
        State q_leaf   = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(root);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // root --w--> ( p_star1(L) + e_star1(R), p_star1(R) + e_star1(L) )
             DLF left_subtree  { one * p_star1 * Subtree_Tag::LEFT, one * e_star1 * Subtree_Tag::RIGHT };
             DLF right_subtree { one * p_star1 * Subtree_Tag::RIGHT, one * e_star1 * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(root, working_qubit, transition);
        }

        { // p_star1 --w--> ( p_1(L), p_1(R))
             DLF left_subtree  { one * p_1 * Subtree_Tag::LEFT };
             DLF right_subtree { one * p_1 * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(p_star1, working_qubit, transition);
        }

        { // p_1 --w--> ( q_leaf(L), 0 q_leaf(R) )
             DLF left_subtree  { one  * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { zero * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(p_1, working_qubit, transition);
        }

        { // e_star1 --w--> ( e_1(L), e_1(R))
             DLF left_subtree  { one * e_1 * Subtree_Tag::LEFT };
             DLF right_subtree { one * e_1 * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(e_star1, working_qubit, transition);
        }

        { // e_1 --w--> ( 0 q_leaf(L), q_leaf(R))
             DLF left_subtree  { zero * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { one  * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(e_1, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_UMA3) {
        State root      = 0;
        State q_no_swap = 1;
        State q_swap    = 2;
        State q_id_1    = 3;
        State q_leaf    = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(root);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());

        { // root --w--> ( p_no_swap(L), q_swap(R) )
             DLF left_subtree  { one * q_no_swap * Subtree_Tag::LEFT};
             DLF right_subtree { one * q_swap    * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(root, working_qubit, transition);
        }

        { // q_no_swap --w--> ( q_id_1(L), q_id_1(R) )
             DLF left_subtree  { one * q_id_1 * Subtree_Tag::LEFT};
             DLF right_subtree { one * q_id_1 * Subtree_Tag::RIGHT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_swap, working_qubit, transition);
        }

        { // q_swap --w--> ( q_id_1(R), q_id_1(L) )
             DLF left_subtree  { one * q_id_1 * Subtree_Tag::RIGHT};
             DLF right_subtree { one * q_id_1 * Subtree_Tag::LEFT};
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_swap, working_qubit, transition);
        }

        { // q_id_1 --w--> ( q_leaf(L), q_leaf(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_id_1, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_MAJ1) {
        State root   = 0;
        State q_id   = 1;
        State q_p1   = 2;  // 投影掉符號為「1」的分支
        State q_e1   = 3;  // Extract branch "1"
        State q_leaf = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(root);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // root --w--> ( q_id(L), q_id(R) )
             DLF left_subtree  { one * q_id * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_id * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(root, working_qubit, transition);
        }

        { // q_id --w--> ( q_p1(L) + q_e1(R), q_p1(R) + q_e1(L) )
             DLF left_subtree  { one * q_p1 * Subtree_Tag::LEFT, one * q_e1 * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q_p1 * Subtree_Tag::RIGHT, one * q_e1 * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_id, working_qubit, transition);
        }

        { // q_p1 --w--> ( q_leaf, 0 q_leaf )
             DLF left_subtree  { one  * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { zero * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_p1, working_qubit, transition);
        }

        { // q_e1 --w--> ( 0 q_leaf, q_leaf )
             DLF left_subtree  { zero * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { one  * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_e1, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_MAJ_RESULT) {
        State q_a = 0;
        State q_b = 1;
        State q_c = 2;
        State q_d = 3;
        State q_e = 4;
        State q_f = 5;
        State q_g = 6;
        State q_h = 7;
        State q_i = 8;
        State q_j = 9;

        u64 number_of_states = 10;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_a);
        builder.mark_state_final(q_j);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // q_a --w--> ( q_b(L) + q_c(R), q_e(L) + q_d(R) )
             DLF left_subtree  { one * q_b * Subtree_Tag::LEFT, one * q_c * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q_e * Subtree_Tag::LEFT, one * q_d * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_a, working_qubit, transition);
        }

        { // q_b --w--> ( q_f(L), q_f(R) )
             DLF left_subtree  { one * q_f * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_f * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_b, working_qubit, transition);
        }

        { // q_c --w--> ( q_h(R), q_h(L) )
             DLF left_subtree  { one * q_h * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q_h * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_c, working_qubit, transition);
        }

        { // q_d --w--> ( q_f(L), q_i(R) )
             DLF left_subtree  { one * q_f * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_i * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_d, working_qubit, transition);
        }

        { // q_e --w--> ( q_h(R), q_g(L) )
             DLF left_subtree  { one * q_h * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q_g * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_e, working_qubit, transition);
        }

        { // q_f --w--> ( q_j(L), 0 q_j(R) )
             DLF left_subtree  {  one * q_j * Subtree_Tag::LEFT };
             DLF right_subtree { zero * q_j * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_f, working_qubit, transition);
        }

        { // q_g --w--> ( q_j(L), 0 q_j(R) )
             DLF left_subtree  {  one * q_j * Subtree_Tag::RIGHT };
             DLF right_subtree { zero * q_j * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_g, working_qubit, transition);
        }

        { // q_h --w--> ( 0 q_j(L), q_j(R) )
             DLF left_subtree  { zero * q_j * Subtree_Tag::LEFT };
             DLF right_subtree { one  * q_j * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_h, working_qubit, transition);
        }

        { // q_i --w--> ( 0 q_j(R), q_j(L) )
             DLF left_subtree  { zero * q_j * Subtree_Tag::RIGHT };
             DLF right_subtree { one  * q_j * Subtree_Tag::LEFT  };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_i, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);

        do_on_debug({
             result.debug_data = new WTT::Debug_Data();
             result.debug_data->state_names[q_a] = "A";
             result.debug_data->state_names[q_b] = "B";
             result.debug_data->state_names[q_c] = "C";
             result.debug_data->state_names[q_d] = "D";
             result.debug_data->state_names[q_e] = "E";
             result.debug_data->state_names[q_f] = "F";
             result.debug_data->state_names[q_g] = "G";
             result.debug_data->state_names[q_h] = "H";
             result.debug_data->state_names[q_i] = "I";
             result.debug_data->state_names[q_j] = "J";
        });

        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_MAJ_RESULT_12) {
        State q0_a0  = 0;
        State q1_p0  = 1;
        State q1_p1  = 2;
        State p0_p0  = 3;
        State p1_p1  = 4;
        State q_leaf = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q0_a0);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // q0_a0 --w--> ( q1_p0(L) + q1_p1(R), q1_p0(R) + q1_p1(L) )
             DLF left_subtree  { one * q1_p0 * Subtree_Tag::LEFT, one * q1_p1 * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q1_p0 * Subtree_Tag::RIGHT, one * q1_p1 * Subtree_Tag::LEFT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q0_a0, working_qubit, transition);
        }

        { // q1_p0 --w--> ( p0_p0(L), p0_p0(R) )
             DLF left_subtree  { one * p0_p0 * Subtree_Tag::LEFT  };
             DLF right_subtree { one * p0_p0 * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q1_p0, working_qubit, transition);
        }

        { // q1_p1 --w--> ( p0_p0(L), p0_p0(R) )
             DLF left_subtree  { one * p1_p1 * Subtree_Tag::RIGHT };
             DLF right_subtree { one * p1_p1 * Subtree_Tag::LEFT  };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q1_p1, working_qubit, transition);
        }

        { // p0_p0 --w--> ( q_leaf(L), 0 q_leaf(R) )
             DLF left_subtree  {  one * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { zero * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(p0_p0, working_qubit, transition);
        }

        { // p1_p1 --w--> ( q_leaf(L), 0 q_leaf(R) )
             DLF left_subtree  { zero * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { one  * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(p1_p1, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ADDER_MIDDLE) {
        State q_wait     = 0;
        State q_ac       = 1;
        State q_c_swap   = 2;
        State q_c_noswap = 3;
        State q_leaf     = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0;
        Internal_Symbol alt_working_qubit = 1;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_wait);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());

        { // q_wait --w--> ( q_wait(L), q_wait(R) )
             DLF left_subtree  { one * q_wait * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_wait * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_wait, working_qubit, transition);
        }

        { // q_wait --W'--> ( q_idle1(L), q_idle1(R) )
             DLF left_subtree  { one * q_ac * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_ac * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_wait, alt_working_qubit, transition);
        }

        { // q_ac --w--> ( q_noswap(L), q_swap(R) )
             DLF left_subtree  { one * q_c_noswap * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_c_swap   * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_ac, working_qubit, transition);
        }

        { // q_no_swap --w--> ( q_leaf(L), q_leaf(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_c_noswap, working_qubit, transition);
        }

        { // q_swap --w--> ( q_leaf(R), q_leaf(L) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q_leaf * Subtree_Tag::LEFT  };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_c_swap, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);

        if (DEBUG) {
            result.debug_data = new WTT::Debug_Data;
            auto& state_names = result.debug_data->state_names;
            state_names[q_wait]     = "Wait";
            state_names[q_ac]       = "Decide(A)";
            state_names[q_c_swap]   = "Swap";
            state_names[q_c_noswap] = "nowap";
            state_names[q_leaf]     = "leaf";
        }

        return result;
    }

    if (name == Predefined_WTT_Names::TEST_FIXED_ID1) {
        State q_root = 0;
        State q_leaf = 1;

        u64 number_of_states = 2;

        Internal_Symbol working_qubit      = 0;
        Internal_Symbol working_qubit_stop = 1;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_root);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());

        { // q_wait --w--> ( q_leaf(L), q_leaf(R) )
             DLF left_subtree  { one * q_leaf * Subtree_Tag::LEFT  };
             DLF right_subtree { one * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_root, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);
        return result;
    }

    if (name == Predefined_WTT_Names::ECC_BOX1) {
        State q_read_decision = 0;
        State q_no_swap_0     = 1;
        State q_no_swap_1     = 2;
        State q_no_swap_2     = 3;
        State q_swap_0        = 4;
        State q_swap_1        = 5;
        State q_swap_2        = 6;
        State q_leaf          = 7;

        u64 number_of_states = 8;

        Internal_Symbol working_qubit = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_read_decision);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());

        auto add_transition = [&builder, &one, &working_qubit](State source, State left_target, Subtree_Tag left_tag, State right_target, Subtree_Tag right_tag) {
             DLF left_subtree  { one*left_target*left_tag };
             DLF right_subtree { one*right_target*right_tag };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(source, working_qubit, transition);
        };

        add_transition(q_read_decision, q_no_swap_0, Subtree_Tag::LEFT, q_swap_0, Subtree_Tag::RIGHT);

        add_transition(q_no_swap_0, q_no_swap_1, Subtree_Tag::LEFT, q_no_swap_1, Subtree_Tag::RIGHT);
        add_transition(q_no_swap_1, q_no_swap_2, Subtree_Tag::LEFT, q_no_swap_2, Subtree_Tag::RIGHT);
        add_transition(q_no_swap_2, q_leaf, Subtree_Tag::LEFT, q_leaf, Subtree_Tag::RIGHT);

        add_transition(q_swap_0, q_swap_1, Subtree_Tag::LEFT, q_swap_1, Subtree_Tag::RIGHT);
        add_transition(q_swap_1, q_swap_2, Subtree_Tag::LEFT, q_swap_2, Subtree_Tag::RIGHT);
        add_transition(q_swap_2, q_leaf, Subtree_Tag::RIGHT, q_leaf, Subtree_Tag::LEFT);

        WTT result = builder.build(number_of_states);

        if (DEBUG) {
            result.debug_data = new WTT::Debug_Data;
            auto& state_names = result.debug_data->state_names;

            state_names[q_read_decision] = "decide";
            state_names[q_no_swap_0]     = "no_swp0";
            state_names[q_no_swap_1]     = "no_swp1";
            state_names[q_no_swap_2]     = "no_swp2";
            state_names[q_swap_0]        = "swp0";
            state_names[q_swap_1]        = "swp1";
            state_names[q_swap_2]        = "swp2";
            state_names[q_leaf]          = "leaf";
        }

        return result;
     }

     if (name == Predefined_WTT_Names::ECC_BOX2) {
        State q_idle0   = 0;
        State q_idle1   = 1;
        State q_decide  = 2;
        State q_no_swap = 3;
        State q_swap    = 4;
        State q_leaf    = 5;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_idle0);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());

        auto add_transition = [&builder, &one, &working_qubit](State source, State left_target, Subtree_Tag left_tag, State right_target, Subtree_Tag right_tag) {
             DLF left_subtree  { one*left_target*left_tag };
             DLF right_subtree { one*right_target*right_tag };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(source, working_qubit, transition);
        };

        add_transition(q_idle0, q_idle1, Subtree_Tag::LEFT, q_idle1, Subtree_Tag::RIGHT);

        add_transition(q_idle1, q_decide, Subtree_Tag::LEFT, q_decide, Subtree_Tag::RIGHT);
        add_transition(q_decide, q_no_swap, Subtree_Tag::LEFT, q_swap, Subtree_Tag::RIGHT);

        add_transition(q_swap, q_leaf, Subtree_Tag::RIGHT, q_leaf, Subtree_Tag::LEFT);
        add_transition(q_no_swap, q_leaf, Subtree_Tag::LEFT, q_leaf, Subtree_Tag::RIGHT);

        WTT result = builder.build(number_of_states);

        if (DEBUG) {
            result.debug_data = new WTT::Debug_Data;
            auto& state_names = result.debug_data->state_names;

            state_names[q_idle0]   = "idle0";
            state_names[q_idle1]   = "idle1";
            state_names[q_decide]  = "decide";
            state_names[q_no_swap] = "no_swp";
            state_names[q_swap]    = "swp";
            state_names[q_leaf]    = "leaf";
        }

        return result;

    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_BC_CNOT) {
        State q_decide  = 0;
        State q_swap    = 1;
        State q_no_swap = 2;
        State q_leaf    = 3;

        u64 number_of_states = 4;

        Internal_Symbol working_qubit = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_decide);
        builder.mark_state_final(q_leaf);

        Def_Coef one  (ACN::ONE());

        auto add_transition = [&builder, &one, &working_qubit](State source, State left_target, Subtree_Tag left_tag, State right_target, Subtree_Tag right_tag) {
             DLF left_subtree  { one*left_target*left_tag };
             DLF right_subtree { one*right_target*right_tag };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(source, working_qubit, transition);
        };

        add_transition(q_decide, q_no_swap, Subtree_Tag::LEFT, q_swap, Subtree_Tag::RIGHT);
        add_transition(q_swap, q_leaf, Subtree_Tag::RIGHT, q_leaf, Subtree_Tag::LEFT);
        add_transition(q_no_swap, q_leaf, Subtree_Tag::LEFT, q_leaf, Subtree_Tag::RIGHT);

        WTT result = builder.build(number_of_states);

        if (DEBUG) {
            result.debug_data = new WTT::Debug_Data;
            auto& state_names = result.debug_data->state_names;

            state_names[q_decide]  = "decide";
            state_names[q_no_swap] = "no_swp";
            state_names[q_swap]    = "swp";
            state_names[q_leaf]    = "leaf";
        }

        return result;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_BC_RZ) {
        State q_idle    = 0;
        State q_apply   = 1;
        State q_leaf    = 2;

        u64 number_of_states = 3;

        Internal_Symbol working_qubit = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_idle);
        builder.mark_state_final(q_leaf);

        Def_Coef one          (ACN::ONE());
        Def_Coef omega1       (ACN(1, 0, 0, 0, 0));
        Def_Coef minus_omega3 (ACN(0, 0, 0, -1, 0));

        auto add_transition = [&builder, &one, &working_qubit](State source, const Def_Linear_Form& left_successor, const Def_Linear_Form& right_successor) {
             DLF left_subtree  { left_successor };
             DLF right_subtree { right_successor };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(source, working_qubit, transition);
        };

        add_transition(q_idle,  one*q_apply*Subtree_Tag::LEFT, one*q_apply*Subtree_Tag::RIGHT);
        add_transition(q_apply, minus_omega3*q_leaf*Subtree_Tag::LEFT, omega1*q_leaf*Subtree_Tag::RIGHT);

        WTT result = builder.build(number_of_states);

        if (DEBUG) {
            result.debug_data = new WTT::Debug_Data;
            auto& state_names = result.debug_data->state_names;

            state_names[q_idle]  = "idle";
            state_names[q_apply] = "apply";
            state_names[q_leaf]  = "leaf";
        }

        return result;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_BC_X) {
        State q_idle    = 0;
        State q_apply   = 1;
        State q_leaf    = 2;

        u64 number_of_states = 3;

        Internal_Symbol working_qubit = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_idle);
        builder.mark_state_final(q_leaf);

        Def_Coef one          (ACN::ONE());

        auto add_transition = [&builder, &one, &working_qubit](State source, const Def_Linear_Form& left_successor, const Def_Linear_Form& right_successor) {
             DLF left_subtree  { left_successor };
             DLF right_subtree { right_successor };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(source, working_qubit, transition);
        };

        add_transition(q_idle,  one*q_apply*Subtree_Tag::LEFT, one*q_apply*Subtree_Tag::RIGHT);
        add_transition(q_apply, one*q_leaf*Subtree_Tag::RIGHT, one*q_leaf*Subtree_Tag::LEFT);

        WTT result = builder.build(number_of_states);

        if (DEBUG) {
            result.debug_data = new WTT::Debug_Data;
            auto& state_names = result.debug_data->state_names;

            state_names[q_idle]  = "idle";
            state_names[q_apply] = "apply-swap";
            state_names[q_leaf]  = "leaf";
        }

        return result;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_BC_H) {
        State q_apply1  = 0;
        State q_apply2  = 1;
        State q_leaf    = 2;

        u64 number_of_states = 3;

        Internal_Symbol working_qubit = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_apply1);
        builder.mark_state_final(q_leaf);

        Def_Coef one         (ACN::ONE());
        Def_Coef invsqrt     (ACN::ONE_OVER_SQRT2());
        Def_Coef neg_invsqrt (-ACN::ONE_OVER_SQRT2());

        { // q_apply1
             DLF left_subtree  { invsqrt * q_apply2 * Subtree_Tag::LEFT,     invsqrt * q_apply2 * Subtree_Tag::RIGHT };
             DLF right_subtree { invsqrt * q_apply2 * Subtree_Tag::LEFT, neg_invsqrt * q_apply2 * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_apply1, working_qubit, transition);
        }

        { // q_apply2
             DLF left_subtree  { invsqrt * q_leaf * Subtree_Tag::LEFT,     invsqrt * q_leaf * Subtree_Tag::RIGHT };
             DLF right_subtree { invsqrt * q_leaf * Subtree_Tag::LEFT, neg_invsqrt * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_apply2, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);

        if (DEBUG) {
            result.debug_data = new WTT::Debug_Data;
            auto& state_names = result.debug_data->state_names;

            state_names[q_apply1] = "apply1";
            state_names[q_apply2] = "apply2";
            state_names[q_leaf]   = "leaf";
        }

        return result;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_BC_S) {
        State q_apply1  = 0;
        State q_apply2  = 1;
        State q_leaf    = 2;

        u64 number_of_states = 3;

        Internal_Symbol working_qubit = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_apply1);
        builder.mark_state_final(q_leaf);

        Def_Coef one          (ACN::ONE());
        Def_Coef complex_unit (ACN(0, 0, 1, 0, 0));

        { // q_apply1
             DLF left_subtree  {          one * q_apply2 * Subtree_Tag::LEFT };
             DLF right_subtree { complex_unit * q_apply2 * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_apply1, working_qubit, transition);
        }

        { // q_apply2
             DLF left_subtree  {          one * q_leaf * Subtree_Tag::LEFT };
             DLF right_subtree { complex_unit * q_leaf * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_apply2, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);

        if (DEBUG) {
            result.debug_data = new WTT::Debug_Data;
            auto& state_names = result.debug_data->state_names;

            state_names[q_apply1] = "apply1";
            state_names[q_apply2] = "apply2";
            state_names[q_leaf]   = "leaf";
        }

        return result;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_H_STAGE) {
        ACN     invsqrt =  ACN::ONE_OVER_SQRT2();
        ACN neg_invsqrt = -ACN::ONE_OVER_SQRT2();

        WTT result = construct_stage_for_hamiltonian(invsqrt, invsqrt, invsqrt, neg_invsqrt);

        return result;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_SQRT_X_STAGE) {
        ACN ll (1, 0,  1, 0, 2);
        ACN lr (1, 0, -1, 0, 2);
        ACN rl (1, 0, -1, 0, 2);
        ACN rr (1, 0,  1, 0, 2);

        WTT result = construct_stage_for_hamiltonian(ll, lr, rl, rr);

        return result;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_S_STAGE) {
        ACN ll (1, 0,  0, 0, 0);
        ACN lr = ACN::ZERO();
        ACN rl = ACN::ZERO();;
        ACN rr (0, 0,  1, 0, 0); // Complex unit

        WTT result = construct_stage_for_hamiltonian(ll, lr, rl, rr);

        return result;
    }
    if (name == Predefined_WTT_Names::HAMILTONIAN_RZZ) {
        WTT cnot = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_BC_CNOT, swta_metadata);
        WTT rz   = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_BC_RZ, swta_metadata);

        WTT rzz_box12 = compose_wtts_sequentially(cnot, rz);
        WTT rzz_box   = compose_wtts_sequentially(rzz_box12, cnot);

        return rzz_box;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_RXX) {
        WTT hadamard = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_BC_H, swta_metadata);
        WTT rzz_box  = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_RZZ, swta_metadata);

        WTT rxx_box12 = compose_wtts_sequentially(hadamard, rzz_box);
        WTT rxx_box   = compose_wtts_sequentially(rxx_box12, hadamard);
        return rxx_box;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_RYY) {
        WTT s        = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_BC_S, swta_metadata);
        WTT rxx_box  = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_RXX, swta_metadata);

        WTT ryy_box12 = compose_wtts_sequentially(s, rxx_box);
        WTT ryy_box   = compose_wtts_sequentially(ryy_box12, s);
        return ryy_box;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_UZZ) {
        WTT rzz_box  = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_RZZ, swta_metadata);
        WTT x_component = get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_BC_X, swta_metadata);
        WTT uzz_box = compose_wtts_sequentially(rzz_box, x_component);
        return uzz_box;
    }

    if (name == Predefined_WTT_Names::HAMILTONIAN_LAST_X_STAGE) {
        State q_wait  = 0;
        State q_apply = 1;
        State q_leaf  = 2;

        u64 number_of_states = 3;

        Internal_Symbol working_qubit      = 0;
        Internal_Symbol working_qubit_stop = 1;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        WTT_Builder builder (metadata);

        builder.mark_state_initial(q_wait);
        builder.mark_state_final(q_leaf);

        Def_Coef one          (ACN::ONE());

        {
             DLF left_subtree  { one * q_wait * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_wait * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_wait, working_qubit, transition);
        }

        {
             DLF left_subtree  { one * q_apply * Subtree_Tag::LEFT };
             DLF right_subtree { one * q_apply * Subtree_Tag::RIGHT };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_wait, working_qubit_stop, transition);
        }

        {
             DLF left_subtree  { one * q_leaf * Subtree_Tag::RIGHT };
             DLF right_subtree { one * q_leaf * Subtree_Tag::LEFT  };
             auto transition = synthetize_wtt_transition(left_subtree, right_subtree);
             builder.add_transition(q_apply, working_qubit, transition);
        }

        WTT result = builder.build(number_of_states);

        if (DEBUG) {
            result.debug_data = new WTT::Debug_Data;
            auto& state_names = result.debug_data->state_names;

            state_names[q_wait]  = "wait";
            state_names[q_apply] = "apply-swap";
            state_names[q_leaf]  = "leaf";
        }

        return result;
    }

    throw std::runtime_error("Unknown WTT. " + std::to_string(static_cast<u64>(name)));
}


// 依照預設的名稱產生 SWTA 範例：根據指定名稱決定狀態數與內部符號，
// 為每個狀態填入以線性形式表示的轉移，並標記初始與接受葉節點，
// 方便測試或示範時直接取得可運行的樹自動機配置。
SWTA get_predefined_swta(Predefined_SWTA_Names name) {
    using DLF = std::vector<Def_Linear_Form>;
    using ACN = Algebraic_Complex_Number;

    if (name == Predefined_SWTA_Names::BV_EXAMPLE_10STAR_PRE) {
        SWTA::Metadata metadata {.number_of_internal_symbols = 2, .number_of_colors = 1};
        SWTA::Transition_Builder builder (metadata);

        State q0    = 0;
        State q1    = 1;
        State q_bot = 2;  // 接受任何高度的空樹

        Internal_Symbol sym_w = 0;
        Internal_Symbol sym_a = 1;

        Color c = 0;

        { // 轉移 q0 ----> w(q0, q_bot)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ONE()) * q0};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q0, sym_w, c, transition);
        }

        { // 轉移 q0 ----> a(q_bot, q1)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * q1};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q0, sym_a, c, transition);
        }

        { // 轉移 q_bot ----> a(q_bot, q_bot)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_bot, sym_w, c, transition);
             builder.add_transition(q_bot, sym_a, c, transition);
        }

        Bit_Set leaf_states (3, {q1, q_bot});
        std::vector<State> initial_states {q0};
        auto transition_fn = builder.build(3);

        SWTA result (transition_fn, initial_states, leaf_states);
        return result;
    }

    if (name == Predefined_SWTA_Names::BV_EXAMPLE_10STAR_POST) {
        State q_g   = 0;
        State q_h   = 1;
        State q_c   = 2;
        State q_bot = 3;

        Internal_Symbol sym_w = 0;
        Internal_Symbol sym_a = 1;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        Color c = 0;

        { // g ----> w(0 q_bot, h)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * q_h};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_g, sym_w, c, transition);
        }

        { // g ----> a(0 q_bot, c)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * q_c};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_g, sym_a, c, transition);
        }

        { // h ----> w(g, 0 q_bot)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ONE()) * q_g};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_h, sym_w, c, transition);
        }

        { // h ----> a(0 q_bot, c)
             DLF left_subtree  {Def_Coef(Algebraic_Complex_Number::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(Algebraic_Complex_Number::ONE()) * q_c};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_h, sym_a, c, transition);
        }

        builder.add_bot_state_transitions(q_bot);

        std::vector<State> initial_states ({q_g});
        Bit_Set leaf_states (4, {q_bot, q_c});
        auto transition_fn = builder.build(4);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::BV_EXAMPLE_10STAR_RESULT || name == Predefined_SWTA_Names::TEST_BV_EXAMPLE_AFTER_STEP3) {
        State q_gamma   = 0;
        State q_delta   = 1;
        State q_epsilon = 2;
        State q_mu      = 3;
        State q_sigma   = 4;
        State q_bot     = 5;

        u64 state_cnt = 6;

        Internal_Symbol sym_w = 0;
        Internal_Symbol sym_a = 1;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        Color c = 0;

        ACN one_half (1, 0, 0, 0, 2);

        { // gamma ----> w(1/2 delta + 1/2 epsilon, 1/2 delta - 1/2 epsilon)
             DLF left_subtree  {Def_Coef(one_half) * q_delta, Def_Coef(one_half) * q_epsilon};
             DLF right_subtree {Def_Coef(one_half) * q_delta, Def_Coef(-one_half) * q_epsilon};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_gamma, sym_w, c, transition);
        }

        { // gamma ----> a(0 q_bot, q_sigma)
             DLF left_subtree  {Def_Coef(ACN::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(ACN::ONE()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_gamma, sym_a, c, transition);
        }

        { // delta ----> w(q_gamma, 0 q_bot)
             DLF left_subtree  {Def_Coef(ACN::ONE()) * q_gamma};
             DLF right_subtree {Def_Coef(ACN::ZERO()) * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_delta, sym_w, c, transition);
        }

        { // delta ----> a(0 q_bot, q_sigma)
             DLF left_subtree  {Def_Coef(ACN::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(ACN::ONE()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_delta, sym_a, c, transition);
        }

        { // epsilon ----> w(q_mu, 0 q_bot)
             DLF left_subtree  {Def_Coef(ACN::ONE()) * q_mu};
             DLF right_subtree {Def_Coef(ACN::ZERO()) * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_epsilon, sym_w, c, transition);
        }

        { // epsilon ----> a(0 q_bot, - q_sigma)
             DLF left_subtree  {Def_Coef(ACN::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(-ACN::ONE()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_epsilon, sym_a, c, transition);
        }

        { // mu ----> w(1/2 epsilon + 1/2 delta, 1/2 delta - 1/2 epsilon)
             DLF left_subtree  {Def_Coef(one_half) * q_epsilon, Def_Coef(one_half) * q_delta};
             DLF right_subtree {Def_Coef(one_half) * q_epsilon, Def_Coef(-one_half) * q_delta};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_mu, sym_w, c, transition);
        }

        { // mu ----> a(0 q_bot, q_sigma)
             DLF left_subtree  {Def_Coef(ACN::ZERO()) * q_bot};
             DLF right_subtree {Def_Coef(-ACN::ONE()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_mu, sym_a, c, transition);
        }

        builder.add_bot_state_transitions(q_bot);

        std::vector<State> initial_states ({q_gamma});
        Bit_Set leaf_states (state_cnt, {q_bot, q_sigma});
        auto transition_fn = builder.build(state_cnt);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::TRIVIAL_BOT) {
        State q_bot   = 0;

        Internal_Symbol sym_w = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        Color c = 0;

        builder.add_bot_state_transitions(q_bot);

        std::vector<State> initial_states ({q_bot});
        Bit_Set leaf_states (1, {q_bot});
        auto transition_fn = builder.build(1);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::TRIVIAL_ONES) {
        State q_one   = 0;

        Internal_Symbol sym_w = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        Color c = 0;

        { // q_one ----> w(q_one, q_one)
             DLF left_subtree  {Def_Coef(ACN::ONE()) * q_one};
             DLF right_subtree {Def_Coef(ACN::ONE()) * q_one};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_one, sym_w, c, transition);
        }

        std::vector<State> initial_states ({q_one});
        Bit_Set leaf_states (1, {q_one});
        auto transition_fn = builder.build(1);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::TEST_BV_EXAMPLE_AFTER_STEP1) {
        State q_alpha  = 0;
        State q_beta   = 1;

        Internal_Symbol working_qubit = 0, ancilla = 1;
        Color c = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        { // q_alpha --w--> (1/sqrt(2) q_alpha, 1/sqrt(2) q_alpha)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_alpha};
             DLF right_subtree {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_alpha};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_alpha, working_qubit, c, transition);
        }

        { // q_alpha --a--> (1/sqrt(2) q_beta, 1/sqrt(2) q_beta)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_beta};
             DLF right_subtree {Def_Coef(-ACN::ONE_OVER_SQRT2()) * q_beta};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_alpha, ancilla, c, transition);
        }

        std::vector<State> initial_states ({q_alpha});
        Bit_Set leaf_states (2, {q_beta});
        auto transition_fn = builder.build(2);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::TEST_BV_EXAMPLE_AFTER_STEP2) {
        State q_gamma   = 0;
        State q_delta   = 1;
        State q_epsilon = 2;
        State q_mu      = 3;
        State q_sigma   = 4;

        u64 number_of_states = 5;

        Internal_Symbol working_qubit = 0, ancilla = 1;
        Color c = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        { // q_gamma --w--> (1/sqrt(2) q_delta, 1/sqrt(2) q_epsilon)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_delta};
             DLF right_subtree {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_epsilon};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_gamma, working_qubit, c, transition);
        }

        { // q_delta --w--> (1/sqrt(2) q_beta, 1/sqrt(2) q_beta)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_gamma};
             DLF right_subtree {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_gamma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_delta, working_qubit, c, transition);
        }

        { // q_epsilon --w--> (1/sqrt(2) q_mu, 1/sqrt(2) q_mu)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_mu};
             DLF right_subtree {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_mu};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_epsilon, working_qubit, c, transition);
        }

        { // q_mu --w--> (1/sqrt(2) q_epsilon, 1/sqrt(2) q_delta)
             DLF left_subtree  {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_epsilon};
             DLF right_subtree {Def_Coef(ACN::ONE_OVER_SQRT2()) * q_delta};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_mu, working_qubit, c, transition);
        }

        { // q_gamma --a--> (1/sqrt(2) q_sigma, -1/sqrt(2) q_sigma)
             DLF left_subtree  {Def_Coef( ACN::ONE_OVER_SQRT2()) * q_sigma};
             DLF right_subtree {Def_Coef(-ACN::ONE_OVER_SQRT2()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_gamma, ancilla, c, transition);
        }

        { // q_delta --a--> (1/sqrt(2) q_sigma, -1/sqrt(2) q_sigma)
             DLF left_subtree  {Def_Coef( ACN::ONE_OVER_SQRT2()) * q_sigma};
             DLF right_subtree {Def_Coef(-ACN::ONE_OVER_SQRT2()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_delta, ancilla, c, transition);
        }

        { // q_epsilon --a--> (-1/sqrt(2) q_sigma, 1/sqrt(2) q_sigma)
             DLF left_subtree  {Def_Coef(-ACN::ONE_OVER_SQRT2()) * q_sigma};
             DLF right_subtree {Def_Coef( ACN::ONE_OVER_SQRT2()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_epsilon, ancilla, c, transition);
        }

        { // q_mu --a--> (-1/sqrt(2) q_sigma, 1/sqrt(2) q_sigma)
             DLF left_subtree  {Def_Coef(-ACN::ONE_OVER_SQRT2()) * q_sigma};
             DLF right_subtree {Def_Coef( ACN::ONE_OVER_SQRT2()) * q_sigma};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_mu, ancilla, c, transition);
        }

        std::vector<State> initial_states ({q_gamma});
        Bit_Set leaf_states (number_of_states, {q_sigma});
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::GROVER_ALL_BASIS) {
        State q_no_one_seen = 0;
        State q_all_zeros   = 1;
        State q_bot         = 2;
        State q_anc         = 3;
        State q_last_anc    = 4;
        State q_leaf        = 5;

        Internal_Symbol sym_w = 0;

        u64 number_of_states = 6;

        Internal_Symbol working_qubit = 0, ancilla = 1, last_working_qubit = 2;
        Color color = 0;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 3, .number_of_colors = 1 };
        SWTA::Transition_Builder builder (metadata);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // q_no_one_seen --w--> (q_no_one_seen, q_all_zeros)
             DLF left_subtree  {one * q_no_one_seen};
             DLF right_subtree {one * q_all_zeros};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_one_seen, working_qubit, color, transition);
        }

        { // q_no_one_seen --W--> (q_last_anc, 0 q_bot)
             DLF left_subtree  {one  * q_last_anc};
             DLF right_subtree {zero * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_no_one_seen, last_working_qubit, color, transition);
        }

        { // q_all_zeros --w--> (q_anc, 0 q_bot)
             DLF left_subtree  {one  * q_anc};
             DLF right_subtree {zero * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_all_zeros, working_qubit, color, transition);
        }

        { // q_anc --a--> (q_all_zeros, q_all_zeros)
             DLF left_subtree  {one * q_all_zeros};
             DLF right_subtree {one * q_all_zeros};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_anc, ancilla, color, transition);
        }

        { // q_all_zeros --W'--> (q_last_anc, 0 q_bot)
             DLF left_subtree  {one  * q_last_anc};
             DLF right_subtree {zero * q_bot};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_all_zeros, last_working_qubit, color, transition);
        }

        { // q_last_anc --a--> (q_leaf, q_leaf)
             DLF left_subtree  {one * q_leaf};
             DLF right_subtree {one * q_leaf};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q_last_anc, ancilla, color, transition);
        }

        builder.add_bot_state_transitions(q_bot);

        std::vector<State> initial_states ({q_no_one_seen});
        Bit_Set leaf_states (number_of_states, {q_bot, q_leaf});
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::TEST_ADDER_ALL_3BASIS) {
        State q3 = 0;
        State q2 = 1;
        State q1 = 2;
        State q0 = 3;

        u64 number_of_states = 4;

        Internal_Symbol working_qubit = 0;
        Color one_left  = 0;
        Color one_right = 1;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 2 };
        SWTA::Transition_Builder builder (metadata);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        { // q3 --one_left--> (q2, 0q2)
             DLF left_subtree  {one  * q2};
             DLF right_subtree {zero * q2};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q3, working_qubit, one_left, transition);
        }

        { // q3 --one_right--> (0q2, q2)
             DLF left_subtree  {zero * q2};
             DLF right_subtree {one  * q2};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q3, working_qubit, one_right, transition);
        }

        { // q2 --one_left--> (q1, 0q1)
             DLF left_subtree  {one  * q1};
             DLF right_subtree {zero * q1};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q2, working_qubit, one_left, transition);
        }

        { // q2 --one_right--> (0q1, q1)
             DLF left_subtree  {zero * q1};
             DLF right_subtree {one  * q1};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q2, working_qubit, one_right, transition);
        }

        { // q1 --one_left--> (q0, 0 q0)
             DLF left_subtree  {one  * q0};
             DLF right_subtree {zero * q0};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q1, working_qubit, one_left, transition);
        }

        { // q1 --one_right--> (0 q0, q0)
             DLF left_subtree  {zero * q0};
             DLF right_subtree {one  * q0};
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(q1, working_qubit, one_right, transition);
        }

        std::vector<State> initial_states ({q3});
        Bit_Set leaf_states (number_of_states, { q0 });
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        return result;
    }

    if (name == Predefined_SWTA_Names::ADDER_POST) {
        State aQ_bQ_cQ = 0;

        State aQ_bQ_c0 = 1;
        State aQ_bQ_c1 = 2;

        State aQ_b0_c0_parity0 = 3;
        State aQ_b0_c0_parity1 = 4;

        State aQ_b0_c1_parity0 = 5;
        State aQ_b0_c1_parity1 = 6;

        State aQ_b1_c0_parity0 = 7;
        State aQ_b1_c0_parity1 = 8;

        State aQ_b1_c1_parity0 = 9;
        State aQ_b1_c1_parity1 = 10;

        State aQ_b0_c0_parity0_fin = 11;
        State aQ_b0_c0_parity1_fin = 12;

        State aQ_b0_c1_parity0_fin = 13;
        State aQ_b0_c1_parity1_fin = 14;

        State aQ_b1_c0_parity0_fin = 15;
        State aQ_b1_c0_parity1_fin = 16;

        State aQ_b1_c1_parity0_fin = 17;
        State aQ_b1_c1_parity1_fin = 18;

        State c0_last = 19;
        State c1_last = 20;

        State ida  = 21;
        State idb  = 22;
        State id_last_b  = 23;
        State id_last_a  = 24;
        State id_last_c  = 25;
        State leaf = 26;

        u64 number_of_states = 27;

        Internal_Symbol working_qubit = 0;
        Internal_Symbol working_qubit_stop = 1;

        Color one_left  = 0;
        Color one_right = 1;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 2 };
        SWTA::Transition_Builder builder (metadata);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        DLF zero_ida_branch     = { zero * ida };
        DLF zero_idb_branch     = { zero * idb };
        DLF zero_ida_fin_branch = { zero * id_last_a };
        DLF zero_idb_fin_branch = { zero * id_last_b };
        DLF zero_leaf_branch    = { zero * leaf };
        DLF one_leaf_branch     = { one * leaf };

        Alternative_Id_Helper id_helper = {
            .ida = ida,
            .idb = idb,
            .ida_branch = zero_ida_branch,
            .idb_branch = zero_idb_branch,
            .ida_fin_branch = zero_ida_fin_branch,
            .idb_fin_branch = zero_idb_fin_branch,
        };

        auto add_transition = [&builder, &working_qubit](State source, Color color, const Def_Linear_Form& left_successor, const Def_Linear_Form& right_successor) {
            DLF left_subtree   {left_successor};
            DLF right_subtree  {right_successor};
            auto transition = synthetize_swta_transition(left_subtree, right_subtree);
            builder.add_transition(source, working_qubit, color, transition);
        };

        auto add_fin_transition = [&builder, &working_qubit_stop](State source, Color color, const Def_Linear_Form& left_successor, const Def_Linear_Form& right_successor) {
            DLF left_subtree   {left_successor};
            DLF right_subtree  {right_successor};
            auto transition = synthetize_swta_transition(left_subtree, right_subtree);
            builder.add_transition(source, working_qubit_stop, color, transition);
        };

        // ------------------------- aQ_bQ_cQ -------------------------
        add_transition(aQ_bQ_cQ, one_left,  one*aQ_bQ_c0, zero*idb);
        add_transition(aQ_bQ_cQ, one_right, zero*idb, one*aQ_bQ_c1);

        // ------------------------- aQ_bQ_c0 -------------------------
        add_transition(aQ_bQ_c0, one_left,  one*aQ_b0_c0_parity0, one*aQ_b0_c0_parity1);
        add_transition(aQ_bQ_c0, one_right, one*aQ_b1_c0_parity0, one*aQ_b1_c0_parity1);

        add_fin_transition(aQ_bQ_c0, one_left,  one*aQ_b0_c0_parity0_fin, one*aQ_b0_c0_parity1_fin);
        add_fin_transition(aQ_bQ_c0, one_right, one*aQ_b1_c0_parity0_fin, one*aQ_b1_c0_parity1_fin);

        // ------------------------- aQ_bQ_c1 -------------------------
        add_transition(aQ_bQ_c1, one_left,  one*aQ_b0_c1_parity0, one*aQ_b0_c1_parity1);
        add_transition(aQ_bQ_c1, one_right, one*aQ_b1_c1_parity0, one*aQ_b1_c1_parity1);

        add_fin_transition( aQ_bQ_c1, one_left, one*aQ_b0_c1_parity0_fin, one*aQ_b0_c1_parity1_fin);
        add_fin_transition(aQ_bQ_c1, one_right, one*aQ_b1_c1_parity0_fin, one*aQ_b1_c1_parity1_fin);

        auto ets = idb; // Empty tree successor

        // ------------------------- aQ_b0_c0_parity0 -------------------------
        add_transition(aQ_b0_c0_parity0, one_left,  one*aQ_bQ_c0, zero*ets); // Parity(0, 0, 0) == 0
        add_transition(aQ_b0_c0_parity0, one_right, zero*ets,     zero*ets); // Parity(1, 0, 0) != 0

        // ------------------------- aQ_b0_c0_parity1 -------------------------
        add_transition(aQ_b0_c0_parity1, one_left,  zero*ets, zero*ets);     // Parity(0, 0, 0) != 1
        add_transition(aQ_b0_c0_parity1, one_right, zero*ets, one*aQ_bQ_c0); // Parity(1, 0, 0) == 1

        // ------------------------- aQ_b0_c1_parity0 -------------------------
        add_transition(aQ_b0_c1_parity0, one_left,  zero*ets, zero*ets);     // Parity(0, 0, 1) != 0
        add_transition(aQ_b0_c1_parity0, one_right, zero*ets, one*aQ_bQ_c1); // Parity(1, 0, 1) == 0

        // ------------------------- aQ_b0_c1_parity1 -------------------------
        add_transition(aQ_b0_c1_parity1, one_left,  one*aQ_bQ_c0, zero*ets);  // Parity(0, 0, 1) == 1
        add_transition(aQ_b0_c1_parity1, one_right, zero*ets, zero*ets);      // Parity(1, 0, 1) != 1

        // ------------------------- aQ_b1_c0_parity0 -------------------------
        add_transition(aQ_b1_c0_parity0, one_left,  zero*ets, zero*ets);      // Parity(0, 1, 0) != 0
        add_transition(aQ_b1_c0_parity0, one_right, zero*ets, one*aQ_bQ_c1);  // Parity(1, 1, 0) == 0

        // ------------------------- aQ_b1_c0_parity1 -------------------------
        add_transition(aQ_b1_c0_parity1, one_left,  one*aQ_bQ_c0, zero*ets);  // Parity(0, 1, 0) == 1
        add_transition(aQ_b1_c0_parity1, one_right, zero*ets, zero*ets);      // Parity(1, 1, 0) != 1

        // ------------------------- aQ_b1_c0_parity1 -------------------------
        add_transition(aQ_b1_c0_parity1, one_left,  one*aQ_bQ_c0, zero*ets);  // Parity(0, 1, 0) == 1
        add_transition(aQ_b1_c0_parity1, one_right, zero*ets, zero*ets);      // Parity(1, 1, 0) != 1

        // ------------------------- aQ_b1_c1_parity0 -------------------------
        add_transition(aQ_b1_c1_parity0, one_left,  one*aQ_bQ_c1, zero*ets);  // Parity(0, 1, 1) == 0
        add_transition(aQ_b1_c1_parity0, one_right, zero*ets, zero*ets);      // Parity(1, 1, 1) != 0

        // ------------------------- aQ_b1_c1_parity1 -------------------------
        add_transition(aQ_b1_c1_parity1, one_left,  zero*ets, zero*ets);      // Parity(0, 1, 1) != 1
        add_transition(aQ_b1_c1_parity1, one_right, zero*ets, one*aQ_bQ_c1);  // Parity(1, 1, 1) == 1

        //
        //  FIN states
        //

        // ------------------------- aQ_b0_c0_parity0_fin -------------------------
        add_transition(aQ_b0_c0_parity0_fin, one_left,  one*c0_last,    zero*id_last_c); // Parity(0, 0, 0) == 0
        add_transition(aQ_b0_c0_parity0_fin, one_right, zero*id_last_c, zero*id_last_c); // Parity(1, 0, 0) != 0

        // ------------------------- aQ_b0_c0_parity1_fin -------------------------
        add_transition(aQ_b0_c0_parity1_fin, one_left,  zero*id_last_c, zero*id_last_c); // Parity(0, 0, 0) != 1
        add_transition(aQ_b0_c0_parity1_fin, one_right, zero*id_last_c, one*c0_last);    // Parity(1, 0, 0) == 1

        // ------------------------- aQ_b0_c1_parity0_fin -------------------------
        add_transition(aQ_b0_c1_parity0_fin, one_left,  zero*id_last_c, zero*id_last_c); // Parity(0, 0, 1) != 0
        add_transition(aQ_b0_c1_parity0_fin, one_right, zero*id_last_c, one*c1_last);    // Parity(1, 0, 1) == 0

        // ------------------------- aQ_b0_c1_parity1_fin -------------------------
        add_transition(aQ_b0_c1_parity1_fin, one_left,  one*c0_last,    zero*id_last_c); // Parity(0, 0, 1) == 1
        add_transition(aQ_b0_c1_parity1_fin, one_right, zero*id_last_c, zero*id_last_c); // Parity(1, 0, 1) != 1

        // ------------------------- aQ_b1_c0_parity0_fin -------------------------
        add_transition(aQ_b1_c0_parity0_fin, one_left,  zero*id_last_c, zero*id_last_c);  // Parity(0, 1, 0) != 0
        add_transition(aQ_b1_c0_parity0_fin, one_right, zero*id_last_c, one*c1_last);     // Parity(1, 1, 0) == 0

        // ------------------------- aQ_b1_c0_parity1_fin -------------------------
        add_transition(aQ_b1_c0_parity1_fin, one_left,  one*c0_last,    zero*id_last_c);  // Parity(0, 1, 0) == 1
        add_transition(aQ_b1_c0_parity1_fin, one_right, zero*id_last_c, zero*id_last_c);  // Parity(1, 1, 0) != 1

        // ------------------------- aQ_b1_c1_parity0_fin -------------------------
        add_transition(aQ_b1_c1_parity0_fin, one_left,  one*c1_last,    zero*id_last_c);  // Parity(0, 1, 1) == 0
        add_transition(aQ_b1_c1_parity0_fin, one_right, zero*id_last_c, zero*id_last_c);  // Parity(1, 1, 1) != 0

        // ------------------------- aQ_b1_c1_parity1_fin -------------------------
        add_transition(aQ_b1_c1_parity1_fin, one_left,  zero*id_last_c, zero*id_last_c);  // Parity(0, 1, 1) != 1
        add_transition(aQ_b1_c1_parity1_fin, one_right, zero*id_last_c, one*c1_last);     // Parity(1, 1, 1) == 1

        add_transition(c0_last, one_left, one*leaf,  zero*leaf);
        add_transition(c1_last, one_left, zero*leaf, one*leaf);

        // The id state
        add_transition(ida, one_left,  zero*idb, zero*idb);
        add_transition(ida, one_right, zero*idb, zero*idb);

        add_transition(idb, one_left,  zero*ida, zero*ida);
        add_transition(idb, one_right, zero*ida, zero*ida);

        add_fin_transition(idb, one_left,  zero*id_last_a, zero*id_last_a);
        add_fin_transition(idb, one_right, zero*id_last_a, zero*id_last_a);

        // --------------- id_last_b ---------------
        add_transition(id_last_b, one_left,  zero*id_last_a, zero*id_last_a);
        add_transition(id_last_b, one_right, zero*id_last_a, zero*id_last_a);

        add_transition(id_last_a, one_left,  zero*id_last_c, zero*id_last_c);
        add_transition(id_last_a, one_right, zero*id_last_c, zero*id_last_c);

        add_transition(id_last_c, one_left, zero*leaf, zero*leaf);

        std::vector<State> initial_states ({aQ_bQ_cQ});
        Bit_Set leaf_states (number_of_states, { leaf });
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        do_on_debug({
            result.debug_data = new SWTA::Debug_Data();
            auto& state_names = result.debug_data->state_names;

            state_names[aQ_bQ_cQ]     = "(a b c)";

            state_names[aQ_bQ_c0]     = "(a b 0)";
            state_names[aQ_bQ_c1]     = "(a b 1)";

            state_names[aQ_b0_c0_parity0]     = "(a 0 0).p0";
            state_names[aQ_b0_c0_parity1]     = "(a 0 0).p1";

            state_names[aQ_b0_c1_parity0]     = "(a 0 1).p0";
            state_names[aQ_b0_c1_parity1]     = "(a 0 1).p1";

            state_names[aQ_b1_c0_parity0]     = "(a 1 0).p0";
            state_names[aQ_b1_c0_parity1]     = "(a 1 0).p1";

            state_names[aQ_b1_c1_parity0]     = "(a 1 1).p0";
            state_names[aQ_b1_c0_parity1]     = "(a 1 0).p1";

            state_names[aQ_b0_c0_parity0_fin] = "F(a 0 0).p0";
            state_names[aQ_b0_c0_parity1_fin] = "F(a 0 0).p1";

            state_names[aQ_b0_c1_parity0_fin] = "F(a 0 0).p0";
            state_names[aQ_b0_c1_parity1_fin] = "F(a 0 0).p1";

            state_names[aQ_b1_c0_parity0_fin] = "F(a 1 0).p0";
            state_names[aQ_b1_c0_parity1_fin] = "F(a 1 0).p1";

            state_names[aQ_b1_c1_parity0_fin] = "F(a 1 1).p0";
            state_names[aQ_b1_c1_parity1_fin] = "F(a 1 1).p1";

            state_names[c0_last] = "c0_last";
            state_names[c1_last] = "c1_last";
            state_names[ida]      = "IDa";
            state_names[idb]      = "IDb";
            state_names[id_last_b] = "ID_last_b";
            state_names[id_last_a] = "ID_last_a";
            state_names[id_last_c] = "ID_last_c";
            state_names[leaf]    = "leaf";
        });

        return result;
    }

    if (name == Predefined_SWTA_Names::ADDER_PRE) {
        State q_c  = 0;
        State q_a  = 1;
        State q_b  = 2;
        State id   = 3;
        State ac   = 4;
        State c    = 5;
        State leaf = 6;

        u64 number_of_states = 7;

        Internal_Symbol working_qubit      = 0;
        Internal_Symbol working_qubit_stop = 1;

        Color one_left  = 0;
        Color one_right = 1;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 2 };
        SWTA::Transition_Builder builder (metadata);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        DLF zero_id_branch =   { zero * id   };
        DLF zero_leaf_branch = { zero * leaf };
        DLF one_leaf_branch =  { one * leaf  };

        auto add_left_transition = [&one, &builder, &zero_id_branch, &working_qubit, &one_left, &id](State source, State destination) {
             if (destination == id) {
                  auto transition = synthetize_swta_transition(zero_id_branch, zero_id_branch);
                  builder.add_transition(source, working_qubit, one_left, transition);
                  return;
             }
             DLF left_subtree  {one  * destination };
             auto transition = synthetize_swta_transition(left_subtree, zero_id_branch);
             builder.add_transition(source, working_qubit, one_left, transition);

        };

        auto add_right_transition = [&one, &builder, &zero_id_branch, &working_qubit, &one_right, &id](State source, State destination) {
             if (destination == id) {
                  auto transition = synthetize_swta_transition(zero_id_branch, zero_id_branch);
                  builder.add_transition(source, working_qubit, one_right, transition);
                  return;
             }
             DLF right_subtree  {one  * destination };
             auto transition = synthetize_swta_transition(zero_id_branch, right_subtree);
             builder.add_transition(source, working_qubit, one_right, transition);
        };

        auto add_left_fin_transition = [&one, &builder, &zero_id_branch, &working_qubit_stop, &one_left, &id](State source, State destination) {
             if (destination == id) {
                  auto transition = synthetize_swta_transition(zero_id_branch, zero_id_branch);
                  builder.add_transition(source, working_qubit_stop, one_left, transition);
                  return;
             }
             DLF left_subtree  {one  * destination };
             auto transition = synthetize_swta_transition(left_subtree, zero_id_branch);
             builder.add_transition(source, working_qubit_stop, one_left, transition);

        };

        auto add_right_fin_transition = [&one, &builder, &zero_id_branch, &working_qubit_stop, &one_right, &id](State source, State destination) {
             if (destination == id) {
                  auto transition = synthetize_swta_transition(zero_id_branch, zero_id_branch);
                  builder.add_transition(source, working_qubit_stop, one_right, transition);
                  return;
             }
             DLF right_subtree  {one  * destination };
             auto transition = synthetize_swta_transition(zero_id_branch, right_subtree);
             builder.add_transition(source, working_qubit_stop, one_right, transition);
        };

        add_left_transition(q_c, q_b);
        add_right_transition(q_c, q_b);

        add_left_transition(q_b, q_a);
        add_right_transition(q_b, q_a);
        add_left_fin_transition(q_b, ac);
        add_right_fin_transition(q_b, ac);

        add_left_transition(q_a, q_b);
        add_right_transition(q_a, q_b);

        add_left_transition(ac, c);
        add_right_transition(ac, c);

        {
             auto transition = synthetize_swta_transition(one_leaf_branch, zero_leaf_branch);
             builder.add_transition(c, working_qubit, one_left, transition);
        }

        { // The last read bit is always 0, the postcondition must handle different values of c_n using internal symbols
             // auto transition = synthetize_swta_transition(zero_leaf_branch, zero_leaf_branch);
             // builder.add_transition(c, working_qubit, one_right, transition);
        }

        // The id state
        add_left_transition(id, id);
        add_right_transition(id, id);
        add_left_fin_transition(id, id);
        add_right_fin_transition(id, id);

        std::vector<State> initial_states ({q_c});
        Bit_Set leaf_states (number_of_states, { leaf, id });
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        do_on_debug({
            result.debug_data = new SWTA::Debug_Data();
            auto& state_names = result.debug_data->state_names;

            state_names[q_a]  = "Qa";
            state_names[q_b]  = "Qb";
            state_names[q_c]  = "Qc";
            state_names[id]   = "ID";
            state_names[ac]   = "AC";
            state_names[c]    = "C";
            state_names[leaf] = "leaf";

        });

        return result;
    }

    if (name == Predefined_SWTA_Names::ECC_PRE) {
        State init        =  0;
        State q_w_init    =  1;
        State q_anc_init  =  2;
        State q_w         =  3;
        State q_anc       =  4;
        State s_w_init    =  5;
        State s_anc_init  =  6;
        State s_w         =  7;
        State s_anc       =  8;

        u64 number_of_states = 9;

        Internal_Symbol working_qubit      = 0;
        Internal_Symbol working_qubit_stop = 1;

        Color color_LL = 0;
        Color color_LR = 1;
        Color color_RL = 2;
        Color color_RR = 3;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 4 };
        SWTA::Transition_Builder builder (metadata);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        auto add_transition = [&builder, &working_qubit](State source, Color color, const Def_Linear_Form& left_successor, const Def_Linear_Form& right_successor) {
             DLF left_subtree  { left_successor };
             DLF right_subtree { right_successor };
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(source, working_qubit, color, transition);
        };

        auto add_stop_transition = [&builder, &working_qubit_stop](State source, Color color, const Def_Linear_Form& left_successor, const Def_Linear_Form& right_successor) {
             DLF left_subtree  { left_successor };
             DLF right_subtree { right_successor };
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(source, working_qubit_stop, color, transition);
        };

        add_transition(init, color_LL, one*q_anc_init, one*s_anc_init);

        add_transition(q_anc_init, color_LL, one*q_w_init, zero*q_w_init); // Ancillas are initialized to zero
        add_transition(s_anc_init, color_LL, one*s_w_init, zero*s_w_init);

        // add_transition(q_w_init, color_LL, one*q_anc, zero*q_anc);
        add_transition(q_w_init, color_LR, one*q_anc, zero*q_anc);
        add_transition(q_w_init, color_RL, zero*q_anc, one*q_anc);
        // add_transition(q_w_init, color_RR, zero*q_anc, one*q_anc);

        // add_transition(s_w_init, color_LL, one*s_anc, zero*s_anc);
        add_transition(s_w_init, color_LR, zero*s_anc, one*s_anc);
        add_transition(s_w_init, color_RL, one*s_anc, zero*s_anc);
        // add_transition(s_w_init, color_RR, zero*s_anc, one*s_anc);

        add_transition(q_anc, color_LL, one*q_w, zero*q_w);
        add_transition(s_anc, color_LL, one*s_w, zero*s_w);

        // add_transition(q_w, color_LL, one*q_anc, zero*q_anc);
        add_transition(q_w, color_LR, one*q_anc, zero*q_anc);
        add_transition(q_w, color_RL, zero*q_anc, one*q_anc);
        // add_transition(q_w, color_RR, zero*q_anc, one*q_anc);

        // add_transition(s_w, color_LL, one*s_anc, zero*s_anc);
        add_transition(s_w, color_LR, zero*s_anc, one*s_anc);
        add_transition(s_w, color_RL, one*s_anc, zero*s_anc);
        // add_transition(s_w, color_RR, zero*s_anc, one*s_anc);

        add_stop_transition(q_anc_init, color_LL, one*q_w_init, zero*q_w_init);
        add_stop_transition(s_anc_init, color_LL, one*s_w_init, zero*s_w_init);

        add_stop_transition(q_anc, color_LL, one*q_w, zero*q_w);
        add_stop_transition(s_anc, color_LL, one*s_w, zero*s_w);

        std::vector<State> initial_states ({init});
        Bit_Set leaf_states (number_of_states, { q_w, s_w });
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        do_on_debug({
            result.debug_data = new SWTA::Debug_Data();
            auto& state_names = result.debug_data->state_names;

            state_names[init]       = "init";
            state_names[q_w_init]   = "q_w_init";
            state_names[q_anc_init] = "q_anc_init";
            state_names[q_w]        = "q_w";
            state_names[q_anc]      = "q_anc";
            state_names[s_w_init]   = "s_w_init";
            state_names[s_anc_init] = "s_anc_init";
            state_names[s_w]        = "s_w";
            state_names[s_anc]      = "s_anc";
        });

        return result;
    }

    if (name == Predefined_SWTA_Names::ECC_POST) {
        State init_state    =  0;

        State q_00_w        =  1;
        State q_00_anc      =  2;
        State q_10_w        =  3;
        State q_10_anc      =  4;
        State q_11_anc      =  5;
        State q_01_anc      =  6;
        State q_01_anc_init =  7;
        State q_00_anc_init =  8;
        State q_00_w_init   =  9;
        State q_10_w_init   = 10;

        State s_00_w        = 11;
        State s_00_anc      = 12;
        State s_10_w        = 13;
        State s_10_anc      = 14;
        State s_11_anc      = 15;
        State s_01_anc      = 16;
        State s_01_anc_init = 17;
        State s_00_anc_init = 18;
        State s_00_w_init   = 19;
        State s_10_w_init   = 20;

        u64 number_of_states = 21;

        Internal_Symbol working_qubit = 0;

        Color color_LL = 0;
        Color color_LR = 1;
        Color color_RL = 2;
        Color color_RR = 3;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 1, .number_of_colors = 4 };
        SWTA::Transition_Builder builder (metadata);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        auto add_transition = [&builder, &working_qubit](State source, Color color, const Def_Linear_Form& left_successor, const Def_Linear_Form& right_successor) {
             DLF left_subtree  { left_successor };
             DLF right_subtree { right_successor };
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(source, working_qubit, color, transition);
        };

        // Init state
        add_transition(init_state, color_LL, one*q_00_anc_init, one*s_01_anc_init);

        add_transition(q_00_anc_init, color_LL, one*q_00_w_init, zero*q_00_w_init);
        add_transition(q_01_anc_init, color_LL, one*q_10_w_init, zero*q_00_w_init);

        add_transition(s_00_anc_init, color_LL, one*s_00_w_init, zero*s_00_w_init);
        add_transition(s_01_anc_init, color_LL, one*s_10_w_init, zero*s_00_w_init);

        // add_transition(q_00_w_init, color_LL, one*q_00_anc, zero*q_11_anc);
        add_transition(q_00_w_init, color_LR, one*q_00_anc, zero*q_11_anc);
        add_transition(q_00_w_init, color_RL, zero*q_00_anc, one*q_11_anc);
        // add_transition(q_00_w_init, color_RR, zero*q_00_anc, one*q_11_anc);
        // add_transition(q_10_w_init, color_LL, one*q_10_anc, zero*q_01_anc);
        add_transition(q_10_w_init, color_LR, one*q_10_anc, zero*q_01_anc);
        add_transition(q_10_w_init, color_RL, zero*q_10_anc, one*q_01_anc);
        // add_transition(q_10_w_init, color_RR, zero*q_10_anc, one*q_01_anc);

        // add_transition(s_00_w_init, color_LL, one*s_00_anc,  zero*s_11_anc);
        add_transition(s_00_w_init, color_LR, zero*s_00_anc, one*s_11_anc);
        add_transition(s_00_w_init, color_RL, one*s_00_anc,  zero*s_11_anc);
        // add_transition(s_00_w_init, color_RR, zero*s_00_anc, one*s_11_anc);
        // add_transition(s_10_w_init, color_LL, one*s_10_anc, zero*s_01_anc);
        add_transition(s_10_w_init, color_LR, zero*s_10_anc, one*s_01_anc);
        add_transition(s_10_w_init, color_RL, one*s_10_anc, zero*s_01_anc);
        // add_transition(s_10_w_init, color_RR, zero*s_10_anc, one*s_01_anc);

        // add_transition(q_00_w, color_LL, one*q_00_anc, zero*q_11_anc);
        add_transition(q_00_w, color_LR, one*q_00_anc, zero*q_11_anc);
        add_transition(q_00_w, color_RL, zero*q_00_anc, one*q_11_anc);
        // add_transition(q_00_w, color_RR, zero*q_00_anc, one*q_11_anc);
        // add_transition(q_10_w, color_LL, one*q_10_anc, zero*q_01_anc);
        add_transition(q_10_w, color_LR, one*q_10_anc, zero*q_01_anc);
        add_transition(q_10_w, color_RL, zero*q_10_anc, one*q_01_anc);
        // add_transition(q_10_w, color_RR, zero*q_10_anc, one*q_01_anc);

        // add_transition(s_00_w, color_LL, one*s_00_anc, zero*s_11_anc);
        add_transition(s_00_w, color_LR, zero*s_00_anc, one*s_11_anc);
        add_transition(s_00_w, color_RL, one*s_00_anc, zero*s_11_anc);
        // add_transition(s_00_w, color_RR, zero*s_00_anc, one*s_11_anc);
        // add_transition(s_10_w, color_LL, one*s_10_anc, zero*s_01_anc);
        add_transition(s_10_w, color_LR, zero*s_10_anc, one*s_01_anc);
        add_transition(s_10_w, color_RL, one*s_10_anc, zero*s_01_anc);
        // add_transition(s_10_w, color_RR, zero*s_10_anc, one*s_01_anc);

        add_transition(q_00_anc, color_LL,  one*q_00_w, zero*q_00_w);
        add_transition(q_01_anc, color_LL,  one*q_10_w, zero*q_10_w);
        add_transition(q_10_anc, color_LL, zero*q_00_w,  one*q_00_w);
        add_transition(q_11_anc, color_LL, zero*q_10_w,  one*q_10_w);

        add_transition(s_00_anc, color_LL,  one*s_00_w, zero*s_00_w);
        add_transition(s_01_anc, color_LL,  one*s_10_w, zero*s_10_w);
        add_transition(s_10_anc, color_LL, zero*s_00_w,  one*s_00_w);
        add_transition(s_11_anc, color_LL, zero*s_10_w,  one*s_10_w);

        std::vector<State> initial_states ({init_state});
        Bit_Set leaf_states (number_of_states, { q_00_w, q_10_w, s_00_w, s_10_w }); // We have both 00 and 10 states here, becasue the last qubit can affect only the last ancilla
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        do_on_debug({
            result.debug_data = new SWTA::Debug_Data();
            auto& state_names = result.debug_data->state_names;

            state_names[init_state]    = "init";
            state_names[q_00_w]        = "q_00_w";
            state_names[q_00_anc]      = "q_00_anc";
            state_names[q_10_w]        = "q_10_w";
            state_names[q_10_anc]      = "q_10_anc";
            state_names[q_11_anc]      = "q_11_anc";
            state_names[q_01_anc]      = "q_01_anc";
            state_names[q_01_anc_init] = "q_01_anc_init";
            state_names[q_00_anc_init] = "q_00_anc_init";
            state_names[q_00_w_init]   = "q_00_w_init";
            state_names[q_10_w_init]   = "q_10_w_init";
        });

        return result;
    }

    if (name == Predefined_SWTA_Names::HAMILTONIAN_ALL_BASIS) {
        State q_loop      = 0;
        State q_loop_last = 1;
        State q_leaf      = 2;

        u64 number_of_states = 3;

        Internal_Symbol working_qubit      = 0;
        Internal_Symbol working_qubit_stop = 1;

        Color color_L = 0;
        Color color_R = 1;

        SWTA::Metadata metadata = { .number_of_internal_symbols = 2, .number_of_colors = 2 };

        SWTA::Transition_Builder builder (metadata);

        Def_Coef one  (ACN::ONE());
        Def_Coef zero (ACN::ZERO());

        auto add_transition = [&builder, &working_qubit](State source, Color color, const Def_Linear_Form& left_successor, const Def_Linear_Form& right_successor) {
             DLF left_subtree  { left_successor };
             DLF right_subtree { right_successor };
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(source, working_qubit, color, transition);
        };

        auto add_stop_transition = [&builder, &working_qubit_stop](State source, Color color, const Def_Linear_Form& left_successor, const Def_Linear_Form& right_successor) {
             DLF left_subtree  { left_successor };
             DLF right_subtree { right_successor };
             auto transition = synthetize_swta_transition(left_subtree, right_subtree);
             builder.add_transition(source, working_qubit_stop, color, transition);
        };

        add_transition(q_loop, color_L, one*q_loop, zero*q_loop);
        add_transition(q_loop, color_R, zero*q_loop, one*q_loop);

        add_stop_transition(q_loop, color_L, one*q_loop_last, zero*q_loop_last);
        add_stop_transition(q_loop, color_R, zero*q_loop_last, one*q_loop_last);

        add_transition(q_loop_last, color_L, one*q_leaf, zero*q_leaf);
        add_transition(q_loop_last, color_R, zero*q_leaf, one*q_leaf);

        std::vector<State> initial_states ({q_loop});
        Bit_Set leaf_states (number_of_states, { q_leaf });
        auto transition_fn = builder.build(number_of_states);
        SWTA result (transition_fn, initial_states, leaf_states);

        do_on_debug({
            result.debug_data = new SWTA::Debug_Data();
            result.debug_data->state_names[q_loop] = "q_loop";
            result.debug_data->state_names[q_leaf] = "leaf";
            result.debug_data->state_names[q_loop_last] = "last_bit";
        });

        return result;
    }

    throw std::runtime_error("No definition for the predefined SWTA: " + std::to_string(static_cast<u64>(name)));
}
