#include "swta.hpp"
#include "arith.hpp"
#include "bit_set.hpp"
#include "nfa.hpp"
#include "basics.hpp"

#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>

#define DEBUG_FRONTIER_CONSTRUCTION 0
#define DEBUG_BUILD_AFFINE_PROGRAM  0
#define DEBUG_COLOR_EQUIVALENCE     0

/*
 * operator<< for Linear_Form
 * 用途：
 *   將一個線性形式 (Linear_Form) 以可讀的方式輸出至串流中。線性形式是狀態與
 *   複數係數的組合，例如 (3)*q0 + (-2)*q1 等。
 * 怎麼算：
 *   逐一遍歷 form.components 內的每個分量（component），印出各分量的係數與
 *   狀態索引，並用空格分隔。輸出格式為 Linear_Form{ (coef)*q<state> ... }。
 */
std::ostream& operator<<(std::ostream& os, const Linear_Form& form) {
    os << "Linear_Form{ ";
    for (u64 i = 0; i < form.components.size(); i++) {
        auto& component = form.components[i];
        os << "(" << component.coef << ")" << "*q" << component.state;
        if (i != form.components.size() - 1) {
            os << " ";
        }
    }
    os << "}";
    return os;
}

/*
 * operator<< for SWTA::Transition
 * 用途：
 *   打印單一 SWTA 轉移，清楚地分別顯示左枝與右枝的線性形式。
 * 怎麼算：
 *   轉移由左右兩個線性形式組成。此函數直接輸出這兩個形式的內容。
 */
std::ostream& operator<<(std::ostream& os, const SWTA::Transition& swta_transition) {
    os << "LEFT: " << swta_transition.left << ", "
       << "RIGHT: " << swta_transition.right;
    return os;
}

/*
 * operator<< for SWTA
 * 用途：
 *   完整地輸出整個樹字轉移自動機 (SWTA) 的結構，包括初始態、接收態與所有轉移。
 * 怎麼算：
 *   先印出初始狀態集與葉片轉移狀態集，然後逐狀態、逐符號、逐顏色遍歷轉移表，
 *   將所有存在的轉移打印出來。
 */
std::ostream& operator<<(std::ostream& os, const SWTA& swta) {
    os << "SWTA {\n";
    std::cout << "  initial states: " << swta.initial_states << "\n";
    std::cout << "  final  states:  " << swta.states_with_leaf_transitions.into_vector() << "\n";
    for (u64 state = 0; state < swta.transitions.size(); state++) {
        auto& transitions_from_state = swta.transitions[state];

        for (Internal_Symbol sym = 0; sym < swta.number_of_internal_symbols(); sym++) {
            const auto& transitions_along_symbol = transitions_from_state[sym];
                for (Color color = 0; color < swta.number_of_colors(); color++) {
                auto& transition = transitions_along_symbol[color];

                if (!transition.is_present()) continue;

                os << "  q" << state << " --symbol=" << sym << ", color=" << color << "--> [ " << transition << "]\n";
            }
        }
    }
    os << "}";
    return os;
}


/*
 * operator<< for Linear_Form::Component
 * 用途：
 *   輸出線性形式內單一分量，顯示為「係數 * q<狀態編號>」的形式。
 * 怎麼算：
 *   直接把分量的係數與狀態編號依序列出。
 */
std::ostream& operator<<(std::ostream& os, const Linear_Form::Component& component) {
     os << component.coef << " * q" << component.state;
     return os;
}


/*
 * write_norm_with_subtree_info
 * 用途：
 *   將一個線性形式寫入至串流，並在每個分量後面附加子樹標記（例如 "(L)" 或 "(R)"）。
 *   用於 WTT 轉移的列印，以區分左右子樹的貢獻。
 * 怎麼算：
 *   如果 needs_leading_plus 為真，先輸出一個前導加號。然後逐一遍歷分量，
 *   在每個分量後追加 subtree_info 字符串，並用加號分隔。
 */
void write_norm_with_subtree_info(std::ostream& target, const Linear_Form& form, const char* subtree_info, bool needs_leading_plus) {
    if (needs_leading_plus) {
        target << " + ";
    }

    for (u64 i = 0; i < form.size(); i++) {
        target << form.components[i] << subtree_info;
        if (i < form.size() - 1) target << " + ";
    }
}

/*
 * operator<< for WTT::Transition
 * 用途：
 *   格式化輸出一個加權樹轉移器 (WTT) 的單一轉移。WTT 轉移由四個線性形式組成：
 *   ll、lr、rl、rr，分別對應左子樹及右子樹的轉移。
 * 怎麼算：
 *   分別輸出左子樹部份（ll + lr）與右子樹部份（rl + rr），各部份內的分量
 *   用 write_norm_with_subtree_info 添加子樹標記。
 */
std::ostream& operator<<(std::ostream& os, const WTT::Transition& wtt_transition) {
    os << "LEFT SUBTREE: ";
    bool needs_plus = false; // True, if anything has been written to the output stream
    if (!wtt_transition.ll.empty()) {
        write_norm_with_subtree_info(os, wtt_transition.ll, "(L)", needs_plus);
        needs_plus = true;
    }

    if (!wtt_transition.lr.empty()) {
        write_norm_with_subtree_info(os, wtt_transition.lr, "(R)", needs_plus);
        needs_plus = true;
    }

    if (wtt_transition.ll.empty() && wtt_transition.lr.empty()) {
        os << "0";
    }

    os << "; RIGHT SUBTREE: ";
    needs_plus = false; // Reset

    if (!wtt_transition.rl.empty()) {
        write_norm_with_subtree_info(os, wtt_transition.rl, "(L)", needs_plus);
        needs_plus = true;
    }

    if (!wtt_transition.rr.empty()) {
        write_norm_with_subtree_info(os, wtt_transition.rr, "(R)", needs_plus);
        needs_plus = true;
    }

    if (wtt_transition.rl.empty() && wtt_transition.rr.empty()) {
        os << "0";
    }

    return os;
}

/*
 * operator<< for WTT
 * 用途：
 *   完整地輸出一個加權樹轉移器的結構，包括初始態、葉片態與所有內部轉移。
 * 怎麼算：
 *   先列印初始狀態集與葉片狀態集，然後逐狀態逐符號遍歷轉移表，
 *   印出所有存在的轉移。
 */
std::ostream& operator<<(std::ostream& os, const WTT& wtt) {
    os << "WTT {\n";
    os << "  initial states: " << wtt.initial_states << "\n";
    os << "  leaf states: " << wtt.states_with_leaf_transitions.into_vector() << "\n";

    for (u64 state = 0; state < wtt.number_of_states(); state++) {
        const std::vector<WTT::Transition>& transitions_from_state = wtt.transitions[state];
        for (Internal_Symbol sym = 0; sym < wtt.number_of_internal_symbols(); sym++) {
            auto& transition = transitions_from_state[sym];
            if (!transition.is_present()) continue;
            os << "  " << state << "--(sym=" << sym << ")-->: " << transitions_from_state[sym] << "\n";
        }
    }

    os << "}";
    return os;
}

/*
 * extend_form_with_product_and_node_discoveries
 * 用途：
 *   計算兩個線性形式的外積，同時發現新的狀態對 (inner_state, outer_state)，
 *   並在 worklist 中註冊此狀態對。結果累加至 destination。
 * 怎麼算：
 *   對 first 與 second 中的每對分量 (inner_comp, outer_comp) 計算其係數乘積，
 *   將對應狀態對提交至 worklist_state 進行發現與編號。若目標狀態已存在於
 *   destination 中，則直接加總係數；否則新增分量。
 */
void extend_form_with_product_and_node_discoveries(Linear_Form& destination, const Linear_Form& first, const Linear_Form& second, Worklist_Construction_Context<State_Pair>& worklist_state) {
    for (auto& outer_comp : second.components) {
        for (auto& inner_comp : first.components) {

            const State_Pair* imm_state = nullptr;
            {
                State_Pair state = { .first = inner_comp.state, .second = outer_comp.state };
                imm_state = worklist_state.mark_discovery(state);
            }

            Algebraic_Complex_Number coef = inner_comp.coef * outer_comp.coef;

            bool already_present = false;
            for (auto& component : destination.components) {
                if (component.state == imm_state->handle) {
                    component.coef += coef;
                    already_present = true;
                    break;
                }
            }

            if (already_present) {
                continue;  // We are done here
            }

            Linear_Form::Component new_component(coef, imm_state->handle);
            destination.components.push_back(new_component);
        }
    }
}


/*
 * compose_wtts_sequentially
 * 用途：
 *   順序合成兩個 WTT：第一個 WTT 的輸出作為第二個 WTT 的輸入。結果是一個新的
 *   WTT，代表先執行 first 再執行 second 的複合轉換。
 * 怎麼算：
 *   以狀態對 (first_state, second_state) 作為積自動機的狀態。初始狀態是兩者的
 *   初始狀態對。對每個轉移，使用 extend_form_with_product_and_node_discoveries
 *   計算四個方向（ll、lr、rl、rr）的合成。接收態則是兩者都是葉片態的狀態對。
 */
WTT compose_wtts_sequentially(WTT& first, WTT& second) {

    const u64 num_of_internal_symbols = first.number_of_internal_symbols();
    assert(num_of_internal_symbols == second.number_of_internal_symbols());

    SWTA::Metadata metadata = {.number_of_internal_symbols = num_of_internal_symbols, .number_of_colors = 1};
    WTT_Builder builder(metadata);

    Worklist_Construction_Context<State_Pair> worklist_state;

    State_Pair initial_state = { .first = first.initial_states[0], .second = second.initial_states[0] };
    auto imm_initial_state = worklist_state.mark_discovery(initial_state);
    builder.mark_state_initial(imm_initial_state->handle);

    while (worklist_state.has_more_to_explore()) {
        auto state_pair = worklist_state.extract();

        auto is_first_leaf  = first.states_with_leaf_transitions.get_bit_value(state_pair->first);
        auto is_second_leaf = second.states_with_leaf_transitions.get_bit_value(state_pair->second);

        if (is_first_leaf && is_second_leaf) {
            builder.mark_state_final(state_pair->handle);
        }

        for (Internal_Symbol internal_symbol = 0; internal_symbol < num_of_internal_symbols; internal_symbol++) {
            WTT::Transition& first_transition   = first.transitions[state_pair->first][internal_symbol];
            WTT::Transition& second_transitions = second.transitions[state_pair->second][internal_symbol];

            Linear_Form ll, lr, rl, rr;
            { // ll
                extend_form_with_product_and_node_discoveries(ll, first_transition.ll, second_transitions.ll, worklist_state);
                extend_form_with_product_and_node_discoveries(ll, first_transition.rl, second_transitions.lr, worklist_state);
            }

            { // lr
                extend_form_with_product_and_node_discoveries(lr, first_transition.lr, second_transitions.ll, worklist_state);
                extend_form_with_product_and_node_discoveries(lr, first_transition.rr, second_transitions.lr, worklist_state);
            }

            { // rl
                extend_form_with_product_and_node_discoveries(rl, first_transition.ll, second_transitions.rl, worklist_state);
                extend_form_with_product_and_node_discoveries(rl, first_transition.rl, second_transitions.rr, worklist_state);
            }

            { // rr
                extend_form_with_product_and_node_discoveries(rr, first_transition.lr, second_transitions.rl, worklist_state);
                extend_form_with_product_and_node_discoveries(rr, first_transition.rr, second_transitions.rr, worklist_state);
            }

            WTT::Transition resulting_transition (ll, lr, rl, rr);

            builder.add_transition(state_pair->handle, internal_symbol, resulting_transition);
        }
    }

    WTT result = builder.build(worklist_state.handles.size());

    if (DEBUG) {
        result.debug_data = new WTT::Debug_Data;
        auto& state_names = result.debug_data->state_names;

        auto make_label = [](const WTT& wtt, u64 state) {
            if (!wtt.debug_data) return std::to_string(state);

            if (wtt.debug_data->state_names.contains(state)) return wtt.debug_data->state_names.at(state);
            return std::to_string(state);
        };

        for (auto& [state, handle] : worklist_state.handles) {
            std::string first_label = make_label(first, state.first);
            std::string second_label = make_label(second, state.second);
            std::string label = "(" + first_label + ", " + second_label + ")";

            state_names[handle] = label;
        }
    }
    return result;
}



/*
 * compute_post (SWTA 過載版本)
 * 用途：
 *   給定一個宏態（SWTA 狀態集合）及轉移符號和顏色，計算該宏態經過此轉移後抵達
 *   的下一個宏態。
 * 怎麼算：
 *   對源宏態中的所有狀態，查詢對應符號與顏色的 SWTA 轉移。若任何狀態不存在該轉移，
 *   則下一宏態為空。否則，將所有轉移的左右形式中的所有目標狀態並集起來作為新宏態。
 */
Macrostate compute_post(const Macrostate* macrostate, const SWTA& swta, Color color, Internal_Symbol symbol) {
    Macrostate post(swta.number_of_states());

    for (State state : macrostate->state_names) {
        auto& transitions_from_state = swta.transitions[state];
        auto& transitions_for_sym    = transitions_from_state[symbol];
        auto& transition             = transitions_for_sym[color];

        if (!transition.is_present()) {
            if (DEBUG_FRONTIER_CONSTRUCTION) {
                std::cout << "The state: " << swta.debug_data->state_names[state] << " has no successor!\n";
            }
            post.state_set.clear();
            break;
        }

        for (auto& component : transition.left.components) {
            post.state_set.set_bit(component.state);
        }

        for (auto& component : transition.right.components) {
            post.state_set.set_bit(component.state);
        }
    }

    post.populate_state_names_from_set();

    return post;
}

/*
 * compute_post (WTT 過載版本)
 * 用途：
 *   同樣計算下一宏態，但針對 WTT。WTT 轉移由四個線性形式 (ll, lr, rl, rr) 組成。
 * 怎麼算：
 *   與 SWTA 版本類似，但需要從四個方向的形式中蒐集目標狀態。
 */
Macrostate compute_post(const Macrostate* macrostate, const WTT& wtt, Color color, Internal_Symbol symbol) {
    Macrostate post(wtt.number_of_states());

    for (State state : macrostate->state_names) {
        auto& transitions_from_state = wtt.transitions[state];
        auto& transitions_for_symbol  = transitions_from_state[symbol];

        if (!transitions_for_symbol.is_present()) {
            post.state_set.clear();
            break;
        }

        for (auto& component : transitions_for_symbol.ll.components) {
            post.state_set.set_bit(component.state);
        }

        for (auto& component : transitions_for_symbol.lr.components) {
            post.state_set.set_bit(component.state);
        }

        for (auto& component : transitions_for_symbol.rl.components) {
            post.state_set.set_bit(component.state);
        }

        for (auto& component : transitions_for_symbol.rr.components) {
            post.state_set.set_bit(component.state);
        }
    }

    post.populate_state_names_from_set();

    return post;
}

/*
 * dump_discovered_transitions
 * 用途：
 *   調試輔助函數，將目前發現的所有狀態列印到控制台。
 * 怎麼算：
 *   遍歷 transitions 映射表，印出所有已知狀態編號。
 */
void dump_discovered_transitions(const std::map<State, std::vector<std::vector<State>>>& transitions) {
    std::cout << "Known states: ";
    for (const auto& [state, transitions_from_state]: transitions) {
        std::cout << state << ", ";
    }
    std::cout << "\n";
}

/*
 * initialize_frontier_with_initial_states
 * 用途：
 *   初始化邊界自動機 (frontier automaton) 的構造，將初始宏態加入工作列表與 NFA 建造器。
 * 怎麼算：
 *   若 root < 0，則以所有初始狀態作為初始宏態；否則以指定的 root 狀態作為起點。
 *   將這個初始宏態標記為初始狀態。
 */
void initialize_frontier_with_initial_states(Worklist_Construction_Context<Macrostate>& worklist_state, NFA_Builder& builder, const std::vector<State>& initial_states, u64 total_number_of_states, s64 root) {
    if (root < 0) {
        Macrostate initial_macrostate (total_number_of_states, initial_states);
        worklist_state.mark_discovery(initial_macrostate);
    } else { // Start the construction from the provided root
        std::vector<State> root_states ({static_cast<State>(root)});
        Macrostate initial_macrostate (total_number_of_states, root_states);
        worklist_state.mark_discovery(initial_macrostate);
    }

    builder.mark_state_initial(0);
}

/*
 * make_macrostate_label_using_debug_data
 * 用途：
 *   根據調試資料（state_names 映射表），為一個宏態製作人類可讀的標籤字符串。
 * 怎麼算：
 *   遍歷宏態內的所有狀態，查詢 state_names 映射以取得名稱；若無對應名稱則用狀態編號。
 *   將所有名稱以 {...} 格式連接。
 */
std::string make_macrostate_label_using_debug_data(const Macrostate& macrostate, std::map<State, std::string>& state_names) {
    std::stringstream result;
    result << "{";
    u64 cnt = 0;
    for (auto state: macrostate.state_names) {
        cnt += 1;
        std::string state_label = state_names.contains(state) ? state_names.at(state) : std::to_string(state);
        result << state_label;
        if (cnt < macrostate.state_names.size()) {
            result << ",";
        }
    }
    result << "}";
    return result.str();
}

/*
 * write_macrostate_debug_names (template)
 * 用途：
 *   調試輔助函數，列印宏態內所有狀態的人類可讀名稱。
 * 怎麼算：
 *   遍歷宏態中的所有狀態，查詢 tts.debug_data->state_names 並依序輸出。
 */
template <typename Tree_Transition_System>
void write_macrostate_debug_names(const Macrostate* macrostate, const Tree_Transition_System& tts) {
    for (auto state: macrostate->state_names) {
        std::cout << tts.debug_data->state_names[state] << ", ";
    }
    std::cout << "\n";
}

/*
 * build_frontier_automaton (template)
 * 用途：
 *   建構邊界自動機 (frontier automaton)：一個 NFA，其字母表為顏色集合，狀態為
 *   樹轉移系統的宏態。接受的語言為該系統中所有合法邊界對應的顏色序列。
 * 怎麼算：
 *   自初始宏態開始，使用工作列表 BFS 遍歷。對每個宏態與每種顏色、符號組合，
 *   計算 post 宏態。若 post 非空，則在 NFA 中新增轉移並將 post 放入工作列表。
 *   若宏態中所有狀態都可以進行葉片轉移，則標記為接收態。
 */
template <typename Tree_Transition_System>
NFA build_frontier_automaton(const Tree_Transition_System& tts, s64 root) {

    Worklist_Construction_Context<Macrostate> worklist_state;
    NFA_Builder builder(tts.number_of_colors());

    initialize_frontier_with_initial_states(worklist_state, builder, tts.initial_states, tts.number_of_states(), root);

    u64 color_cnt = tts.number_of_colors();

    while (worklist_state.has_more_to_explore()) {
        auto macrostate = worklist_state.extract();

        if (DEBUG_FRONTIER_CONSTRUCTION) {
            std::cout << "Processing: ";
            write_macrostate_debug_names(macrostate, tts);
        }

        if (tts.states_with_leaf_transitions.is_superset(macrostate->state_set)) { // All of the states in macrostate can make a leaf transition
            builder.mark_state_final(macrostate->handle);
        }

        for (Internal_Symbol internal_symbol = 0; internal_symbol < tts.number_of_internal_symbols(); internal_symbol++) {
            for (Color color = 0; color < color_cnt; color++) {
                Macrostate imm_post = compute_post(macrostate, tts, color, internal_symbol);

                if (DEBUG_FRONTIER_CONSTRUCTION) {
                    std::cout << "Successor color=" << color << ", symbol=" << internal_symbol << ":";
                    write_macrostate_debug_names(&imm_post, tts);
                }

                if (imm_post.empty()) {
                    continue;
                }

                auto post = worklist_state.mark_discovery(imm_post);
                builder.add_transition(macrostate->handle, color, post->handle);
            }
        }
    }

    NFA result = builder.build(worklist_state.handles.size());

    if (DEBUG) {
        result.debug_data = new NFA::Debug_Data;

        if (tts.debug_data == nullptr) {
            auto& state_names = tts.debug_data->state_names;

            for (auto& [macrostate, handle] : worklist_state.handles) {
                result.debug_data->state_names.emplace(handle, macrostate.to_string());
            }
        } else {
            for (auto& [macrostate, handle] : worklist_state.handles) {
                auto macrostate_label = make_macrostate_label_using_debug_data(macrostate, tts.debug_data->state_names);
                result.debug_data->state_names.emplace(handle, macrostate_label);
            }
        }
    };

    return result;
}

template
NFA build_frontier_automaton<SWTA>(const SWTA& tts, s64 root = -1);

template
NFA build_frontier_automaton<WTT>(const WTT& tts, s64 root = -1);

/*
 * WTT::does_state_accept_trees_for_any_colored_sequence
 * 用途：
 *   判定 WTT 從指定狀態開始是否能接受所有可能顏色序列對應的樹（全稱性檢測）。
 * 怎麼算：
 *   構建 frontier automaton，將其確定化，然後檢查所有狀態是否都是接收態。
 */
bool WTT::does_state_accept_trees_for_any_colored_sequence(State state) const {
    NFA accepted_colored_sequences_abstraction = build_frontier_automaton(*this);
    NFA determinized_abstraction = accepted_colored_sequences_abstraction.determinize();

    return determinized_abstraction.is_every_state_accepting();
}

enum class State_Universality_Status : u8 {
    UNKNOWN = 0,
    UNIVERSAL = 1,
    NONUNIVERSAL = 2,
};

/*
 * can_component_be_removed
 * 用途：
 *   判定一個線性形式的分量是否可以安全移除。若該分量係數為零且其目標狀態
 *   可接受任意顏色序列，則可移除。
 * 怎麼算：
 *   先查詢快取以判定該狀態是否具全稱性。若係數為零且狀態全稱，則返回 true。
 */
bool can_component_be_removed(Linear_Form::Component& component, const WTT& wtt, std::vector<State_Universality_Status>& cache) {

    if (cache[component.state] != State_Universality_Status::UNKNOWN) {
        bool is_coef_zero = component.coef.is_zero();

        if (is_coef_zero && cache[component.state] == State_Universality_Status::UNIVERSAL) {
            return true;
        }
        return false;
    }

    bool is_coef_zero = component.coef.is_zero();
    bool is_state_universal = wtt.does_state_accept_trees_for_any_colored_sequence(component.state);

    cache[component.state] = is_state_universal ? State_Universality_Status::UNIVERSAL : State_Universality_Status::NONUNIVERSAL;

    return is_coef_zero && is_state_universal;
}

/*
 * remove_zeros_from_form
 * 用途：
 *   從線性形式中移除可移除的零分量，藉由兩指針交換法壓縮列表。
 * 怎麼算：
 *   使用 zero_idx 指向可移除的分量，nonzero_idx 指向應保留的分量。從雙端掃過，
 *   不斷交換，將保留分量壓縮至前端，最後截斷陣列。
 */
void remove_zeros_from_form(const WTT& wtt, Linear_Form& form, std::vector<State_Universality_Status>& cache) {
    s64 nonzero_idx = form.size() - 1;
    s64 zero_idx    = 0;

    if (nonzero_idx == zero_idx) { // there is only one element
        auto& component = form.components[zero_idx];

        if (can_component_be_removed(component, wtt, cache)) {
            form.components.clear();
        }
        return;
    }

    while (zero_idx < nonzero_idx) {
        // Search for the next zero slot that needs to be filled
        for (; zero_idx < form.size(); zero_idx++) {
            if (can_component_be_removed(form.components[zero_idx], wtt, cache)) {
                break;
            }
        }

        // Search for the next coef to fill the slot from the back
        for (; nonzero_idx >= 0; nonzero_idx--) {
            if (!can_component_be_removed(form.components[zero_idx], wtt, cache)) {
                break;
            }
        }
        if (zero_idx < nonzero_idx) break;

        form.components[zero_idx].swap(form.components[nonzero_idx]);
    }

    if (zero_idx >= form.size()) {
        form.components.resize(nonzero_idx);
    }
}

/*
 * WTT::remove_zeros_from_transitions
 * 用途：
 *   遍歷 WTT 的所有轉移，逐一對四個方向的線性形式（ll、lr、rl、rr）呼叫
 *   remove_zeros_from_form，以精簡轉移表。
 * 怎麼算：
 *   建立狀態全稱性快取，然後逐符號、逐轉移進行去零處理。
 */
void WTT::remove_zeros_from_transitions() {
    std::vector<State_Universality_Status> cache;
    cache.resize(this->number_of_states());

    for (Internal_Symbol internal_symbol = 0; internal_symbol < this->transitions.size(); internal_symbol++) {
        std::vector<Transition>& transitions_for_symbol = this->transitions[internal_symbol];

        for (auto& transition : transitions_for_symbol) {
            remove_zeros_from_form(*this, transition.ll, cache);
            remove_zeros_from_form(*this, transition.lr, cache);
            remove_zeros_from_form(*this, transition.rl, cache);
            remove_zeros_from_form(*this, transition.rr, cache);
        }
    }
}

/*
 * compose_swta_transition_with_wtt
 * 用途：
 *   將一個 SWTA 轉移與 WTT 轉移合成，產生新的 SWTA 轉移。此操作代表先套用
 *   WTT 轉換再進行 SWTA 轉移。
 * 怎麼算：
 *   左形式 = SWTA.left × WTT.ll + SWTA.right × WTT.lr
 *   右形式 = SWTA.right × WTT.rr + SWTA.left × WTT.rl
 *   使用 extend_form_with_product_and_node_discoveries 計算外積。
 */
SWTA::Transition compose_swta_transition_with_wtt(const SWTA::Transition& swta_transition, const WTT::Transition& wtt_transition, Worklist_Construction_Context<State_Pair>& worklist_state) {
    SWTA::Transition result;

    extend_form_with_product_and_node_discoveries(result.left, swta_transition.left,  wtt_transition.ll, worklist_state);
    extend_form_with_product_and_node_discoveries(result.left, swta_transition.right, wtt_transition.lr, worklist_state);

    extend_form_with_product_and_node_discoveries(result.right, swta_transition.right, wtt_transition.rr, worklist_state);
    extend_form_with_product_and_node_discoveries(result.right, swta_transition.left, wtt_transition.rl, worklist_state);

    return result;
}

/*
 * operator<< for SWTA::Transition_Builder
 * 用途：
 *   輸出 SWTA 轉移建造器的內部狀態，用於調試。
 * 怎麼算：
 *   遍歷 builder.transitions，印出所有儲存的轉移。
 */
std::ostream& operator<<(std::ostream& out, const SWTA::Transition_Builder& builder) {
    std::cout << "Builder (stored transitions): ";
    for (auto& [handle, transitions] : builder.transitions) {
        out << transitions << "\n";
    }

    return out;
}

/*
 * apply_wtt_to_swta
 * 用途：
 *   將一個 WTT 套用到 SWTA，生成新 SWTA。結果是一個積自動機，狀態為 (SWTA_state, WTT_state) 對。
 * 怎麼算：
 *   初始狀態為 (SWTA.initial, WTT.initial)。對每個積狀態與每個符號、顏色組合，
 *   使用 compose_swta_transition_with_wtt 計算新轉移。接收態為兩者都是葉片態的狀態對。
 */
SWTA apply_wtt_to_swta(const SWTA& swta, const WTT& wtt) {
    Worklist_Construction_Context<State_Pair> worklist_state;

    Bit_Set leaf_states (1);
    std::vector<State> initial_states;
    SWTA::Transition_Builder transition_builder ( swta.get_metadata() );


    State_Pair initial_state = { .first = swta.initial_states[0], .second = wtt.initial_states[0] };
    auto imm_initial_state = worklist_state.mark_discovery(initial_state);
    initial_states.push_back(imm_initial_state->handle);

    while (worklist_state.has_more_to_explore()) {
        auto product_state = worklist_state.extract();

        bool is_swta_state_leaf = swta.states_with_leaf_transitions.get_bit_value(product_state->first);
        bool is_wtt_state_leaf  = wtt.states_with_leaf_transitions.get_bit_value(product_state->second);

        if (is_swta_state_leaf && is_wtt_state_leaf) {
            leaf_states.grow_and_set_bit(product_state->handle);
        }

        for (Internal_Symbol sym = 0; sym < swta.number_of_internal_symbols(); sym++) {
            const auto& wtt_transition = wtt.transitions[product_state->second][sym];

            for (Color color = 0; color < swta.number_of_colors(); color++) {
                const auto& swta_transition = swta.transitions[product_state->first][sym][color];

                if (!swta_transition.is_present()) continue;
                if (!wtt_transition.is_present())  continue;

                auto result_form = compose_swta_transition_with_wtt(swta_transition, wtt_transition, worklist_state);
                transition_builder.add_transition(product_state->handle, sym, color, result_form);
            }
        }

    }

    leaf_states.grow(worklist_state.handles.size());
    auto transition_fn = transition_builder.build(worklist_state.handles.size());
    SWTA result (transition_fn, initial_states, leaf_states);

    if (DEBUG) {
        result.debug_data = new SWTA::Debug_Data;
        for (auto& [state, handle] : worklist_state.handles) {
            std::string swta_label = swta.debug_data->state_names.contains(state.first) ? swta.debug_data->state_names.at(state.first) : std::to_string(state.first);
            std::string wtt_label = wtt.debug_data->state_names.contains(state.second) ? wtt.debug_data->state_names.at(state.second) : std::to_string(state.second);

            result.debug_data->state_names[handle] = "(" + swta_label + ", " + wtt_label + ")";
        }
    }


    return result;
}

/*
 * put_form_into_matrix
 * 用途：
 *   將一個線性形式填入矩陣的指定行，使得形式中每個分量 (coef, state) 對應
 *   矩陣的 (source_state, state) 位置。
 * 怎麼算：
 *   遍歷形式內的所有分量，將其係數寫入矩陣相應位置。
 */
void put_form_into_matrix(ACN_Matrix& matrix, State source_state, const Linear_Form& form) {
    for (const auto& component : form.components) {
        matrix.set(source_state, component.state, component.coef);
    }
}

/*
 * extract_transition_matrices_from_swta
 * 用途：
 *   從 SWTA 中提取所有轉移，將其表示為矩陣形式以供 Affine_Program 使用。
 *   對每個 (符號, 顏色) 對，分別提取左、右方向的轉移矩陣。
 * 怎麼算：
 *   逐符號、逐顏色遍歷 SWTA 狀態，對每個狀態取其轉移的左、右線性形式，
 *   使用 put_form_into_matrix 填入對應矩陣。
 */
std::pair<Affine_Program<Branch_Selector>::Symbol_Handles, Affine_Program<Branch_Selector>::Symbol_Store>
extract_transition_matrices_from_swta(const SWTA& swta) {
    Affine_Program<Branch_Selector>::Symbol_Store symbol_store;
    std::map<Branch_Selector, u64> symbol_handles;

    for (Internal_Symbol internal_sym = 0; internal_sym < swta.number_of_internal_symbols(); internal_sym++) {
        for (Color color = 0; color < swta.number_of_colors(); color++) {
            Branch_Selector left_symbol  (internal_sym, color, Subtree_Tag::LEFT);
            Branch_Selector right_symbol (internal_sym, color, Subtree_Tag::RIGHT);

            ACN_Matrix left_matrix  (swta.number_of_states(), swta.number_of_states());
            ACN_Matrix right_matrix (swta.number_of_states(), swta.number_of_states());

            for (State state = 0; state < swta.number_of_states(); state++) {
                auto& transition = swta.get_transition(state, internal_sym, color);

                put_form_into_matrix(left_matrix,  state, transition.left);
                put_form_into_matrix(right_matrix, state, transition.right);
            }

            const auto [left_insert_pos, left_was_inserted]   = symbol_handles.emplace(left_symbol,  symbol_handles.size());
            const auto [right_insert_pos, right_was_inserted] = symbol_handles.emplace(right_symbol, symbol_handles.size());

            symbol_store.push_back( {left_symbol, left_matrix} );
            symbol_store.push_back( {right_symbol, right_matrix} );
        }
    }

    return {symbol_handles, symbol_store};
}


struct AP_State_Info {
    Macrostate macrostate; // What states are active in a branch
    State color_sym_abstraction_state;

    AP_State_Info(u64 swta_state_cnt, State color_sym_abstr_state) :
        macrostate(swta_state_cnt),
        color_sym_abstraction_state(color_sym_abstr_state) {}

    AP_State_Info(const Macrostate& macrostate, State color_sym_abstr_state) :
        macrostate(macrostate),
        color_sym_abstraction_state(color_sym_abstr_state) {}

    bool operator<(const AP_State_Info& other) const {
        INSERT_LEX_LT_CODE(color_sym_abstraction_state, other.color_sym_abstraction_state);
        return this->macrostate < other.macrostate;
    }
};

struct AP_Build_Context {
    const SWTA& swta;
    const Color_Symbol_Abstraction               color_sym_abstraction;
    Affine_Program_Builder<Branch_Selector>      builder;
    Worklist_Construction_Context<AP_State_Info> worklist_state;

    AP_Build_Context (const SWTA& swta, const Color_Symbol_Abstraction& abstraction, const Affine_Program<Branch_Selector>::Symbol_Handles& handles, const Affine_Program<Branch_Selector>::Symbol_Store& store) :
        swta(swta),
        color_sym_abstraction(abstraction),
        builder(handles, store)
    {}
};

/*
 * add_transition_and_add_to_worklist_if_new
 * 用途：
 *   在親和程式 (Affine Program) 中新增一條轉移，若目標狀態為新發現狀態則加入工作列表。
 * 怎麼算：
 *   查詢 post 是否已在 worklist_state 的 handles 中。若新，則分配編號並入工作列表；
 *   若已存，則使用既有編號。最後在 builder 中紀錄這條轉移。
 */
void add_transition_and_add_to_worklist_if_new(const AP_State_Info& source, AP_State_Info& post, Branch_Selector& symbol_info, AP_Build_Context& context) {

    post.macrostate.handle = context.worklist_state.handles.size();
    auto [insert_pos, was_inserted] = context.worklist_state.handles.emplace(post, post.macrostate.handle);

    if (was_inserted) {
        context.worklist_state.worklist.push_back(&insert_pos->first);
        const_cast<AP_State_Info*>(&insert_pos->first)->macrostate.populate_state_names_from_set();
    } else {
        post.macrostate.handle = insert_pos->second;
    }

    u64 symbol_handle = context.builder.symbol_handles.at(symbol_info);
    context.builder.add_transition(source.macrostate.handle, symbol_handle, post.macrostate.handle);
}

/*
 * extend_affine_program_with_post
 * 用途：
 *   給定親和程式中的源狀態及符號/顏色資訊，計算左右兩個下一宏態，並新增轉移。
 * 怎麼算：
 *   先檢查源宏態內所有狀態是否都定義了該符號/顏色的轉移。若有缺失則返回。
 *   然後查詢顏色-符號抽象自動機以確定下一抽象狀態。最後對左右方向分別計算
 *   下一宏態並新增轉移。
 */
void extend_affine_program_with_post(const AP_State_Info& source_ap_state, Branch_Selector& symbol_info, AP_Build_Context& context) {
    for (State state : source_ap_state.macrostate.state_names) {
        auto& transition = context.swta.get_transition(state, symbol_info.symbol, symbol_info.color);

        if (!transition.is_present()) return; // There is no post to compute, all of the states have to have a transition along the symbol
    }

    Color_Symbol color_symbol = {.color = symbol_info.color, .symbol = symbol_info.symbol};
    u64 abstraction_color_sym_handle = context.color_sym_abstraction.symbol_handles.at(color_symbol);
    auto& abstraction_successors = context.color_sym_abstraction.abstraction.transitions[source_ap_state.color_sym_abstraction_state][abstraction_color_sym_handle];
    if (abstraction_successors.empty()) return;

    State abstraction_post = abstraction_successors[0];

    AP_State_Info left_post (context.swta.number_of_states(), abstraction_post);
    AP_State_Info right_post (context.swta.number_of_states(), abstraction_post);

    for (State state : source_ap_state.macrostate.state_names) {
        auto& transition = context.swta.get_transition(state, symbol_info.symbol, symbol_info.color);

        for (auto& component : transition.left.components) {
            left_post.macrostate.state_set.set_bit(component.state);
        }

        for (auto& component : transition.right.components) {
            right_post.macrostate.state_set.set_bit(component.state);
        }
    }

    symbol_info.tag = Subtree_Tag::LEFT;
    add_transition_and_add_to_worklist_if_new(source_ap_state, left_post, symbol_info, context);

    symbol_info.tag = Subtree_Tag::RIGHT;
    add_transition_and_add_to_worklist_if_new(source_ap_state, right_post, symbol_info, context);
}

/*
 * build_affine_program
 * 用途：
 *   從 SWTA 與顏色-符號抽象構造親和程式 (Affine Program)。親和程式以矩陣與
 *   向量表示 SWTA 的轉移，以便進行數值等價性檢查。
 * 怎麼算：
 *   先從 SWTA 提取轉移矩陣，然後自初始狀態開始進行工作列表 BFS。對每個親和
 *   狀態及每個符號/顏色組合，計算左右下一宏態。若該宏態內所有狀態都是葉片態
 *   且顏色符號抽象中也到達接收態，則標記為接收態。
 */
Affine_Program<Branch_Selector> build_affine_program(const SWTA& swta, const Color_Symbol_Abstraction& color_sym_abstraction) {
    auto [symbol_handles, symbol_store] = extract_transition_matrices_from_swta(swta);
    AP_Build_Context build_context (swta, color_sym_abstraction, std::move(symbol_handles), std::move(symbol_store));

    {
        Macrostate init_macrostate (swta.number_of_states(), swta.initial_states);
        init_macrostate.handle = 0;
        AP_State_Info initial_ap_state (init_macrostate, build_context.color_sym_abstraction.abstraction.initial_states[0]);
        auto [insert_pos, was_inserted] = build_context.worklist_state.handles.emplace(initial_ap_state, 0);

        build_context.worklist_state.worklist.push_back(&insert_pos->first);
        build_context.builder.mark_state_initial(0);
    }

    while (build_context.worklist_state.has_more_to_explore()) {
        auto current_ap_state = build_context.worklist_state.extract();

        if (DEBUG_BUILD_AFFINE_PROGRAM) {
            std::cout << "Processing: ";
            for (auto state : current_ap_state->macrostate.state_names) {
                std::cout << build_context.swta.debug_data->state_names[state] << ", ";
            }
            std::cout << "; handle=" << current_ap_state->macrostate.handle << "\n";
        }

        bool swta_branch_contains_only_leaves = build_context.swta.states_with_leaf_transitions.is_superset(current_ap_state->macrostate.state_set);
        bool swta_accepts_in_remaining_branches = build_context.color_sym_abstraction.abstraction.final_states.get_bit_value(current_ap_state->color_sym_abstraction_state);

        if (swta_accepts_in_remaining_branches && swta_accepts_in_remaining_branches) {
            build_context.builder.mark_state_final(current_ap_state->macrostate.handle);
        }

        for (Internal_Symbol symbol = 0; symbol < swta.number_of_internal_symbols(); symbol++) {
            for (Color color = 0; color < swta.number_of_colors(); color++) {
                Branch_Selector symbol_info(symbol, color, Subtree_Tag::NONE);

                extend_affine_program_with_post(*current_ap_state, symbol_info, build_context);
            }
        }
    }

    Affine_Program program = build_context.builder.build(build_context.worklist_state.handles.size());

    do_on_debug({
        program.debug_data = new Affine_Program<Branch_Selector>::Debug_Data;
        for (auto& [ap_state, handle] : build_context.worklist_state.handles) {
            program.debug_data->state_names.emplace(handle, ap_state.macrostate.to_string());
        }
    });

    return program;
}

/*
 * write_affine_program_into_dot (template)
 * 用途：
 *   將親和程式輸出為 DOT 格式，可用 Graphviz 視覺化。
 * 怎麼算：
 *   編寫 DOT 圖描述：節點、初始指針、邊（標記為符號），接收態用綠色表示。
 */
template <typename T>
void write_affine_program_into_dot(std::ostream& stream, const Affine_Program<T>& program) {
    stream << "digraph Affine_Program {\n";

    stream << "  qInit [shape=none, label=\"\"]\n";
    for (State state = 0; state < program.number_of_states(); state++) {
        stream << "  q" << state;
        if (program.debug_data != nullptr && program.debug_data->state_names.contains(state)) {
            const char* color = program.final_states.get_bit_value(state) ? "green" : "black";
            stream << " [label=\"" << program.debug_data->state_names.at(state) << "\", color=" << color << "]";
        }
        stream << "\n";
    }

    stream << "  qInit -> q" << program.initial_state << "\n";

    for (State state = 0; state < program.number_of_states(); state++) {
        for (u64 symbol = 0; symbol < program.transition_fn[state].size(); symbol++) {
            typename Affine_Program<T>::Transition_Symbol symbol_desc = program.symbol_store.at(symbol);
            for (State destination : program.transition_fn[state][symbol]) {
                stream << "  q" << state << " -> q" << destination << " [label=\"" << symbol_desc.info << "\"]\n";
            }
        }
    }

    stream << "}\n";
}


/*
 * operator<< for Branch_Selector
 * 用途：
 *   輸出分支選擇符 (Branch_Selector)，以緊湊格式表示「顏色-符號-方向」。
 * 怎麼算：
 *   直接組合顏色、符號與方向標記輸出。
 */
std::ostream& operator<<(std::ostream& os, const Branch_Selector& info) {
    const char* tag_name = info.tag == Subtree_Tag::LEFT ? "L" : "R";
    // os << "(symbol=" << info.symbol << ", color=" << info.color << ", " << tag_name << ")";
    os << info.color << "-" << info.symbol << "-" << tag_name;
    return os;
}

/*
 * build_color_language_abstraction
 * 用途：
 *   構建 SWTA 的顏色語言抽象：一個 DFA，其字母表為顏色集合，語言為該 SWTA
 *   所有合法邊界對應的顏色序列。
 * 怎麼算：
 *   先建構邊界自動機 (NFA)，再確定化並完備化，得到結果 DFA。
 */
NFA build_color_language_abstraction(const SWTA& swta) {
    auto abstraction_nfa = build_frontier_automaton(swta);
    // abstraction_nfa.write_dot(std::cout);
    auto abstraction_dfa = abstraction_nfa.determinize();
    abstraction_dfa.complete();
    return std::move(abstraction_dfa);
}

/*
 * are_two_swtas_color_equivalent
 * 用途：
 *   判定兩個 SWTA 是否「顏色等價」。顏色等價意味著它們接受相同的顏色語言並且
 *   轉移矩陣在代數意義上也相等（通過最終非零檢查）。
 * 怎麼算：
 *   1. 先檢查兩者的顏色語言是否相等（若不等則返回 false）。
 *   2. 若相等，分別構建親和程式。
 *   3. 建構兩個親和程式的積自動機。
 *   4. 運行非零最終狀態檢查：若積無法到達非零最終狀態，則等價；否則不等價。
 */
bool are_two_swtas_color_equivalent(const SWTA& first, const SWTA& second) {
    NFA first_swta_abstraction  = build_color_language_abstraction(first);

    NFA second_swta_abstraction = build_color_language_abstraction(second);
    bool are_colored_languages_equivalent = are_two_complete_dfas_equivalent(first_swta_abstraction, second_swta_abstraction);

    if (!are_colored_languages_equivalent) {
        do_on_debug({
            std::cout << "Provided SWTAs do not have equal colored languages:" << "\n";
            first_swta_abstraction.write_dot(std::cout);
            std::cout << "\n------------------------" << "\n";
            second_swta_abstraction.write_dot(std::cout);
        });
        return false;
    }

    auto first_color_sym_abstraction = build_color_internal_symbol_abstraction(first);
    auto first_program  = build_affine_program(first, first_color_sym_abstraction);

    auto second_color_sym_abstraction = build_color_internal_symbol_abstraction(second);
    auto second_program = build_affine_program(second, second_color_sym_abstraction);

    auto product = build_colored_product_of_affine_programs(first_program, second_program, first_swta_abstraction);

    Underlying_SWTA_Info first_swta_info (first);
    Underlying_SWTA_Info second_swta_info (second);
    Underlying_SWTA_Info_Pair swta_info_pair (first_swta_info, second_swta_info);

    bool are_equivalent = !does_affine_program_reach_nonzero_final_states(product, swta_info_pair);

    return are_equivalent;
}

struct State_Delta_Info {
    u32 new_vector_count = 0;  // 該狀態新發現多少個向量
    s32 new_vector_start = 0;  // 新發現的向量在矩陣中的起始位置
};

// 調試用途：記錄向量傳播過程的詳細資訊
struct Propagation_Info {
    State source;
    State target;
    ACN_Matrix source_row;
    ACN_Matrix propagated_row;
    ACN_Matrix propagated_row_after_insert;
    ACN_Matrix symbol_matrix;
    Branch_Product_Sym symbol_info;
};

struct Affine_Program_Propagation_Context {
    const Affine_Program<Branch_Product_Sym>& program;
    ACN_Matrix                          final_vector;
    const Underlying_SWTA_Info_Pair&    underlying_swtas;
    std::vector<State_Delta_Info>       state_deltas;
    std::vector<State>                  worklist;
    std::vector<ACN_Matrix>             state_vector_spaces;
    u64                                 state_space_dimension;
    std::vector<Propagation_Info>       propagation_log;
    s64                                 final_state_with_nonzero = -1;
    bool                                check_final_state_early = true;
};

struct Propagation_Stats {
    u32 propagation_cnt = 0;
};

/*
 * filter_propagations_for_those_that_affect_state
 * 用途：
 *   從傳播記錄中篩選出影響到特定最終狀態的傳播路徑（反向追蹤）。
 * 怎麼算：
 *   從給定狀態開始，逆序遍歷傳播日誌，若傳播的目標在有趣集合中，
 *   則將源狀態加入有趣集合，並記錄此傳播。
 */
std::vector<Propagation_Info*> filter_propagations_for_those_that_affect_state(std::vector<Propagation_Info>& propagations, State state_of_interest) {
    std::set<State> interesting_states;
    interesting_states.insert(state_of_interest);

    std::vector<Propagation_Info*> filtered_propagations;

    for (auto it = propagations.rbegin(); it != propagations.rend(); ++it) {
        auto& prop_info = *it;

        if (interesting_states.contains(prop_info.target)) {
            interesting_states.insert(prop_info.source);
            filtered_propagations.push_back(&prop_info);
        }
    }

    return filtered_propagations;
}

/*
 * write_propagation_info
 * 用途：
 *   列印單一傳播記錄，包括源狀態、目標狀態、轉移符號及進入/離開行向量。
 * 怎麼算：
 *   格式化輸出傳播資訊，若提供 state_names 則使用人類可讀名稱；
 *   並將行向量轉換為狀態形式列印。
 */
void write_propagation_info(std::map<State, std::string> state_names, Propagation_Info& info, Affine_Program_Propagation_Context& ctx, bool write_state_names = false) {
    std::cout << "Propagated from "
              << info.source;
    if (write_state_names) {
        std::cout << " aka " << state_names.at(info.source);
    }
    std::cout << " to " << info.target;
    if (write_state_names) {
        std::cout << " aka " << state_names.at(info.target);
    }
    std::cout <<" along " << info.symbol_info << "\n";

    auto translate_into_states = [&state_names, &ctx](ACN_Matrix& row) {
        for (u32 state = 0; state < row.width; state++) {
            if (!row.at(0, state).is_zero()) {
                std::cout << "  :> Q" << state << " = " << row.at(0, state) << "\n";
            }
        }
    };

    std::cout << "  row entering pipe : \n";
    translate_into_states(info.source_row);
    std::cout << "\n";

    std::cout << "  row exiting pipe : \n";
    translate_into_states(info.propagated_row);
    std::cout << "\n";

    // std::cout << "  row exiting pipe(0):" << info.propagated_row << "\n";
    // std::cout << "  row exiting pipe(1):" << info.propagated_row_after_insert << "\n";
    // std::cout << "  pipe matrix:" << propagation_ptr->symbol_matrix << "\n";
}

/*
 * dump_propagations（指針版本）
 * 用途：
 *   列印一組傳播記錄（指針向量版本）。
 * 怎麼算：
 *   遍歷向量，對每個傳播呼叫 write_propagation_info。
 */
void dump_propagations(std::vector<Propagation_Info*> propagations, Affine_Program_Propagation_Context& context) {
    std::map<State, std::string>& state_names = context.program.debug_data->state_names;
    for (auto propagation_ptr : propagations) {
        write_propagation_info(state_names, *propagation_ptr, context);
    }
}

/*
 * dump_propagations（值版本）
 * 用途：
 *   列印一組傳播記錄（值向量版本）。
 * 怎麼算：
 *   遍歷向量，對每個傳播呼叫 write_propagation_info。
 */
void dump_propagations(std::vector<Propagation_Info> propagations, Affine_Program_Propagation_Context& context) {
    std::map<State, std::string>& state_names = context.program.debug_data->state_names;
    for (auto propagation_ptr : propagations) {
        write_propagation_info(state_names, propagation_ptr, context);
    }
}

/*
 * does_state_vector_space_contain_nonzeros
 * 用途：
 *   檢查某狀態的列向量空間與最終向量相乘是否產生非零值。用於判定最終狀態是否可達。
 * 怎麼算：
 *   計算 state_matrix × final_vector，檢查結果是否含非零元素。
 */
bool does_state_vector_space_contain_nonzeros(Affine_Program_Propagation_Context& context, const ACN_Matrix& final_vector, State state) {
    auto& current_state_matrix  = context.state_vector_spaces[state];
    auto pontential_leaf_values = current_state_matrix * final_vector; // [n x 1] vector

    if (!pontential_leaf_values.contains_only_zeros()) { // A final state is reachable with non-zero value
        context.final_state_with_nonzero = state;
        return true;
    };

    return false;
}

/*
 * propagate_vector_from_state_matrix_to_successors
 * 用途：
 *   將狀態的新發現行向量進行矩陣乘法並傳播至後繼狀態。若後繼產生新的行向量，
 *   則將後繼加入工作列表。
 * 怎麼算：
 *   對源狀態新發現的每一行，乘以轉移符號對應的轉移矩陣，將結果插入後繼的行
 *   列向量空間（以行梯形形式）。若後繼為接收態，則進行早期非零檢查。
 */
void propagate_vector_from_state_matrix_to_successors(Affine_Program_Propagation_Context& context, State_Delta_Info& old_state_info, State current_state, u64 symbol) {
    auto& successors_along_symbol = context.program.transition_fn[current_state][symbol];

    u64 last_row_idx_exc = old_state_info.new_vector_start + old_state_info.new_vector_count;

    auto& current_state_matrix  = context.state_vector_spaces[current_state];

    const auto& transition_matrix = context.program.symbol_store[symbol].matrix;

    for (State successor : successors_along_symbol) {
        auto& successor_delta = context.state_deltas[successor];
        if (successor_delta.new_vector_start >= context.state_space_dimension) {
            continue;  // The matrix for the successor is already saturated
        }

        for (u64 row_idx = old_state_info.new_vector_start; row_idx < last_row_idx_exc; row_idx++) {
            auto& target_matrix = context.state_vector_spaces[successor];

            auto row = current_state_matrix.extract_row(row_idx); // @Optimize: maybe we do not need to do a copy here since we already will be allocating during matrix multiply
            auto propagated_row = row * transition_matrix;

            if (propagated_row.contains_only_zeros()) continue;

            s64 new_row_position = add_row_to_row_echelon_matrix_no_copy(target_matrix, propagated_row);

            if (true) {
                Propagation_Info propagation_info = {
                    .source = current_state,
                    .target = successor,
                    .source_row = row,
                    .propagated_row = row * transition_matrix,
                    .propagated_row_after_insert = propagated_row,
                    .symbol_matrix = transition_matrix,
                    .symbol_info = context.program.symbol_store[symbol].info
                };

                context.propagation_log.push_back(propagation_info);

                if (context.check_final_state_early) {
                    bool is_target_state_final = context.program.final_states.get_bit_value(successor);
                    if (is_target_state_final) {
                        if (does_state_vector_space_contain_nonzeros(context, context.final_vector, successor)) {
                            return;
                        }
                    }
                }
            };

            if (new_row_position < 0) {
                continue;
            }

            if (new_row_position == context.state_space_dimension) {
                successor_delta.new_vector_start = context.program.number_of_states();
                break;  // The matrix became satured, no need to try propagate anything
            }

            if (successor_delta.new_vector_count == 0) {
                context.worklist.push_back(successor);
            }

            successor_delta.new_vector_count += 1;
        }
    }
}

/*
 * does_affine_program_reach_nonzero_final_states
 * 用途：
 *   **核心演算法**：檢查積親和程式是否能到達非零的最終狀態。此檢查是色等價性
 *   判定的關鍵步驟。若無法到達，表示兩 SWTA 代數等價。
 * 怎麼算：
 *   1. 初始化最終向量：第一個 SWTA 的接收態對應 +1，第二個對應 -1。
 *   2. 初始化源狀態的向量空間及工作列表。
 *   3. BFS 遍歷親和程式：對每個狀態，將其行向量傳播至後繼狀態的行梯形矩陣。
 *   4. 對每個接收態，檢查其向量空間與最終向量相乘是否非零。
 *   5. 若任意接收態產生非零值，則返回 true（不等價）；否則返回 false（等價）。
 */
bool does_affine_program_reach_nonzero_final_states(const Affine_Program<Branch_Product_Sym>& program, const Underlying_SWTA_Info_Pair& swta_pair_info) {
    u64 swta_state_cnt = swta_pair_info.first_swta_info.state_cnt + swta_pair_info.second_swta_info.state_cnt;

    ACN_Matrix final_vector (swta_state_cnt, 1);
    for (State final_state : swta_pair_info.first_swta_info.final_states) {
        Algebraic_Complex_Number acn_one = Algebraic_Complex_Number::ONE();
        final_vector.set(final_state, 0, acn_one);
    }
    for (State final_state : swta_pair_info.second_swta_info.final_states) {
        Algebraic_Complex_Number acn_one = Algebraic_Complex_Number::ONE();
        final_vector.set(final_state + swta_pair_info.first_swta_info.state_cnt, 0, acn_one);
    }

    Propagation_Stats stats;
    Affine_Program_Propagation_Context context(program, final_vector, swta_pair_info);
    {  // Initialize the context
        context.state_space_dimension = swta_state_cnt;
        context.state_vector_spaces.reserve(program.number_of_states());

        for (State state = 0; state < program.number_of_states(); state++) {
            context.state_vector_spaces.push_back(ACN_Matrix(context.state_space_dimension, context.state_space_dimension));
        }
        context.state_deltas.resize(program.number_of_states());
    }

    { // Push the initial state with the initial vector into the worklist
        context.worklist.push_back(program.initial_state);
        context.state_deltas[program.initial_state].new_vector_count = 1;

        auto& initial_state_matrix = context.state_vector_spaces[program.initial_state];
        Algebraic_Complex_Number init_vector_value = Algebraic_Complex_Number::ONE();
        initial_state_matrix.set(0, swta_pair_info.first_swta_info.initial_states[0], init_vector_value);

        State second_initial_state_with_offset = swta_pair_info.second_swta_info.initial_states[0] + swta_pair_info.first_swta_info.state_cnt;
        initial_state_matrix.set(0, second_initial_state_with_offset, -init_vector_value);
    }

    while (!context.worklist.empty()) {
        stats.propagation_cnt += 1;
        // do_on_debug({
            // std::cout << "Processing to " <<  context.worklist.back() << " state matrix:\n " <<  context.state_vector_spaces[context.worklist.back()] << "\n";
        // });

        State current_state = context.worklist.back();
        State_Delta_Info old_delta_info = context.state_deltas[current_state];
        context.worklist.pop_back();

        // Mark in state deltas that we have already processe the new vectors
        context.state_deltas[current_state].new_vector_start += old_delta_info.new_vector_count;
        context.state_deltas[current_state].new_vector_count = 0;

        if (program.final_states.get_bit_value(current_state)) {
            bool contains_nonzeros = does_state_vector_space_contain_nonzeros(context, final_vector, current_state);
            if (contains_nonzeros) {
                goto terminate;
            }
        }

        for (u64 symbol = 0; symbol < program.number_of_symbols(); symbol++) {
            propagate_vector_from_state_matrix_to_successors(context, old_delta_info, current_state, symbol);
            if (context.final_state_with_nonzero > -1) goto terminate;
        }
    }

    for (State state = 0; state < program.number_of_states(); state++) {
        if (program.final_states.get_bit_value(state)) {
            bool contains_nonzeros = does_state_vector_space_contain_nonzeros(context, final_vector, state);
            if (contains_nonzeros) {
                goto terminate;
            }
        }
    }

    terminate:
        bool reaches_nonzero = context.final_state_with_nonzero >= 0;

        do_on_debug({std::cout << "Performed " << context.propagation_log.size() << " propagations.\n"; });

        if (reaches_nonzero && DEBUG_COLOR_EQUIVALENCE) {
            auto filtered_propagations = filter_propagations_for_those_that_affect_state(context.propagation_log, context.final_state_with_nonzero);
            dump_propagations(filtered_propagations, context);
        }

        return reaches_nonzero;
}

struct Colored_Product_State {
    State in_first;
    State in_second;
    State in_frontier;

    s64 handle = -1;

    bool operator<(const Colored_Product_State& other) const {
        if (this->in_first < other.in_first) return true;
        if (this->in_first > other.in_first) return false;

        if (this->in_second < other.in_second) return true;
        if (this->in_second > other.in_second) return false;

        return this->in_frontier < other.in_frontier;
    }
};

/*
 * compose_transition_matrices
 * 用途：
 *   將兩個轉移矩陣進行區塊對角合成，用於積親和程式的轉移矩陣建構。
 * 怎麼算：
 *   新矩陣為區塊對角形式：左上為 first_matrix，右下為 second_matrix，
 *   其餘位置為 0。
 */
ACN_Matrix compose_transition_matrices(const ACN_Matrix& first_matrix, const ACN_Matrix& second_matrix) {
    u64 width  = first_matrix.width + second_matrix.width;
    u64 height = first_matrix.height + second_matrix.height;
    ACN_Matrix result(height, width);

    for (u64 first_matrix_row_idx = 0; first_matrix_row_idx < first_matrix.height; first_matrix_row_idx++) {
        for (u64 first_matrix_col_idx = 0; first_matrix_col_idx < first_matrix.width; first_matrix_col_idx++) {
            auto& elem = first_matrix.at(first_matrix_row_idx, first_matrix_col_idx);
            result.set(first_matrix_row_idx, first_matrix_col_idx, elem);
        }
    }

    for (u64 second_matrix_row_idx = 0; second_matrix_row_idx < second_matrix.height; second_matrix_row_idx++) {
        for (u64 second_matrix_col_idx = 0; second_matrix_col_idx < second_matrix.width; second_matrix_col_idx++) {
            auto& elem = second_matrix.at(second_matrix_row_idx, second_matrix_col_idx);
            result.set(second_matrix_row_idx + first_matrix.height, second_matrix_col_idx + first_matrix.width, elem);
        }
    }

    return result;
}

/*
 * group_transition_symbols_by_color
 * 用途：
 *   將親和程式符號儲存中的符號按顏色分組，方便按顏色查找匹配符號。
 * 怎麼算：
 *   遍歷所有符號，根據其顏色欄位將指標放入相應顏色的組中。
 */
std::vector<std::vector<const Affine_Program<Branch_Selector>::Transition_Symbol*>>
group_transition_symbols_by_color(const Affine_Program<Branch_Selector>::Symbol_Store& symbol_store, u64 color_cnt) {
    std::vector<std::vector<const Affine_Program<Branch_Selector>::Transition_Symbol*>> result;
    result.resize(color_cnt);

    for (auto& symbol : symbol_store) {
        result[symbol.info.color].push_back(&symbol);
    }

    return result;
}

struct Product_Alphabet_Prep {
    std::vector<std::vector<const Affine_Program<Branch_Selector>::Transition_Symbol*>> first_transition_syms_by_color;
    std::vector<std::vector<const Affine_Program<Branch_Selector>::Transition_Symbol*>> second_transition_syms_by_color;
    Affine_Program<Branch_Product_Sym>::Symbol_Handles symbol_handles;
    Affine_Program<Branch_Product_Sym>::Symbol_Store   symbol_store;
};

/*
 * prepare_transition_symbols_for_product_automaton
 * 用途：
 *   為積親和程式準備符號字母表與轉移矩陣。配對兩個親和程式中同色且同方向的
 *   分支選擇符，合成新的積符號及其對應的區塊對角轉移矩陣。
 * 怎麼算：
 *   先按顏色分組兩個親和程式的符號。對每種顏色，配對左右方向相同的符號對，
 *   使用 compose_transition_matrices 合成其轉移矩陣。
 */
Product_Alphabet_Prep prepare_transition_symbols_for_product_automaton(const Affine_Program<Branch_Selector>& first_ap, const Affine_Program<Branch_Selector>& second_ap, u64 color_cnt) {
    auto selectors_in_first_by_color  = group_transition_symbols_by_color(first_ap.symbol_store, color_cnt);
    auto selectors_in_second_by_color = group_transition_symbols_by_color(second_ap.symbol_store, color_cnt);

    Affine_Program<Branch_Product_Sym>::Symbol_Handles handles;
    Affine_Program<Branch_Product_Sym>::Symbol_Store   store;

    for (Color color = 0; color < color_cnt; color++) {
        for (auto& first_ap_symbol : selectors_in_first_by_color[color]) {
            for (auto& second_ap_symbol : selectors_in_second_by_color[color]) {
                if (first_ap_symbol->info.tag != second_ap_symbol->info.tag) continue;

                Branch_Product_Sym sym = {
                    .color = color,
                    .first_sym = first_ap_symbol->info.symbol, .first_tag = first_ap_symbol->info.tag,
                    .second_sym = second_ap_symbol->info.symbol, .second_tag = second_ap_symbol->info.tag
                };

                handles.emplace(sym, handles.size());

                ACN_Matrix matrix = compose_transition_matrices(first_ap_symbol->matrix, second_ap_symbol->matrix);
                Affine_Program<Branch_Product_Sym>::Transition_Symbol symbol_to_store = {.info = sym, .matrix = matrix};
                store.push_back(symbol_to_store);
            }
        }
    }

    Product_Alphabet_Prep result(selectors_in_first_by_color, selectors_in_second_by_color, handles, store);
    return result;
}


struct Colored_Product_Build_Context {
    Worklist_Construction_Context<Colored_Product_State> worklist_info;
    const Affine_Program<Branch_Selector>& first_ap;
    const Affine_Program<Branch_Selector>& second_ap;
    const NFA& frontier;
    Affine_Program_Builder<Branch_Product_Sym> builder;

    Colored_Product_Build_Context(
        const Affine_Program<Branch_Selector>& first_ap,
        const Affine_Program<Branch_Selector>& second_ap,
        const NFA& frontier,
        const Affine_Program<Branch_Product_Sym>::Symbol_Handles& handles,
        const Affine_Program<Branch_Product_Sym>::Symbol_Store& store) :
        first_ap(first_ap),
        second_ap(second_ap),
        frontier(frontier),
        builder(handles, store) {}
};

/*
 * product_build_step_along_color
 * 用途：
 *   給定積狀態及顏色，在兩個親和程式中查找相應轉移，並在積自動機中添加所有
 *   可能的後繼狀態。
 * 怎麼算：
 *   取出兩親和程式中該顏色且相同方向的分支符號對，查其轉移，產生所有狀態對組合。
 */
void product_build_step_along_color(Colored_Product_Build_Context& context, Product_Alphabet_Prep& alphabet_prep, const Colored_Product_State* current_state, Color color, State frontier_nfa_state) {
    for (auto first_branch_selector_ptr : alphabet_prep.first_transition_syms_by_color[color]) {
        for (auto second_branch_selector_ptr : alphabet_prep.second_transition_syms_by_color[color]) {
            if (first_branch_selector_ptr->info.tag != second_branch_selector_ptr->info.tag) {
                // We have to go along the same branches in the tree, otherwise we would be comparing leaves from completely different paths
                continue;
            }

            Branch_Product_Sym product_symbol = {
                .color = color,
                .first_sym = first_branch_selector_ptr->info.symbol, .first_tag = first_branch_selector_ptr->info.tag,
                .second_sym = second_branch_selector_ptr->info.symbol, .second_tag = second_branch_selector_ptr->info.tag
            };
            u64 product_symbol_handle = context.builder.symbol_handles.at(product_symbol);

            u64 first_branch_selector_handle  = context.first_ap.symbol_handles.at(first_branch_selector_ptr->info);
            u64 second_branch_selector_handle = context.second_ap.symbol_handles.at(second_branch_selector_ptr->info);

            auto& first_ap_successors  = context.first_ap.transition_fn[current_state->in_first][first_branch_selector_handle];
            auto& second_ap_successors = context.second_ap.transition_fn[current_state->in_second][second_branch_selector_handle];

            for (auto first_successor : first_ap_successors) {
                for (auto second_successor : second_ap_successors) {
                    Colored_Product_State imm_product_state = {.in_first = first_successor, .in_second = second_successor, .in_frontier = frontier_nfa_state};
                    auto dest_state = context.worklist_info.mark_discovery(imm_product_state);

                    context.builder.add_transition(current_state->handle, product_symbol_handle, dest_state->handle);
                }
            }
        }
    }
}

/*
 * build_colored_product_of_affine_programs
 * 用途：
 *   構建兩個親和程式在邊界自動機同步下的積自動機。此積用於檢查兩個 SWTA 是否
 *   代數等價。
 * 怎麼算：
 *   狀態為三元組 (第一親和狀態, 第二親和狀態, 邊界NFA狀態)。自三重初始狀態開始進行
 *   工作列表BFS。對每個三元組及每種顏色，在邊界NFA中進行轉移，同時在兩親和程式中
 *   查找對應轉移。接收態則三者皆為接收態。
 */
Affine_Program<Branch_Product_Sym>
build_colored_product_of_affine_programs(const Affine_Program<Branch_Selector>& first_ap, const Affine_Program<Branch_Selector>& second_ap, const NFA& frontier) {
    u64 color_cnt = frontier.alphabet_size();
    auto alphabet_prep = prepare_transition_symbols_for_product_automaton(first_ap, second_ap, color_cnt);

    Colored_Product_Build_Context context (first_ap, second_ap, frontier, alphabet_prep.symbol_handles, alphabet_prep.symbol_store);
    {
        Colored_Product_State imm_initial_state {.in_first = first_ap.initial_state, .in_second = second_ap.initial_state, .in_frontier = frontier.initial_states[0]};
        auto initial_state = context.worklist_info.mark_discovery(imm_initial_state);
        context.builder.mark_state_initial(initial_state->handle);
    }

    while (!context.worklist_info.worklist.empty()) {
        auto current_state = context.worklist_info.extract();

        bool first_accepts    = context.first_ap.final_states.get_bit_value(current_state->in_first);
        bool second_accepts   = context.second_ap.final_states.get_bit_value(current_state->in_second);
        bool frontier_accepts = context.frontier.final_states.get_bit_value(current_state->in_frontier);

        if (first_accepts && second_accepts && frontier_accepts) {
            context.builder.mark_state_final(current_state->handle);
        };

        for (Color color = 0; color < color_cnt; color++) {
            State nfa_post_state = context.frontier.transitions[current_state->in_frontier][color].at(0);
            product_build_step_along_color(context, alphabet_prep, current_state, color, nfa_post_state);
        }
    }

    auto result = context.builder.build(context.worklist_info.handles.size());

    do_on_debug({
        result.debug_data = new Affine_Program<Branch_Product_Sym>::Debug_Data;

        auto& state_names = result.debug_data->state_names;
        for (auto& [product_state, handle] : context.worklist_info.handles) {
            std::string& first_ap_state_name  = first_ap.debug_data->state_names[product_state.in_first];
            std::string& second_ap_state_name = second_ap.debug_data->state_names[product_state.in_second];
            std::string& frontier_state_name  = frontier.debug_data->state_names[product_state.in_frontier];
            state_names[product_state.handle] = "(" + first_ap_state_name + ", " + second_ap_state_name + ", " + frontier_state_name + ")";
        }
    });

    return result;
}


/*
 * operator<< for Branch_Product_Sym
 * 用途：
 *   輸出積符號 (Branch_Product_Sym)，以緊湊格式表示兩個分支選擇符的組合。
 * 怎麼算：
 *   輸出顏色、兩個符號及其方向標記。
 */
std::ostream& operator<<(std::ostream& stream, const Branch_Product_Sym& sym) {
    const char* first_sym_str  = sym.first_tag == Subtree_Tag::LEFT ? "L" : "R";
    const char* second_sym_str = sym.second_tag == Subtree_Tag::LEFT ? "L" : "R";

    stream << "<c=" << sym.color << "-" << sym.first_sym << first_sym_str << "-" << sym.second_sym << second_sym_str << ">";

    return stream;
}

struct Color_Symbol_Build_Context {
    const SWTA& swta;
    Worklist_Construction_Context<Macrostate> worklist_state;
    std::map<Color_Symbol, u64> symbol_handles;
    std::map<State, std::vector<std::vector<State>> > result_transitions;
    Bit_Set final_states;
    u64 number_of_symbols;

    Color_Symbol_Build_Context(const SWTA& swta_ref) :
        swta(swta_ref),
        final_states(0),
        number_of_symbols(swta.number_of_internal_symbols() * swta.number_of_colors()) {}
};

/*
 * compute_post_in_color_sym_abstraction
 * 用途：
 *   在顏色-符號抽象 (Color_Symbol_Abstraction) 中計算下一宏態。此抽象將邊界
 *   自動機進一步細化，字母表為 (顏色, 內部符號) 對。
 * 怎麼算：
 *   對源宏態中的所有狀態查詢 (顏色, 符號) 轉移，若任何狀態缺失該轉移則返回。
 *   否則將所有轉移的左右形式的目標狀態並集作為下一宏態。
 */
void compute_post_in_color_sym_abstraction(Color_Symbol_Build_Context& context, const Macrostate& source, Color_Symbol& color_symbol) {
    auto [insert_pos, was_inserted] = context.symbol_handles.emplace(color_symbol, context.symbol_handles.size());
    u64 symbol_handle = insert_pos->second;

    for (State state : source.state_names) {
        auto& transition = context.swta.transitions[state][color_symbol.symbol][color_symbol.color];
        if (!transition.is_present()) return;  // All of the members need to define a transition
    }

    Macrostate imm_post (context.swta.number_of_states());

    for (State state : source.state_names) {
        auto& transition = context.swta.transitions[state][color_symbol.symbol][color_symbol.color];

        for (auto& component : transition.left.components) {
            imm_post.state_set.set_bit(component.state);
        }

        for (auto& component : transition.right.components) {
            imm_post.state_set.set_bit(component.state);
        }
    }

    auto post = context.worklist_state.mark_discovery(imm_post);
    {
        auto mut_post_ptr = const_cast<Macrostate*>(post);
        if (mut_post_ptr->state_names.empty()) {
            mut_post_ptr->populate_state_names_from_set();
        }
    }

    auto& transitions_from_state = context.result_transitions[source.handle];
    if (transitions_from_state.empty()) {
        transitions_from_state.resize(context.number_of_symbols);
    }

    transitions_from_state [symbol_handle].push_back(post->handle);
}

/*
 * build_color_internal_symbol_abstraction
 * 用途：
 *   為 SWTA 構建顏色-符號抽象：一個 NFA，狀態為宏態，字母表為 (顏色, 內部符號) 對。
 *   此抽象用於親和程式的建構。
 * 怎麼算：
 *   自初始宏態開始進行工作列表BFS。對每個宏態及每個 (顏色, 符號) 組合，計算下一宏態
 *   並添加轉移。若宏態內所有狀態都是葉片態則標記為接收態。
 */
Color_Symbol_Abstraction build_color_internal_symbol_abstraction(const SWTA& swta) {
    Color_Symbol_Build_Context context(swta);

    {
        Macrostate initial_macrostate (context.swta.number_of_states(), swta.initial_states);
        context.worklist_state.mark_discovery(initial_macrostate);
    }

    while (context.worklist_state.has_more_to_explore()) {
        auto current_macrostate = context.worklist_state.extract();

        if (context.swta.states_with_leaf_transitions.is_superset(current_macrostate->state_set)) {
            context.final_states.grow_and_set_bit(current_macrostate->handle);
        }

        for (Internal_Symbol internal_sym = 0; internal_sym < context.swta.number_of_internal_symbols(); internal_sym++) {
            for (Color color = 0; color < context.swta.number_of_colors(); color++) {
                Color_Symbol color_symbol = {.color = color, .symbol = internal_sym};
                compute_post_in_color_sym_abstraction(context, *current_macrostate, color_symbol);
            }
        }
    }

    context.final_states.grow(context.worklist_state.handles.size());

    NFA::Transition_Fn transitions;
    transitions.resize(context.worklist_state.handles.size());

    for (State state = 0; state < context.worklist_state.handles.size(); state++) {
        auto& created_transitions = context.result_transitions[state];
        if (!created_transitions.empty()) {
            transitions[state] = std::move(created_transitions);
        } else {
            transitions[state].resize(context.number_of_symbols);
        }
    }

    NFA resulting_abstraction({0}, context.final_states, transitions);

    do_on_debug({
        resulting_abstraction.debug_data = new NFA::Debug_Data;

        for (auto& [macrostate, handle] : context.worklist_state.handles) {
            resulting_abstraction.debug_data->state_names[handle] = macrostate.to_string();
        }
    });

    Color_Symbol_Abstraction result = {
        .abstraction = resulting_abstraction,
        .symbol_handles = context.symbol_handles,
    };

    return result;
}

/*
 * write_wtt_transition_with_debug_data
 * 用途：
 *   以人類可讀格式列印 WTT 轉移，使用 debug_data 中的狀態名稱。
 * 怎麼算：
 *   分別輸出左右子樹部份，使用提供的狀態名稱映射；若係數為整數則簡化輸出。
 */
void write_wtt_transition_with_debug_data(std::ostream& os, const WTT::Transition& transition, const WTT::Debug_Data* debug_data) {
    auto write_norm_with_debug_info = [&debug_data](std::ostream& target, const Linear_Form& form, const char* subtree_info, bool needs_leading_plus) {
        if (needs_leading_plus) {
            target << " + ";
        }

        for (u64 i = 0; i < form.size(); i++) {
            auto& component = form.components[i];
            std::string state_name = debug_data->state_names.contains(component.state) ? debug_data->state_names.at(component.state) : "q" + std::to_string(component.state);

            if (component.coef.is_integer()) {
                s64 value = mpz_get_si(component.coef.a);
                target << value << "*" << state_name << subtree_info;
            } else {
                target << component.coef << "*" << state_name << subtree_info;
            }

            if (i < form.size() - 1) target << " + ";
        }
    };


    os << "LEFT SUBTREE: ";
    bool needs_plus = false; // True, if anything has been written to the output stream
    if (!transition.ll.empty()) {
        write_norm_with_debug_info(os, transition.ll, "(L)", needs_plus);
        needs_plus = true;
    }

    if (!transition.lr.empty()) {
        write_norm_with_debug_info(os, transition.lr, "(R)", needs_plus);
        needs_plus = true;
    }

    if (transition.ll.empty() && transition.lr.empty()) {
        os << "0";
    }

    os << "; RIGHT SUBTREE: ";
    needs_plus = false; // Reset

    if (!transition.rl.empty()) {
        write_norm_with_debug_info(os, transition.rl, "(L)", needs_plus);
        needs_plus = true;
    }

    if (!transition.rr.empty()) {
        write_norm_with_debug_info(os, transition.rr, "(R)", needs_plus);
        needs_plus = true;
    }

    if (transition.rl.empty() && transition.rr.empty()) {
        os << "0";
    }
}

/*
 * write_wtt_with_debug_data
 * 用途：
 *   使用調試資料以可讀格式輸出整個 WTT 的結構。
 * 怎麼算：
 *   列印初始態、葉片態與所有轉移，使用 write_wtt_transition_with_debug_data
 *   格式化各轉移。
 */
void write_wtt_with_debug_data(std::ostream& os, const WTT& wtt) {
    os << "WTT {\n";
    os << "  initial states: " << wtt.initial_states << "\n";
    os << "  leaf states: " << wtt.states_with_leaf_transitions.into_vector() << "\n";

    u64 transition_cnt = 0;
    for (u64 state = 0; state < wtt.number_of_states(); state++) {
        std::string state_name = wtt.debug_data->state_names.contains(state) ? wtt.debug_data->state_names.at(state) : "q" + std::to_string(state);

        const std::vector<WTT::Transition>& transitions_from_state = wtt.transitions[state];
        for (Internal_Symbol sym = 0; sym < wtt.number_of_internal_symbols(); sym++) {
            auto& transition = transitions_from_state[sym];
            if (!transition.is_present()) continue;
            os << "  " << state_name << "--(sym=" << sym << ")-->: ";
            write_wtt_transition_with_debug_data(os, transitions_from_state[sym], wtt.debug_data);
            os << "\n";
            transition_cnt += 1;
        }
    }

    std::cout << "  ... " << transition_cnt << " transitions total.\n";

    os << "}";
}

/*
 * compose_wtts_horizontally
 * 用途：
 *   水平合成兩個 WTT：將第一個 WTT 的葉片狀態重新連接至第二個 WTT 的初始狀態，
 *   以實現第一個完成後繼續第二個的效果。
 * 怎麼算：
 *   複製第一個 WTT 的結構，將其所有葉片轉移重定向至第二個 WTT 的初始狀態（加上狀態偏移），
 *   然後追加第二個 WTT 的所有狀態與轉移（狀態編號偏移）。新的葉片態為第二個 WTT 的葉片態。
 */
WTT compose_wtts_horizontally(const WTT& first, const WTT& second) {
    WTT result(first);

    Bit_Set original_leaves (first.states_with_leaf_transitions);

    auto redirect_leaves_to_initial_state = [&first, &second](Linear_Form& form, u64 state_offset) {
        for (auto& component : form.components) {
            if (first.states_with_leaf_transitions.get_bit_value(component.state)) {
                component.state = second.initial_states[0] + state_offset;
            }
        }
    };

    u64 state_offset = first.number_of_states();

    auto offset_states_in_form = [](Linear_Form& form, u64 offset) {
        for (auto& component : form.components) {
            component.state += offset;
        }
    };

    for (auto& transitions_from_state : result.transitions) {
        for (auto& transition : transitions_from_state) {
            redirect_leaves_to_initial_state(transition.ll, state_offset);
            redirect_leaves_to_initial_state(transition.lr, state_offset);
            redirect_leaves_to_initial_state(transition.rl, state_offset);
            redirect_leaves_to_initial_state(transition.rr, state_offset);
        }
    }

    result.states_with_leaf_transitions.clear();
    result.states_with_leaf_transitions.grow(first.number_of_states() + second.number_of_states());

    for (State second_state = 0; second_state < second.number_of_states(); second_state++) {
        result.transitions.push_back({});
        result.transitions[state_offset + second_state].resize(second.number_of_internal_symbols());

        for (Internal_Symbol internal_symbol = 0; internal_symbol < second.number_of_internal_symbols(); internal_symbol++) {
            auto transition = second.transitions[second_state][internal_symbol];

            offset_states_in_form(transition.ll, state_offset);
            offset_states_in_form(transition.lr, state_offset);
            offset_states_in_form(transition.rl, state_offset);
            offset_states_in_form(transition.rr, state_offset);

            result.transitions[second_state + state_offset][internal_symbol] = transition;
        }

        if (second.states_with_leaf_transitions.get_bit_value(second_state)) {
            result.states_with_leaf_transitions.set_bit(second_state + state_offset);
        }
    }

    return result;
}

/*
 * evaluate_wtt_on_tree
 * 用途：
 *   遞迴地在一棵具體二元樹上評估 WTT，計算輸出值。樹以平坦向量形式編碼，
 *   其中葉值在向量中。
 * 怎麼算：
 *   若當前狀態是葉片態，直接返回樹的葉值。否則，遞迴評估左右子樹，根據轉移的
 *   四個方向形式（ll、lr、rl、rr）線性組合子樹的評估結果。
 */
std::vector<Algebraic_Complex_Number> evaluate_wtt_on_tree(const WTT& wtt, State root, const std::vector<Algebraic_Complex_Number>& tree, const std::vector<u32>& internal_symbols, u32 sym_idx) {
    if (wtt.states_with_leaf_transitions.get_bit_value(root)) {
        assert (tree.size() == 1);
        return { tree[0] };
    }

    std::vector<Algebraic_Complex_Number> left_half, right_half;
    u64 i = 0;

    for (; i < tree.size() / 2; i++) {
        left_half.push_back(tree[i]);
    }

    for (; i < tree.size(); i++) {
        right_half.push_back(tree[i]);
    }

    auto symbol = internal_symbols[sym_idx];
    auto& transition = wtt.transitions[root][symbol];

    assert (tree.size() >= 2);


    std::vector<Algebraic_Complex_Number> left_result, right_result;
    left_result.resize(tree.size() / 2);
    right_result.resize(tree.size() / 2);

    auto vector_mad = [](std::vector<Algebraic_Complex_Number>& dest, std::vector<Algebraic_Complex_Number>& vector, const Algebraic_Complex_Number& constant) {
        for (u64 idx = 0; idx < dest.size(); idx++) {
            dest[idx] += constant * vector[idx];
        }
    };

    assert (transition.is_present());

    for (auto& component : transition.ll.components) {
        auto component_result = evaluate_wtt_on_tree(wtt, component.state, left_half, internal_symbols, sym_idx+1);
        vector_mad(left_result, component_result, component.coef);
    }

    for (auto& component : transition.lr.components) {
        auto component_result = evaluate_wtt_on_tree(wtt, component.state, right_half, internal_symbols, sym_idx+1);
        vector_mad(left_result, component_result, component.coef);
    }

    for (auto& component : transition.rl.components) {
        auto component_result = evaluate_wtt_on_tree(wtt, component.state, left_half, internal_symbols, sym_idx+1);
        vector_mad(right_result, component_result, component.coef);
    }

    for (auto& component : transition.rr.components) {
        auto component_result = evaluate_wtt_on_tree(wtt, component.state, right_half, internal_symbols, sym_idx+1);
        vector_mad(right_result, component_result, component.coef);
    }

    std::vector<Algebraic_Complex_Number> result (left_result.begin(), left_result.end());
    for (auto& acn : right_result) result.push_back(acn);

    return result;
}
