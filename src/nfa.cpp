#include "nfa.hpp"
#include "basics.hpp"

#include <map>
#include <algorithm>
#include <sstream>


/*
 * operator<< (Macrostate)
 * 用途：
 *   將巨狀態的內容（包含狀態名稱與 handle）以人類可讀形式輸出。
 * 怎麼算：
 *   將 state_names 與 handle 依序寫入串流，格式為 "Macrostate [..] (handle=..)"。
 */
std::ostream& operator<<(std::ostream& os, const Macrostate& macrostate) {
    os << "Macrostate " << macrostate.state_names << " (handle=" << macrostate.handle << ")";
    return os;
}


/*
 * compute_post
 * 用途：
 *   計算給定巨狀態在指定符號下的後繼巨狀態。
 * 怎麼算：
 *   遍歷巨狀態內所有原子狀態，依符號收集轉移後的狀態到 Bit_Set 並同步填入 state_names，最後回傳新巨狀態。
 */
Macrostate compute_post(const Macrostate* macrostate, const NFA& nfa, u64 symbol) {
    Macrostate post (nfa.number_of_states());

    for (NFA::State state : macrostate->state_names) {
        auto& transitions_from_state = nfa.transitions[state];
        auto& transitions_for_symbol = transitions_from_state[symbol];

        for (auto post_state : transitions_for_symbol) {
            post.state_set.set_bit(post_state);
        }
    }

    for (NFA::State state = 0; state < nfa.number_of_states(); state++) {
        if (post.state_set.get_bit_value(state)) {
            post.state_names.push_back(state);
        }
    }

    return post;
}

/*
 * NFA::complete
 * 用途：
 *   使 NFA 完備，確保每個狀態在每個符號下皆有轉移。
 * 怎麼算：
 *   若缺少轉移則新增指向沉沒狀態的邊，必要時建立沉沒狀態並自循環，回傳是否有新增沉沒狀態。
 */
bool NFA::complete() {
    u64 state_cnt = this->number_of_states();
    State sink_state = state_cnt;

    bool was_sink_state_needed = false;
    for (State state = 0; state < state_cnt; state++) {
        for (u64 color = 0; color < this->alphabet_size(); color++) {
            if (this->transitions[state][color].empty()) {
                this->transitions[state][color].push_back(sink_state);
                was_sink_state_needed = true;
            }
        }
    }

    if (was_sink_state_needed) {
        std::vector<std::vector<State>> sink_state_transitions;
        sink_state_transitions.resize(this->alphabet_size());

        for (u64 color = 0; color < this->alphabet_size(); color++) {
            sink_state_transitions[color].push_back(sink_state);
        }

        this->transitions.push_back(sink_state_transitions);
    }

    return was_sink_state_needed;
}

/*
 * NFA::determinize
 * 用途：
 *   將 NFA 轉換為等價的 DFA。
 * 怎麼算：
 *   以子集構造法遍歷巨狀態：從初始巨狀態開始，對每個符號計算 post，為新巨狀態分配 handle 並建立轉移，最終重建 DFA 結構。
 */
NFA NFA::determinize() const {
    std::map<Macrostate, NFA::State> handles;

    std::vector<const Macrostate*> worklist;

    { // Initialize frontier
        Macrostate initial_macrostate (this->number_of_states(), initial_states);
        initial_macrostate.handle = 0;
        auto [insert_pos, was_inserted] = handles.emplace(initial_macrostate, 0);
        worklist.push_back(&insert_pos->first);
    }

    std::map<State, std::vector<std::vector<State>>> resulting_transitions; // Use an associative container here, because we do now know the final number of states
    Bit_Set final_macrostates (0);

    while (!worklist.empty()) {
        auto macrostate = worklist.back();
        worklist.pop_back();

        if (!this->final_states.is_intersection_empty(macrostate->state_set)) {
            // There is at least one state that can make a leaf transition
            final_macrostates.grow_and_set_bit(macrostate->handle);
        }

        for (u64 symbol = 0; symbol < this->alphabet_size(); symbol++) {
            Macrostate post = compute_post(macrostate, *this, symbol);

            if (post.empty()) {
                continue;
            }

            post.handle = handles.size(); // Set this speculatively, before we store it in the handles map

            const auto& [insert_pos, was_inserted] = handles.emplace(post, handles.size());
            if (was_inserted) { // the macrostate is new, we need to explore it
                worklist.push_back(&insert_pos->first);
            } else { // There already is such a macrostate, so our speculation with the handle was incorrect
                post.handle = insert_pos->second;
            }

            auto& transitions_from_this_macrostate = resulting_transitions[macrostate->handle];
            if (transitions_from_this_macrostate.size() < this->alphabet_size()) {
                transitions_from_this_macrostate.resize(this->alphabet_size());
            }

            transitions_from_this_macrostate[symbol].push_back(post.handle);
        }
    }

    std::vector<NFA::Transitions_From_State> ordered_resulting_transitions;
    ordered_resulting_transitions.resize(handles.size());

    for (auto& [state, transitions_from_state] : resulting_transitions) {
        ordered_resulting_transitions[state] = transitions_from_state;
    }

    for (State state = 0; state < handles.size(); state++) {
        if (ordered_resulting_transitions[state].empty()) {
            ordered_resulting_transitions[state].resize(this->alphabet_size());
        }
    }

    NFA result ({0}, final_macrostates, ordered_resulting_transitions);

    do_on_debug({
        result.debug_data = new NFA::Debug_Data;

        for (auto& [macrostate, handle] : handles) {
            std::stringstream string_stream;
            string_stream << "{";

            u64 written_states_cnt = 0;
            for (State state : macrostate.state_names) {
                std::string state_name = (this->debug_data->state_names.contains(state)) ?  this->debug_data->state_names.at(state) : std::to_string(state);
                string_stream << state_name;

                if (written_states_cnt + 1 != macrostate.state_names.size()) {
                    string_stream << ", ";
                }

                written_states_cnt += 1;
            }

            string_stream << "}";

            std::string macrostate_name = string_stream.str();
            result.debug_data->state_names[handle] = macrostate_name;
        }
    });

    return result;
}

/*
 * NFA::is_every_state_accepting
 * 用途：
 *   檢查自動機是否所有狀態皆為接受狀態。
 * 怎麼算：
 *   呼叫 final_states.are_all_bits_set() 確認所有位元皆為 1。
 */
bool NFA::is_every_state_accepting() const {
    return this->final_states.are_all_bits_set();
}

/*
 * NFA::write_dot
 * 用途：
 *   將 NFA 輸出為 Graphviz DOT 描述以供可視化。
 * 怎麼算：
 *   先寫入初始箭頭與節點樣式（區分終態與一般狀態），再列舉所有轉移並以 label 標示符號，最後關閉圖形定義。
 */
void NFA::write_dot(std::ostream& stream) const {
    stream << "digraph NFA {\n";

    for (State initial_state : this->initial_states) {
        stream << "  qInit" << initial_state << " [shape=none, label=\"\"]\n";
    }
    for (State state = 0; state < this->number_of_states(); state++) {
        stream << "  q" << state;
        const char* color = this->final_states.get_bit_value(state) ? "green" : "black";
        if (this->debug_data != nullptr && this->debug_data->state_names.contains(state)) {
            stream << " ["
                   << "label=\"" << this->debug_data->state_names.at(state) << "\", "
                   << "color=\"" << color << "\""
                   << "]";
        } else {
            stream << " ["
                   << "label=\"" << state << "\", "
                   << "color=\"" << color << "\""
                   << "]";
        }
        stream << "\n";
    }

    for (State initial_state : this->initial_states) {
        stream << "  qInit" << initial_state << " -> q" << initial_state << "\n";
    }

    for (State state = 0; state < number_of_states(); state++) {
        for (u64 symbol = 0; symbol < transitions[state].size(); symbol++) {
            for (State destination : transitions[state][symbol]) {
                stream << "  q" << state << " -> q" << destination << " [label=\"" << symbol << "\"]\n";
            }
        }
    }

    stream << "}\n";
}

/*
 * print_diseq_witness
 * 用途：
 *   在等價性檢查失敗時，輸出導致不等價的顏色序列。
 * 怎麼算：
 *   依 pred 指標回溯狀態配對，重建顏色字串並逆序輸出，展示抵達分歧的輸入字。
 */
void print_diseq_witness(Worklist_Construction_Context<State_Pair>& ctx, u64 state) {
    std::vector<State_Pair> states_by_handles;
    states_by_handles.resize(ctx.handles.size());

    for (auto& [state, handle] : ctx.handles) {
        states_by_handles[handle] = state;
    }

    std::vector<u64> color_word;

    State_Pair& current_state = states_by_handles[state];

    while (current_state.pred != -1) {
        u32 color       = current_state.pred >> 32;
        u32 pred_handle = current_state.pred & (0xFFFFFFFF);

        color_word.push_back(color);

        current_state = states_by_handles[pred_handle];
    }

    std::cout << "Diseq witness: ";
    for (auto color_it = color_word.rbegin(); color_it != color_word.rend(); color_it++) {
        std::cout << *color_it << ", ";
    }
    std::cout << "\n";
}

/*
 * are_two_complete_dfas_equivalent
 * 用途：
 *   檢查兩個已完備 DFA 是否接受同一語言。
 * 怎麼算：
 *   以工作清單同步遍歷狀態配對，若某配對終態性不同則回傳 false；否則依每個符號擴展配對直到遍歷完畢，最後回傳 true。
 */
bool are_two_complete_dfas_equivalent(const NFA& first_nfa, NFA& second_nfa) {
    u64 alphabet_size = first_nfa.alphabet_size();
    assert(alphabet_size == second_nfa.alphabet_size());

    Worklist_Construction_Context<State_Pair> context;

    {
        State_Pair initial_pair = {.first = first_nfa.initial_states.at(0), .second = second_nfa.initial_states.at(0), .handle = 0};
        context.mark_discovery(initial_pair);
    }

    while (!context.worklist.empty()) {
        auto current_pair = context.extract();

        bool is_final_in_first  = first_nfa.final_states.get_bit_value(current_pair->first);
        bool is_final_in_second = second_nfa.final_states.get_bit_value(current_pair->second);

        if (is_final_in_first != is_final_in_second) {
            if (DEBUG) {
                std::cout << "First  accepts? " << is_final_in_first << "\n";
                std::cout << "Second accepts? " << is_final_in_second << "\n";
                print_diseq_witness(context, current_pair->handle);
            }
            return false;
        }

        for (u32 color = 0; color < alphabet_size; color++) {
            auto first_post  = first_nfa.transitions[current_pair->first][color][0];
            auto second_post = second_nfa.transitions[current_pair->second][color][0];

            u64 pred = color << 32 | current_pair->handle;

            State_Pair discovery = {.first = first_post, .second = second_post, .pred = pred };
            context.mark_discovery(discovery);
        }
    }

    return true;
}

/*
 * operator<< (State_Pair)
 * 用途：
 *   以易讀格式輸出狀態配對與 handle。
 * 怎麼算：
 *   依序將 first、second 與 handle 填入串流，格式為 "(first, second, handle=...)"。
 */
std::ostream& operator<<(std::ostream& os, const State_Pair& state) {
    os << "(" << state.first << ", " << state.second << ", handle=" << state.handle << ")";
    return os;
}


/*
 * NFA_Builder::build
 * 用途：
 *   根據累積的轉移與初終態資訊建構 NFA 實例。
 * 怎麼算：
 *   推導所需狀態數（必要時從轉移中求最大狀態），為每個狀態配置符號轉移向量並填入集合轉移內容，最後擴充終態位元集合並建立 NFA。
 */
NFA NFA_Builder::build(s64 state_cnt) {
    std::vector<NFA::Transitions_From_State> result_transitions;

    if (state_cnt < 0) {
        for (auto& [state, outgoing_transitions] : this->transitions) {
            for (auto& post_vector : outgoing_transitions) {
                NFA::State post_vector_max = *std::max_element(post_vector.begin(), post_vector.end());
                state_cnt = std::max(state_cnt, static_cast<s64>(post_vector_max));
            }
        }
    }

    result_transitions.resize(state_cnt);

    for (auto& [state, transitions_from_state] : this->transitions) {
        result_transitions[state].resize(this->alphabet_size);

        for (u64 symbol = 0; symbol < this->alphabet_size; symbol++) {
            result_transitions[state][symbol] = std::vector<NFA::State>(transitions_from_state[symbol].begin(), transitions_from_state[symbol].end());
        }
    }

    for (NFA::State state = 0; state < state_cnt; state++) {
        if (result_transitions[state].empty()) {
            result_transitions[state].resize(this->alphabet_size);
        }
    }

    final_states.grow(state_cnt);

    NFA result ({0}, final_states, result_transitions);
    return result;
}

/*
 * NFA_Builder::add_transition
 * 用途：
 *   在建構期間為指定狀態與符號新增轉移目標。
 * 怎麼算：
 *   確保來源狀態的轉移向量尺寸足夠，然後將 destination 插入對應符號的集合以去除重複。
 */
void NFA_Builder::add_transition(NFA::State source_state, u64 symbol, NFA::State destination) {
    auto& transitions_from_state = this->transitions[source_state];
    if (transitions_from_state.empty()) {
        transitions_from_state.resize(this->alphabet_size);
    }

    transitions_from_state[symbol].insert(destination);
}
