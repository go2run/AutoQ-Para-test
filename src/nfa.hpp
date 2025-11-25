#pragma once

#include "arith.hpp"
#include "bit_set.hpp"

#include <map>
#include <set>
#include <string>
#include <sstream>

namespace std {
    template <typename T>
    ostream& operator<<(ostream& os, std::vector<T> vec) {
        os << "[";
        for (u64 idx = 0; idx < vec.size(); idx++) {
            auto& elem = vec[idx];
            os << elem;
            if (idx + 1 != vec.size()) {
                os << ", ";
            }
        }
        os << "]";
        return os;
    }
}

struct NFA {
    using State = u64;
    using Transitions_From_State = std::vector<std::vector<State>>;
    using Transition_Fn = std::vector<Transitions_From_State>;

    struct Debug_Data {
        std::map<State, std::string> state_names;
    };

    std::vector<State> initial_states;
    Bit_Set            final_states;
    Transition_Fn      transitions;
    Debug_Data*        debug_data = nullptr;

    /*
     * NFA()
     * 用途：
     *   建立空的自動機實例，預設無狀態與轉移。
     * 怎麼算：
     *   初始化初始狀態空向量、零長度終態 Bit_Set 與空轉移表。
     */
    NFA() : initial_states({}), final_states(0), transitions({}) {}
    /*
     * NFA(const std::vector<State>&, const std::vector<State>&, const std::vector<Transitions_From_State>&)
     * 用途：
     *   以初始狀態集合、終態集合與轉移表建立 NFA。
     * 怎麼算：
     *   依提供的狀態數設定 final_states 大小並將指定終態位置設為 1，同步儲存轉移表。
     */
    NFA(const std::vector<State>& init_states, const std::vector<State>& fin_states, const std::vector<Transitions_From_State>& aut_transitions) :
        initial_states(init_states),
        final_states(aut_transitions.size()),
        transitions(aut_transitions)
    {
        for (auto state: fin_states) {
            final_states.set_bit(state, true);
        }
    }

    /*
     * NFA(const NFA& other)
     * 用途：
     *   拷貝建構函式，複製另一個 NFA 的結構與除錯資料。
     * 怎麼算：
     *   直接複製初始狀態、終態 Bit_Set 與轉移表；若存在 debug_data 則配置新物件並深拷貝。
     */
    NFA(const NFA& other) :
        initial_states(other.initial_states),
        final_states(other.final_states),
        transitions(other.transitions)
    {
        if (other.debug_data != nullptr) {
            this->debug_data = new Debug_Data(*other.debug_data);
        }
    };

    /*
     * NFA(const std::vector<State>&, const Bit_Set&, const std::vector<Transitions_From_State>&)
     * 用途：
     *   以現成的終態位元集合建立 NFA，避免重新標記終態。
     * 怎麼算：
     *   直接指定初始狀態、終態 Bit_Set 及轉移表。
     */
    NFA(const std::vector<State>& init_states, const Bit_Set& fin_states, const std::vector<Transitions_From_State>& aut_transitions) :
        initial_states(init_states),
        final_states(fin_states),
        transitions(aut_transitions) {}

    /*
     * ~NFA
     * 用途：
     *   解構並釋放可選的除錯資料。
     * 怎麼算：
     *   若 debug_data 非空則 delete，並將指標設為 nullptr。
     */
    ~NFA() {
        if (this->debug_data != nullptr) delete this->debug_data;
        this->debug_data = nullptr;
    }

    /*
     * number_of_states
     * 用途：
     *   取得自動機狀態數量。
     * 怎麼算：
     *   回傳轉移表的外層向量長度。
     */
    u64 number_of_states() const {
        return this->transitions.size();
    }

    /*
     * alphabet_size
     * 用途：
     *   取得字母表大小（每個狀態的符號轉移數）。
     * 怎麼算：
     *   回傳第一個狀態的轉移向量長度。
     */
    u64 alphabet_size() const {
        return this->transitions[0].size();
    }

    NFA determinize() const;
    bool complete();
    bool is_every_state_accepting() const;
    void write_dot(std::ostream& stream) const;
};


struct NFA_Builder {
    std::map<NFA::State, std::vector<std::set<NFA::State>>> transitions;
    std::vector<NFA::State> initial_states;
    Bit_Set final_states;
    u64 alphabet_size;

    /*
     * NFA_Builder(u64 alphabet_size)
     * 用途：
     *   初始化建構器並設定字母表大小。
     * 怎麼算：
     *   記錄 alphabet_size，並以空終態 Bit_Set 起始。
     */
    NFA_Builder(u64 alphabet_size) : alphabet_size(alphabet_size), final_states(0) {};

    /*
     * mark_state_initial
     * 用途：
     *   將指定狀態標記為初始狀態。
     * 怎麼算：
     *   將狀態索引推入 initial_states 向量。
     */
    void mark_state_initial(NFA::State state) {
        this->initial_states.push_back(state);
    }

    /*
     * mark_state_final
     * 用途：
     *   標記指定狀態為終態，必要時擴充 Bit_Set。
     * 怎麼算：
     *   呼叫 grow_and_set_bit 在對應索引設 1，確保容量足夠。
     */
    void mark_state_final(NFA::State state) {
        this->final_states.grow_and_set_bit(state);
    }

    void add_transition(NFA::State source_state, u64 symbol, NFA::State destination);
    NFA build(s64 state_cnt = -1);
};


struct Macrostate {
    using State = u64;

    std::vector<State> state_names; // Accelerates computation of Post
    Bit_Set            state_set;   // To test what states are present in the macrostate
    State              handle;      // Assigned after initialization

    /*
     * Macrostate(u64 state_cnt)
     * 用途：
     *   建立空巨狀態並設定底層 Bit_Set 尺寸。
     * 怎麼算：
     *   初始化空的 state_names，並以 state_cnt 建立 state_set。
     */
    Macrostate(u64 state_cnt) : state_names({}), state_set(state_cnt) {}
    /*
     * Macrostate(u64 size, const std::vector<State>& content)
     * 用途：
     *   以指定狀態集合初始化巨狀態。
     * 怎麼算：
     *   直接保存 state_names 並以 content 建立對應的 Bit_Set。
     */
    Macrostate(u64 size, const std::vector<State>& content) : state_names(content), state_set(size, content) {}

    /*
     * Macrostate(Macrostate&& other)
     * 用途：
     *   移動建構，轉移 other 的內部集合與 handle。
     * 怎麼算：
     *   直接複製成員並保留 other 的資料來源，供之後自行管理。
     */
    Macrostate(Macrostate&& other) : state_names(other.state_names), state_set(other.state_set), handle(other.handle) {}

    /*
     * Macrostate(const Macrostate& other)
     * 用途：
     *   拷貝建構巨狀態。
     * 怎麼算：
     *   逐一複製 state_names、state_set 與 handle。
     */
    Macrostate(const Macrostate& other) : state_names(other.state_names), state_set(other.state_set), handle(other.handle) {}

    /*
     * operator<
     * 用途：
     *   定義巨狀態的排序順序以便放入關聯容器。
     * 怎麼算：
     *   直接比較底層 state_set 的字典序結果。
     */
    bool operator<(const Macrostate& other) const {
        return state_set < other.state_set;
    }

    /*
     * operator=(Macrostate&& other)
     * 用途：
     *   移動指派，接手 other 的內容。
     * 怎麼算：
     *   轉移 handle、state_set 與 state_names，並回傳自身。
     */
    Macrostate& operator=(Macrostate&& other) {
        this->handle      = other.handle;
        this->state_set   = std::move(other.state_set);
        this->state_names = std::move(other.state_names);

        return *this;
    }

    /*
     * operator=(const Macrostate& other)
     * 用途：
     *   拷貝指派，將內容複製自 other。
     * 怎麼算：
     *   複製 handle 與兩個集合成員，回傳自身。
     */
    Macrostate& operator=(const Macrostate& other) {
        this->handle      = other.handle;
        this->state_set   = other.state_set;
        this->state_names = other.state_names;

        return *this;
    }

    /*
     * empty
     * 用途：
     *   判斷巨狀態是否不含任何原子狀態。
     * 怎麼算：
     *   檢查 state_names 是否為空。
     */
    bool empty() const {
        return state_names.empty();
    }

    /*
     * populate_state_names_from_set
     * 用途：
     *   依據 Bit_Set 內容重新填充 state_names。
     * 怎麼算：
     *   清空原列表，掃描 state_set 中為 1 的索引並依序推入 state_names。
     */
    void populate_state_names_from_set() {
        this->state_names.clear();
        for (State state = 0; state < this->state_set.size; state++) {
            if (this->state_set.get_bit_value(state)) {
                this->state_names.push_back(state);
            }
        }
    }

    /*
     * to_string
     * 用途：
     *   生成描述巨狀態內容的字串。
     * 怎麼算：
     *   使用 stringstream 組合 state_names 與 handle，回傳格式化結果。
     */
    std::string to_string() const {
        std::stringstream ss;
        ss << "Macrostate " << state_names << " (handle=" << handle << ")";
        return ss.str();
    }
};

template <typename T>
struct Worklist_Construction_Context {
    std::map<T, u64>      handles;
    std::vector<const T*> worklist;

    /*
     * extract
     * 用途：
     *   從待處理清單取出最後加入的元素並回傳指標。
     * 怎麼算：
     *   取得 worklist 的末尾指標、彈出該元素，再回傳。
     */
    const T* extract() {
        auto elem = worklist.back();
        worklist.pop_back();
        return elem;
    }

    /*
     * mark_discovery
     * 用途：
     *   為新發現的元素分配 handle 並加入工作清單；若已存在則僅回傳舊元素指標。
     * 怎麼算：
     *   將 discovery.handle 設為當前 handles 大小，嘗試插入 map；若成功則推入 worklist，最後回傳 map 中的元素指標。
     */
    const T* mark_discovery(T& discovery) {
        discovery.handle = this->handles.size();
        auto [insert_pos, was_inserted] = this->handles.emplace(discovery, discovery.handle);

        const T* result_ptr = &insert_pos->first;

        if (was_inserted) {
            this->worklist.push_back(result_ptr);
        }

        return result_ptr;
    }

    /*
     * has_more_to_explore
     * 用途：
     *   檢查工作清單是否仍有待處理元素。
     * 怎麼算：
     *   回傳 worklist 是否為非空。
     */
    bool has_more_to_explore() const {
        return !this->worklist.empty();
    }
};

struct State_Pair {
    u64 first, second;
    u64 pred = -1;
    u64 handle = -1;

    /**
     * 字典序排序：
     *   依照 first、second 的順序比較，供在排序容器中安定排列。
     */
    bool operator<(const State_Pair& other) const {
        if (this->first < other.first)
            return true;
        else if (this->first > other.first)
            return false;
        return this->second < other.second;
    }

    bool operator==(const State_Pair& other) const {
        return first == other.first && second == other.second;
    }
};

std::ostream& operator<<(std::ostream& os, const State_Pair& state);
std::ostream& operator<<(std::ostream& os, const Macrostate& macrostate);

bool are_two_complete_dfas_equivalent(const NFA& first_nfa, NFA& other_nfa);
