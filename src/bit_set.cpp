#include "bit_set.hpp"

/*
 * Bit_Set::into_vector
 * 用途：
 *   產出所有設為 1 的位元索引清單，便於外部遍歷。
 * 怎麼算：
 *   從 0 迭代到 size-1，對每個索引呼叫 get_bit_value，為真則推入結果向量。
 */
std::vector<u64> Bit_Set::into_vector() const {
    std::vector<u64> bit_set_contents;

    for (size_t bit_idx = 0; bit_idx < this->size; bit_idx++) {
        if (this->get_bit_value(bit_idx)) {
            bit_set_contents.push_back(bit_idx);
        }
    }

    return bit_set_contents;
}
