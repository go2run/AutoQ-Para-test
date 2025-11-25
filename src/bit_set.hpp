#pragma once

#include "basics.hpp"

#include <vector>
#include <cstring>
#include <cassert>
#include <bit>  // For popcount

struct Bit_Set {
    u64  size; // Number of *bits* stored
    u64* data;

    /*
     * Bit_Set(u64 size, const std::vector<u64>& content)
     * 用途：
     *   以指定的位數建立位集合，並根據輸入的索引向量預先將對應位元設為 1。
     * 怎麼算：
     *   先計算需要的 bucket 數量並配置儲存空間，初始化為 0；
     *   接著遍歷 content，對每個索引呼叫 set_bit 以設定位元。
     */
    Bit_Set(u64 size, const std::vector<u64>& content) : size(size), data(nullptr) {
        u64 bucket_count = this->get_bucket_count();
        this->data = new u64[bucket_count];

        memset(this->data, 0, sizeof(u64)*bucket_count);

        for (auto& bit_to_set : content) {
            this->set_bit(bit_to_set, true);
        }
    }

    /*
     * Bit_Set(u64 size, u64* data_ptr = nullptr)
     * 用途：
     *   建立具有指定位數的 Bit_Set，可選擇外部提供的資料指標。
     * 怎麼算：
     *   若未提供指標且 size > 0，依 bucket 數量配置並清零；
     *   若 size 為 0 則直接將 data 設為 nullptr。
     */
    Bit_Set(u64 size, u64* data_ptr = nullptr) : size(size), data(data_ptr) {
        if (this->data == nullptr && this->size > 0) {
            u64 bucket_count = size / bits_in_bucket();
            bucket_count += (size % bits_in_bucket()) > 0;
            this->data = new u64[bucket_count];

            memset(this->data, 0, sizeof(u64)*bucket_count);
        }

        if (this->size == 0) this->data = nullptr;
    }

    /*
     * Bit_Set(const Bit_Set& other)
     * 用途：
     *   拷貝建構函式，生成與 other 相同大小與內容的 Bit_Set。
     * 怎麼算：
     *   複製 size，依 bucket 數量配置新空間，使用 memcpy 將位元資料逐 bucket 複製。
     */
    Bit_Set(const Bit_Set& other) {
        this->size = other.size;
        this->data = new u64[other.get_bucket_count()];
        std::memcpy(this->data, other.data, sizeof(u64) * other.get_bucket_count());
    }

    /*
     * Bit_Set(Bit_Set&& other)
     * 用途：
     *   移動建構函式，接管 other 的緩衝區以避免重新配置。
     * 怎麼算：
     *   複製 size 與 data 指標，並將 other 的指標置空、大小歸零以避免重複釋放。
     */
    Bit_Set(Bit_Set&& other) {
        this->size = other.size;
        this->data = other.data;

        other.data = nullptr;
        other.size = 0;
    }

    /*
     * operator=(const Bit_Set& other)
     * 用途：
     *   拷貝指派，將當前 Bit_Set 內容更新為 other。
     * 怎麼算：
     *   更新 size，重新配置足夠的 bucket，並複製 other 的位元資料。
     */
    Bit_Set& operator=(const Bit_Set& other) {
        this->size = other.size;
        this->data = new u64[other.get_bucket_count()];
        std::memcpy(this->data, other.data, sizeof(u64) * other.get_bucket_count());

        return *this;
    }

    /*
     * operator=(Bit_Set&& other)
     * 用途：
     *   移動指派，接管 other 的資料指標並釋放自身舊資源。
     * 怎麼算：
     *   更新 size，將 other 的 data 轉移至當前物件，同時清空 other 的指標與大小。
     */
    Bit_Set& operator=(Bit_Set&& other) {
        this->size = other.size;
        other.size = 0;

        this->data = other.data;
        other.data = nullptr;

        return *this;
    }

    /*
     * bits_in_bucket
     * 用途：
     *   取得單一 bucket（u64）可以儲存的位元數。
     * 怎麼算：
     *   回傳 u64 的位元數（8 位元組 * 8 位元）。
     */
    u64 bits_in_bucket() const {
        return sizeof(u64) * 8;
    }

    /*
     * operator==
     * 用途：
     *   判斷兩個 Bit_Set 是否大小相同且內容一致。
     * 怎麼算：
     *   若 size 不同直接返回 false；否則逐 bucket 比對所有儲存值。
     */
    bool operator==(const Bit_Set& other) const {
        if (this->size != other.size) return false;

        u64 bucket_cnt = this->get_bucket_count();

        for (u64 bucket_idx = 0; bucket_idx < bucket_cnt; bucket_idx++) {
            if (this->data[bucket_idx] != other.data[bucket_idx]) return false;
        }

        return true;
    }

    /*
     * popcount
     * 用途：
     *   計算位集合中被設定為 1 的位元總數。
     * 怎麼算：
     *   逐 bucket 使用 std::popcount 計算 1 的數量並累加後回傳。
     */
    u64 popcount() const {
        u64 popcount = 0;
        for (u64 bucket_idx = 0; bucket_idx < this->get_bucket_count(); bucket_idx++) {
            u64 bucket_val = this->data[bucket_idx];
            popcount += std::popcount(bucket_val);
        }
        return popcount;
    }

    /*
     * are_all_bits_set
     * 用途：
     *   檢查所有有效位元是否皆為 1。
     * 怎麼算：
     *   除最後一個 bucket 外應全部為全 1；最後一個 bucket 則需依有效位元建立遮罩比對。
     */
    bool are_all_bits_set() const {
        if (this->get_bucket_count() == 0) return true;

        for (u64 bucket_idx = 0; bucket_idx < this->get_bucket_count() - 1; bucket_idx++) {
            if (this->data[bucket_idx] != ~static_cast<u64>(0)) return false;
        }

        u64 last_bucket = this->data[this->get_bucket_count() - 1];

        u64 last_valid_bit_offset = (this->size - 1) % bits_in_bucket();
        u64 bits_to_zero_out = bits_in_bucket() - (last_valid_bit_offset + 1);

        u64 all_ones = ~static_cast<u64>(0);
        u64 mask = (all_ones >> bits_to_zero_out);

        return last_bucket == mask;
    }

    /*
     * into_vector
     * 用途：
     *   將集合中所有值為 1 的位元索引輸出成向量。
     * 怎麼算：
     *   逐位掃描 size 範圍，若 get_bit_value 為真則將索引加入結果。
     */
    std::vector<u64> into_vector() const;

    /**
     * Lexicographical ordering on the underlying bit representation
     */
    /*
     * operator<
     * 用途：
     *   依底層 bucket 值進行字典序比較，方便放入排序容器。
     * 怎麼算：
     *   假設 size 相同，逐 bucket 比對，第一個不同的 bucket 決定順序。
     */
    bool operator<(const Bit_Set& other) const {
        assert(this->size == other.size);

        for (u64 bucket_idx = 0; bucket_idx < this->get_bucket_count(); bucket_idx++) {
            if (this->data[bucket_idx] < other.data[bucket_idx]) return true;
            if (this->data[bucket_idx] > other.data[bucket_idx]) return false;
        }
        return false; // All buckets have equal value
    }

    /*
     * empty
     * 用途：
     *   判斷集合是否全為 0。
     * 怎麼算：
     *   將所有 bucket 做 OR，若結果為 0 則表示沒有任何位元被設為 1。
     */
    bool empty() const {
        u64 reduction_result = 0;
        for (u64 bucket_idx = 0; bucket_idx < this->get_bucket_count(); bucket_idx++) {
            reduction_result |= this->data[bucket_idx];
        }
        return reduction_result;
    };

    /*
     * is_intersection_empty
     * 用途：
     *   測試與另一個 Bit_Set 的交集是否為空。
     * 怎麼算：
     *   逐 bucket 做 AND，若任何 bucket 結果非零則交集非空；最後一個殘餘 bucket 需遮罩無效位元。
     */
    bool is_intersection_empty(const Bit_Set& other) const {
        u64 min_size = std::min(this->size, other.size);

        u64 chunks_to_compare = min_size / bits_in_bucket();

        for (u64 chunk_idx = 0; chunk_idx < chunks_to_compare; chunk_idx++) {
            if (this->data[chunk_idx] & other.data[chunk_idx]) return false;
        }

        u64 spilled_size = min_size % bits_in_bucket();
        if (spilled_size == 0) return true;


        u64 spilled_contents = this->data[chunks_to_compare] & other.data[chunks_to_compare];
        spilled_contents = spilled_contents & ~((~static_cast<u64>(0)) << spilled_size); // Extract only the relevant bits

        return spilled_contents == 0; // The result of and is 0, meaning that the intersection is empty
    }

    /*
     * is_superset
     * 用途：
     *   判斷當前集合是否為 other 的超集合。
     * 怎麼算：
     *   逐 bucket 檢查 this 的位元包含 other 的所有 1，最後一個 bucket 亦需套用遮罩比對。
     */
    bool is_superset(const Bit_Set& other) const {
        u64 min_size = std::min(this->size, other.size);

        u64 chunks_to_compare = min_size / bits_in_bucket();

        for (u64 chunk_idx = 0; chunk_idx < chunks_to_compare; chunk_idx++) {
            u64 interleaved = this->data[chunk_idx] & other.data[chunk_idx];

            // There is a bit set in other chunk that is not set in this
            if (interleaved != this->data[chunk_idx]) return false;
        }

        u64 spilled_size = min_size % bits_in_bucket();
        if (spilled_size == 0) return true;

        u64 mask =~((~static_cast<u64>(0)) << spilled_size); // Extract only the relevant bits
        u64 spilled_interleaved = this->data[chunks_to_compare] | other.data[chunks_to_compare];

        spilled_interleaved = spilled_interleaved & mask;

        return spilled_interleaved == (this->data[chunks_to_compare] & mask);
    }

    /*
     * calc_target_bucket
     * 用途：
     *   計算指定位元索引所屬的 bucket 編號。
     * 怎麼算：
     *   以 bits_in_bucket() 取整除，回傳商值。
     */
    u64 calc_target_bucket(u64 bit_idx) const {
        u64 bucket = bit_idx / bits_in_bucket();
        return bucket;
    }

    /*
     * set_bit
     * 用途：
     *   設定或清除指定索引的位元。
     * 怎麼算：
     *   先找到目標 bucket 與偏移，若 value 為真則以 OR 置位，否則以 AND 與反遮罩清位。
     */
    void set_bit(u64 bit_idx, bool value = true) {
        assert (bit_idx < this->size);

        u64 target_bucket_idx = this->calc_target_bucket(bit_idx);
        u64 bucket_offset     = bit_idx % bits_in_bucket();

        if (value) {
            this->data[target_bucket_idx] = this->data[target_bucket_idx] | (static_cast<u64>(1u) << bucket_offset);
        } else {
            this->data[target_bucket_idx] = this->data[target_bucket_idx] & ~(static_cast<u64>(1u) << bucket_offset);
        }
    }

   /*
    * grow
    * 用途：
    *   將 Bit_Set 的容量擴充至至少 new_size 位元。
    * 怎麼算：
    *   若現有 bucket 足夠僅更新 size；否則重新配置更大的陣列，複製舊資料並將新增 bucket 清零。
    */
   void grow(u64 new_size) {
        if (new_size <= size) return;

        u64 target_bucket_idx = this->calc_target_bucket(new_size - 1);

        if (target_bucket_idx < this->get_bucket_count()) { // We have enough buckets, no need to alloc
            this->size = new_size; // Note that we can accomodate more bits now
            return;
        }

        // There is not enough space; we need to grow the underlying storage
        u64 old_bucket_cnt = this->get_bucket_count();
        this->size = new_size;
        u64 new_bucket_cnt = this->get_bucket_count();

        u64* old_data = this->data;
        u64* new_data = new u64[new_bucket_cnt];

        for (u64 i = 0; i < old_bucket_cnt; i++) {
            new_data[i] = this->data[i];
        }

        for (u64 i = old_bucket_cnt; i < new_bucket_cnt; i++) {
            new_data[i] = 0;
        }

        this->data = new_data;
        delete[] old_data;
    }

    /*
     * grow_and_set_bit
     * 用途：
     *   在需要時自動擴充容量後設定指定位元。
     * 怎麼算：
     *   計算需要的最小尺寸，呼叫 grow 擴充，再呼叫 set_bit 寫入值。
     */
    void grow_and_set_bit(u64 bit_idx, bool bit_value = true) {
        u64 required_size = bit_idx + 1;
        this->grow(required_size);
        this->set_bit(bit_idx, bit_value);
    }

    /*
     * get_bit_value
     * 用途：
     *   讀取指定索引的位元值。
     * 怎麼算：
     *   找到目標 bucket 與偏移，取出該位元並回傳是否為 1。
     */
    bool get_bit_value(u64 bit_idx) const {
        u64 target_bucket_idx = this->calc_target_bucket(bit_idx);
        u64 bucket_offset     = bit_idx % bits_in_bucket();

        u64 bucket_value = this->data[target_bucket_idx];
        u64 bit_value = (bucket_value & (static_cast<u64>(1u) << bucket_offset)) > 0;

        return bit_value;
    }

    /*
     * set_all
     * 用途：
     *   將所有有效位元批次設定為 0 或 1。
     * 怎麼算：
     *   每個 bucket 寫入全 0 或全 1，最後一個 bucket 需位移清掉超出 size 的無效位元。
     */
    void set_all(bool value = true) {
        u64 bucket_value = value ? (~0) : 0;

        u64 bucket_cnt = this->get_bucket_count();
        for (u64 bucket_idx = 0; bucket_idx < bucket_cnt; bucket_idx++) {
            this->data[bucket_idx] = bucket_value;
        }

        u64 last_bucket_idx = bucket_cnt - 1;
        u64 last_valid_bit_offset = (this->size - 1) % bits_in_bucket();
        u64 bits_to_zero_out = bits_in_bucket() - (last_valid_bit_offset + 1);

        this->data[last_bucket_idx] = (this->data[last_bucket_idx] >> bits_to_zero_out);
    }

    /*
     * clear
     * 用途：
     *   將集合重置為全 0。
     * 怎麼算：
     *   呼叫 set_all(false) 清除所有位元。
     */
    void clear() {
        this->set_all(false);
    }

    /*
     * get_bucket_count
     * 用途：
     *   計算儲存當前 size 所需的 bucket 數。
     * 怎麼算：
     *   以 bits_in_bucket() 取整除並視餘數是否增加一個 bucket。
     */
    u64 get_bucket_count() const {
        u64 bucket_count = size / bits_in_bucket();
        bucket_count += (size %  bits_in_bucket()) > 0;

        return bucket_count;
    }

    /*
     * ~Bit_Set
     * 用途：
     *   解構函式，釋放內部配置的位元儲存空間。
     * 怎麼算：
     *   若 data 非空則 delete[]，並將指標設為 nullptr 以避免懸掛。
     */
    ~Bit_Set() {
        if (this->data != nullptr) delete[] this->data;
        this->data = nullptr;
    }
};
