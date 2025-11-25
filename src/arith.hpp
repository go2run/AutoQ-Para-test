#pragma once

#include <cassert>
#include <cmath>
#include <gmp-x86_64.h>
#include <iostream>
#include <ostream>
#include <vector>

#include "gmp.h"
#include "basics.hpp"

/*
 * Complex_Number
 * 用途：
 *   用浮點數表示複數，儲存實部和虛部。
 * 怎麼算：
 *   real：複數的實部
 *   im：複數的虛部
 */
struct Complex_Number {
    float real, im;
};

/*
 * Fixed_Precision_ACN
 * 用途：
 *   用固定精度整數表示代數複數，儲存 5 個 64 位整數分量。
 * 怎麼算：
 *   代數複數 = (1/sqrt(2))^k * (a + b*Omega + c*Omega^2 + d*Omega^3)
 *   其中 Omega = e^(PI*i/4)，只使用 64 位整數，不使用任意精度。
 */
struct Fixed_Precision_ACN {
    s64 a, b, c, d, k;

    Fixed_Precision_ACN(s64 a, s64 b, s64 c, s64 d, s64 k) :
        a(a), b(b), c(c), d(d), k(k) {}

    bool operator==(const Fixed_Precision_ACN& other) const {
        return (a == other.a) && (b == other.b) && (c == other.c) && (d == other.d) && (k == other.k);
    }
};

std::ostream& operator<<(std::ostream& os, const Complex_Number& number);
std::ostream& operator<<(std::ostream& os, const Fixed_Precision_ACN& number);

/*
 * Algebraic_Complex_Number
 * 用途：
 *   任意精度代數複數表示。複數 C 表示為 5 元組 (a, b, c, d, k)，
 *   滿足 C = (1/sqrt(2))^k * (a + b*Omega + c*Omega^2 + d*Omega^3)，
 *   其中 Omega = e^(PI*i/4)。
 * 怎麼算：
 *   使用 GMP 任意精度整數儲存各分量，透過四元數乘法規則進行運算。
 */
struct Algebraic_Complex_Number {
    mpz_t a, b, c, d, k;

    Algebraic_Complex_Number() {
        mpz_inits(a, b, c, d, k, 0);
    }

    Algebraic_Complex_Number(s64 a, s64 b, s64 c, s64 d, s64 k = 1) {
        mpz_init_set_si(this->a, a);
        mpz_init_set_si(this->b, b);
        mpz_init_set_si(this->c, c);
        mpz_init_set_si(this->d, d);
        mpz_init_set_si(this->k, k);
    }

    Algebraic_Complex_Number(const Algebraic_Complex_Number& other) {
        mpz_init_set(this->a, other.a);
        mpz_init_set(this->b, other.b);
        mpz_init_set(this->c, other.c);
        mpz_init_set(this->d, other.d);
        mpz_init_set(this->k, other.k);
    }

    Algebraic_Complex_Number(Algebraic_Complex_Number&& other) {
        mpz_init_set(this->a, other.a);
        mpz_init_set(this->b, other.b);
        mpz_init_set(this->c, other.c);
        mpz_init_set(this->d, other.d);
        mpz_init_set(this->k, other.k);
    };

    ~Algebraic_Complex_Number() {
        mpz_clears(a, b, c, d, k, 0);
    }

    /*
     * 常用常數：
     * ONE_OVER_SQRT2() - 返回 1/sqrt(2)
     * ONE() - 返回 1
     * ZERO() - 返回 0
     */
    static Algebraic_Complex_Number ONE_OVER_SQRT2() {
        return Algebraic_Complex_Number(1, 0, 0, 0, 1);
    }

    static Algebraic_Complex_Number ONE() {
        return Algebraic_Complex_Number(1, 0, 0, 0, 0);
    }

    static Algebraic_Complex_Number ZERO() {
        return Algebraic_Complex_Number(0, 0, 0, 0, 0);
    }

    /*
     * operator*
     * 用途：
     *   將兩個代數複數相乘，計算結果為新的代數複數。
     * 怎麼算：
     *   使用四元數乘法規則，a + b*Omega + c*Omega^2 + d*Omega^3 的乘法會涉及 Omega 的冪次運算。
     *   結果的 k 值為兩個輸入 k 值之和。計算出新的 a, b, c, d 值後，若結果為 0 則設 k = 0。
     */
    Algebraic_Complex_Number operator*(const Algebraic_Complex_Number& other) const {

        Algebraic_Complex_Number result;

        mpz_t immediate;
        mpz_init(immediate);

        mpz_add(result.k, this->k, other.k);

        { // (this->a*other.a) - (this->b*other.d) - (this->c*other.c) - (this->d*other.b)
            mpz_mul(immediate, this->a, other.a);
            mpz_add(result.a, result.a, immediate);

            mpz_mul(immediate, this->b, other.d);
            mpz_sub(result.a, result.a, immediate);

            mpz_mul(immediate, this->c, other.c);
            mpz_sub(result.a, result.a, immediate);

            mpz_mul(immediate, this->d, other.b);
            mpz_sub(result.a, result.a, immediate);
        }

        { // (this->a*other.b) + (this->b*other.a) - (this->c*other.d) - (this->d*other.c)
            mpz_mul(immediate, this->a, other.b);
            mpz_add(result.b, result.b, immediate);

            mpz_mul(immediate, this->b, other.a);
            mpz_add(result.b, result.b, immediate);

            mpz_mul(immediate, this->c, other.d);
            mpz_sub(result.b, result.b, immediate);

            mpz_mul(immediate, this->d, other.c);
            mpz_sub(result.b, result.b, immediate);
        }

        { // (this->a*other.c) + (this->b*other.b) + (this->c*other.a) - (this->d*other.d)
            mpz_mul(immediate, this->a, other.c);
            mpz_add(result.c, result.c, immediate);

            mpz_mul(immediate, this->b, other.b);
            mpz_add(result.c, result.c, immediate);

            mpz_mul(immediate, this->c, other.a);
            mpz_add(result.c, result.c, immediate);

            mpz_mul(immediate, this->d, other.d);
            mpz_sub(result.c, result.c, immediate);
        }

        { // result_d = (this->a*other.d) + (this->b*other.c) + (this->c*other.b) + (this->d*other.a);
            mpz_mul(immediate, this->a, other.d);
            mpz_add(result.d, result.d, immediate);

            mpz_mul(immediate, this->b, other.c);
            mpz_add(result.d, result.d, immediate);

            mpz_mul(immediate, this->c, other.b);
            mpz_add(result.d, result.d, immediate);

            mpz_mul(immediate, this->d, other.a);
            mpz_add(result.d, result.d, immediate);
        }

        mpz_clear(immediate);

        if (result.is_zero()) {
            mpz_set_si(result.k, 0);
        }

        return result;
    }

    /*
     * rescale
     * 用途：
     *   將代數複數重新縮放至指定的 k 值。若目標 k 值與原 k 值之差為奇數，
     *   需要額外乘以 sqrt(2) 來轉換表達式。
     * 怎麼算：
     *   先計算 k 的差值，用 mpz_mul_2exp 將各分量乘以 2^(差值/2)。
     *   若差值為奇數，執行乘以 sqrt(2) 的轉換：新_a = -(b+d), 新_b = a+c, 新_c = b+d, 新_d = c-a。
     */
    Algebraic_Complex_Number rescale(const mpz_t larger_k) const {
        assert(mpz_cmp(larger_k, this->k) >= 0);

        s64 scale_difference_int = mpz_get_si(larger_k) - mpz_get_si(this->k);
        s64 half_scale_diff = scale_difference_int / 2;

        Algebraic_Complex_Number rescaled;

        mpz_mul_2exp(rescaled.a, this->a, half_scale_diff);
        mpz_mul_2exp(rescaled.b, this->b, half_scale_diff);
        mpz_mul_2exp(rescaled.c, this->c, half_scale_diff);
        mpz_mul_2exp(rescaled.d, this->d, half_scale_diff);

        if (scale_difference_int % 2) { /*
                                        * 若需要，乘以 sqrt(2)
                                        */
            mpz_t imm;
            mpz_init(imm);

            Algebraic_Complex_Number multiplied_by_sqrt2;

            /*
             * 使用 scale_difference 作為中間值
             * 計算：new_a = -rescaled.b - rescaled.d
             */
            mpz_set_ui(imm, 0);
            mpz_sub(imm, imm, rescaled.b);
            mpz_sub(multiplied_by_sqrt2.a, imm, rescaled.d);

            /*
             * 計算：new_b = rescaled.a + rescaled.c
             */
            mpz_add(multiplied_by_sqrt2.b, rescaled.a, rescaled.c);

            /*
             * 計算：new_c = rescaled.b + rescaled.d
             */
            mpz_add(multiplied_by_sqrt2.c, rescaled.b, rescaled.d);

            /*
             * 計算：new_d = rescaled.c - rescaled.a
             */
            mpz_sub(multiplied_by_sqrt2.d, rescaled.c, rescaled.a);

            mpz_clear(imm);

            return multiplied_by_sqrt2;
        }

        return rescaled;
    }

    /*
     * operator-
     * 用途：
     *   對代數複數執行一元負運算，將所有分量反號。
     * 怎麼算：
     *   使用 mpz_neg 對 a, b, c, d 各分量反號，k 值保持不變，返回結果。
     */
    Algebraic_Complex_Number operator-() const {
        Algebraic_Complex_Number result;

        mpz_neg(result.a, this->a);
        mpz_neg(result.b, this->b);
        mpz_neg(result.c, this->c);
        mpz_neg(result.d, this->d);
        mpz_set(result.k, this->k);

        return result;
    }

    /*
     * operator+
     * 用途：
     *   將兩個代數複數相加。由於 k 值不同的數無法直接相加，需要先統一縮放至相同 k 值。
     * 怎麼算：
     *   判斷兩個數的 k 值大小，將 k 值較大的數重新縮放至 k 值較小者。
     *   對縮放後的分量進行元素級加法：a+a', b+b', c+c', d+d'，結果 k 值為較小者。
     *   若結果為 0 則設 k = 0。
     */
    Algebraic_Complex_Number operator+(const Algebraic_Complex_Number& other) const {
        const Algebraic_Complex_Number* smaller = this;
        const Algebraic_Complex_Number* larger  = &other;
        if (mpz_cmp(larger->k, smaller->k) >= 1) {
            std::swap(smaller, larger);
        }

        Algebraic_Complex_Number larger_rescaled = larger->rescale(smaller->k);

        Algebraic_Complex_Number result;

        mpz_add(result.a, smaller->a, larger_rescaled.a);
        mpz_add(result.b, smaller->b, larger_rescaled.b);
        mpz_add(result.c, smaller->c, larger_rescaled.c);
        mpz_add(result.d, smaller->d, larger_rescaled.d);
        mpz_set(result.k, smaller->k);

        if (result.is_zero()) {
            mpz_set_si(result.k, 0);
        }

        return result;
    }

    Algebraic_Complex_Number operator-(const Algebraic_Complex_Number& other) const {
        auto other_negated = -other;
        return *this + other_negated;
    }

    void operator=(const Algebraic_Complex_Number& other) {
        mpz_set(this->a, other.a);
        mpz_set(this->b, other.b);
        mpz_set(this->c, other.c);
        mpz_set(this->d, other.d);
        mpz_set(this->k, other.k);
    }

    void operator+=(const Algebraic_Complex_Number& other) {
        auto addition_result = *this + other;

        mpz_set(this->a, addition_result.a);
        mpz_set(this->b, addition_result.b);
        mpz_set(this->c, addition_result.c);
        mpz_set(this->d, addition_result.d);
        mpz_set(this->k, addition_result.k);
    }

    bool operator==(const Algebraic_Complex_Number& other) const {
        return (mpz_cmp(this->a, other.a) == 0) &&
               (mpz_cmp(this->b, other.b) == 0) &&
               (mpz_cmp(this->c, other.c) == 0) &&
               (mpz_cmp(this->d, other.d) == 0) &&
               (mpz_cmp(this->k, other.k) == 0);
    }

    bool operator!=(const Algebraic_Complex_Number& other) const {
        return !(*this == other);
    }

    /*
     * is_zero
     * 用途：
     *   檢查代數複數是否為 0。
     * 怎麼算：
     *   只需檢查 a, b, c, d 四個分量是否全部為 0，若全為 0 則返回 true，否則返回 false。
     */
    bool is_zero() const {
        return (mpz_cmp_ui(this->a, 0) == 0) &&
               (mpz_cmp_ui(this->b, 0) == 0) &&
               (mpz_cmp_ui(this->c, 0) == 0) &&
               (mpz_cmp_ui(this->d, 0) == 0);
    }

    /*
     * is_integer
     * 用途：
     *   檢查代數複數是否為整數值。
     * 怎麼算：
     *   先檢查虛部分量（b, c, d）是否全為 0，不為實數則直接返回 false。
     *   若為實數，根據 k 值判斷：k <= 0 時，若 k 的絕對值為奇數則乘以 sqrt(2) 不是整數；
     *   k > 0 時，計算 a 的 2 進制尾零位數，若不小於 k/2 則為整數。
     */
    bool is_integer() const {
        bool is_real = (mpz_cmp_ui(this->b, 0) == 0) &&
                       (mpz_cmp_ui(this->c, 0) == 0) &&
                       (mpz_cmp_ui(this->d, 0) == 0);

        if (!is_real) return false;

        s64 k = mpz_get_si(this->k);
        if (k > 0) {
            if (k & static_cast<s64>(1)) return false;
            s64 dividing_power_of_two = k >> 1;

            s64 a = mpz_get_si(this->a);
            if (a == 0) return true;

            s64 a_two_power = 0;

            while (!(a & 1)) {
                a_two_power += 1;
                a >>= 1;
            }

            return (a_two_power >= dividing_power_of_two);
        } else {
            bool multiplying_by_sqrt2 = std::abs(k) & 1;
            return !multiplying_by_sqrt2;
        }
    }

    void swap(Algebraic_Complex_Number& other) noexcept {
        mpz_swap(this->a, other.a);
        mpz_swap(this->b, other.b);
        mpz_swap(this->c, other.c);
        mpz_swap(this->d, other.d);
        mpz_swap(this->k, other.k);
    }

    /*
     * into_approx
     * 用途：
     *   將代數複數轉換為單精度浮點近似值 Complex_Number。
     * 怎麼算：
     *   根據代數複數公式計算實部和虛部的浮點值：
     *   real = a + (b-d)/sqrt(2)，im = c + (b+d)/sqrt(2)，
     *   再乘以 (1/sqrt(2))^k 進行縮放後返回。
     */
    Complex_Number into_approx() const {
        float one_over_sqrt2 = 1.0f / std::sqrt(2);

        float real = static_cast<float>(mpz_get_si(this->a)) + one_over_sqrt2*(static_cast<float>(mpz_get_si(b)) - static_cast<float>(mpz_get_si(this->d)));
        float im   = static_cast<float>(mpz_get_si(this->c)) + one_over_sqrt2*(static_cast<float>(mpz_get_si(b)) + static_cast<float>(mpz_get_si(this->d)));

        real *= std::pow(one_over_sqrt2, mpz_get_si(this->k));
        im   *= std::pow(one_over_sqrt2, mpz_get_si(this->k));

        return {.real = real, .im = im};
    }

    Fixed_Precision_ACN into_fixed_precision() const {
        return Fixed_Precision_ACN(
            mpz_get_si(this->a),
            mpz_get_si(this->b),
            mpz_get_si(this->c),
            mpz_get_si(this->d),
            mpz_get_si(this->k)
        );
    }

    /*
     * normalize
     * 用途：
     *   將代數複數正規化，提取 a, b, c, d 分量中的公因子 2，並更新 k 值以保持等值。
     * 怎麼算：
     *   先檢查是否為 0；若非 0，計算各分量最低位 1 的位置（即 2 進制尾零位數）。
     *   求四個分量的最小尾零位數作為可提取的公因子。
     *   同時，k 可提供最多 k/2 個因子 2 的提升。取兩者最小值進行除法。
     *   最後從 k 中減去 2*normalization_2pow 以補償分量的縮小。
     */
    void normalize() {
        if (this->is_zero()) {
            mpz_set_ui(this->k, 0);
        }

        u64 a_2pow = mpz_scan1(this->a, 0);
        u64 b_2pow = mpz_scan1(this->b, 0);
        u64 c_2pow = mpz_scan1(this->c, 0);
        u64 d_2pow = mpz_scan1(this->d, 0);

        u64 spare_2pow;
        {
            u64 spare_2pow_a = std::min(a_2pow, b_2pow);
            u64 spare_2pow_b = std::min(c_2pow, d_2pow);
            spare_2pow = std::min(spare_2pow_a, spare_2pow_b);
        }

        if (spare_2pow == 0) return; /*
                                       * 沒有公因子 2 整除所有分量（乘以 omega）
                                       */

        s64 available_k = mpz_get_ui(this->k) / 2;
        if (available_k < 0) return;

        u64 normalization_2pow = std::min(spare_2pow, static_cast<u64>(available_k));

        /*
         * 將 (a, b, c, d) 除以約定的 2 的冪次
         */
        mpz_div_2exp(this->a, this->a, normalization_2pow);
        mpz_div_2exp(this->b, this->b, normalization_2pow);
        mpz_div_2exp(this->c, this->c, normalization_2pow);
        mpz_div_2exp(this->d, this->d, normalization_2pow);

        /*
         * 從 k 中減去過程中使用的部分
         */
        mpz_sub_ui(this->k, this->k, 2*normalization_2pow);
    }
};

std::ostream& operator<<(std::ostream& os, const Algebraic_Complex_Number& number);

Algebraic_Complex_Number acn_zero();

/*
 * Direct_ACN
 * 用途：
 *   固定精度代數複數的直接表示形式。複數表示為：
 *   (1/2)^k * (a + b*(1/sqrt(2)) + i*c + i*d*(1/sqrt(2)))
 * 怎麼算：
 *   縮放因子 k 為必要，否則無法表示任意大小/小的數值。
 *   相比 Algebraic_Complex_Number 使用固定 64 位整數以加快計算速度。
 */
struct Direct_ACN {
    s64 a, b, c, d, k;

    bool operator==(const Direct_ACN& other) const {
        return (this->a == other.a) && (this->b == other.b) && (this->c == other.c) && (this->d == other.d) && (this->k == other.k);
    }
};

std::ostream& operator<<(std::ostream& os, const Direct_ACN& number);
Direct_ACN convert_acn_into_direct_repr(const Algebraic_Complex_Number& number);


struct ACN_Matrix {
    u64 height, width;
    Algebraic_Complex_Number* data = nullptr;

    ACN_Matrix(u64 height, u64 width, Algebraic_Complex_Number* data_ptr = nullptr) :
        height(height), width(width), data(data_ptr)
    {
        if (this->data == nullptr) {
            /*
             * 使用零初始化，所有代數複數分量預設為 0
             */
            this->data = new Algebraic_Complex_Number[this->width*this->height];
        }
    }

    ACN_Matrix(const ACN_Matrix& other) {
        this->width  = other.width;
        this->height = other.height;

        this->data   = new Algebraic_Complex_Number[this->width*this->height];

        for (u64 idx = 0; idx < other.width*other.height; idx++) {
            this->data[idx] = other.data[idx];
        }
    }

    ACN_Matrix(ACN_Matrix&& other) {
        this->width  = other.width;
        this->height = other.height;

        this->data = other.data;
        other.data = nullptr;
    }

    /*
     * operator*
     * 用途：
     *   計算兩個代數複數矩陣的乘積。左矩陣的寬度必須等於右矩陣的高度。
     * 怎麼算：
     *   結果矩陣的大小為 (左矩陣高, 右矩陣寬)。
     *   對於每個結果元素 [i][j]，計算第 i 行與第 j 列的代數複數點積：
     *   sum(this[i][k] * other[k][j]) for k = 0 to width-1。
     */
    ACN_Matrix operator*(const ACN_Matrix& other) const {
        assert(this->width == other.height);

        u64 result_height = this->height;
        u64 result_width  = other.width;

        Algebraic_Complex_Number* result_data = new Algebraic_Complex_Number[result_height*result_width];
        u64 target_cell_idx = 0;

        for (u64 row_idx = 0; row_idx < this->height; row_idx++) {
            for (u64 col_idx = 0; col_idx < other.width; col_idx++) {

                Algebraic_Complex_Number dot_product;
                for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
                    Algebraic_Complex_Number fragment = this->at(row_idx, elem_idx) * other.at(elem_idx, col_idx);
                    dot_product += fragment;
                }

                result_data[target_cell_idx] = dot_product;
                target_cell_idx += 1;
            }
        }

        return ACN_Matrix(result_height, result_width, result_data);
    }

    Algebraic_Complex_Number& at(u64 row_idx, u64 col_idx) const {
        return this->data[row_idx*this->width + col_idx];
    }

    void set(u64 row_idx, u64 col_idx, const Algebraic_Complex_Number& value) {
        this->data[row_idx*this->width + col_idx] = value;
    }

    /*
     * find_nonzero_elem_in_row
     * 用途：
     *   在指定列中找到第一個非 0 的元素的列索引（主元）。
     * 怎麼算：
     *   從左到右掃描該列的各元素，返回第一個非 0 元素的列索引。
     *   若整列都是 0 則返回 width（表示找不到）。
     */
    u64 find_nonzero_elem_in_row(u64 row_idx) const {
        assert (row_idx < this->height);

        for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
            Algebraic_Complex_Number& elem = this->at(row_idx, elem_idx);
            if (!elem.is_zero()) {
                return elem_idx;
            }
        }

        return this->width;
    }

    ACN_Matrix extract_row(u64 row_idx) const {
        ACN_Matrix row (1, this->width);
        for (u64 column_idx = 0; column_idx < this->width; column_idx++) {
            auto& elem_value = this->at(row_idx, column_idx);
            row.set(0, column_idx, elem_value);
        }
        return row;
    }


    /*
     * insert_row_at
     * 用途：
     *   在指定位置插入一列（1 列的矩陣）。若 skip_shifting_subsequent_rows 為 false，
     *   則後續列向下移動一個位置。
     * 怎麼算：
     *   從最後一列開始，逐列向下複製已存在的元素以騰出空間。
     *   然後將新列的元素複製到指定位置。
     */
    void insert_row_at(const ACN_Matrix& row, u64 row_idx, bool skip_shifting_subsequent_rows = false) {
        assert(row.width == this->width);

        if (!skip_shifting_subsequent_rows) {
            s64 last_elem_idx  = (this->height - 1) * this->width - 1;
            s64 first_elem_idx = row_idx * this->width;

            for (s64 elem_idx = last_elem_idx; elem_idx >= first_elem_idx; elem_idx--) {
                this->data[elem_idx + this->width] = this->data[elem_idx];
            }
        }

        for (u64 elem_idx = 0; elem_idx < row.width; elem_idx++) {
            this->data[row_idx*this->width + elem_idx] = row.at(0, elem_idx);
        }
    }
    /*
     * subtract_from_ith_row
     * 用途：
     *   從第 row_idx 列中減去指定列，用於高斯消元法。
     *   重要提示：
     *   函數調用中複製 row_to_subtract_coef，因為若只儲存外部數值的參考（例如矩陣元素），
     *   則減法迴圈中可能會修改該係數，導致計算錯誤。
     * 怎麼算：
     *   對每個元素計算：this[row_idx][j] * row_coef - rows_to_subtract[row_to_subtract_idx][j] * row_to_subtract_coef
     */
    void subtract_from_ith_row(u64 row_idx, Algebraic_Complex_Number& row_coef, ACN_Matrix& rows_to_subtract, u64 row_to_subtract_idx, Algebraic_Complex_Number row_to_subtract_coef) const {
        for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
            auto subtractee_weighted = this->at(row_idx, elem_idx) * row_coef;
            auto subtractor_weighted = rows_to_subtract.at(row_to_subtract_idx, elem_idx) * row_to_subtract_coef;

            // std::cout << "Subtractee: " << this->at(row_idx, elem_idx) << "*" << row_coef << " = " << subtractee_weighted << "\n";
            // std::cout << "Subtractor: " << rows_to_subtract.at(row_to_subtract_idx, elem_idx) << "*" << row_to_subtract_coef << " = " << subtractor_weighted << "\n";

            auto result_elem = subtractee_weighted - subtractor_weighted;

            this->data[row_idx*this->width + elem_idx] = result_elem;
        }
    }

    /*
     * contains_only_zeros
     * 用途：
     *   檢查矩陣中的所有元素是否都是 0。
     * 怎麼算：
     *   掃描矩陣的所有元素，若發現任何非 0 元素立即返回 false，否則返回 true。
     */
    bool contains_only_zeros() const {
        for (u64 idx = 0; idx < this->width*this->height; idx++) {
            if (!this->data[idx].is_zero()) return false;
        }
        return true;
    }

    bool operator==(const ACN_Matrix& other) const {
        if (width != other.width || height != other.height) return false;

        for (u64 idx = 0; idx < height*width; idx++) {
            if (data[idx] != other.data[idx]) return false;
        }

        return true;
    }

    ~ACN_Matrix() {
        if (this->data != nullptr) delete[] this->data;
    }
};

std::ostream& operator<<(std::ostream& os, const ACN_Matrix& matrix);

/*
 * add_row_to_row_echelon_matrix
 * 用途：
 *   將給定列添加到行簡化列階梯形矩陣中，輸入列進行複製。
 * 怎麼算：
 *   先複製輸入列，再調用 add_row_to_row_echelon_matrix_no_copy 進行高斯消元。
 */
s64 add_row_to_row_echelon_matrix(ACN_Matrix& matrix, const ACN_Matrix& row);

/*
 * add_row_to_row_echelon_matrix_no_copy
 * 用途：
 *   將給定列添加到行簡化列階梯形矩陣中。輸入列在處理過程中被修改，
 *   最終化簡結果被插入矩陣中。
 * 怎麼算：
 *   執行高斯消元，找到與目標列主元位置相同的矩陣列進行消元，
 *   最終將結果列插入正確的行位置。
 */
s64 add_row_to_row_echelon_matrix_no_copy(ACN_Matrix& matrix, ACN_Matrix& row);

ACN_Matrix square_acn_matrix_from_ints(const std::vector<s64>& ints);
ACN_Matrix row_from_ints(const std::vector<s64>& row_data);
ACN_Matrix column_from_ints(const std::vector<s64>& column_data);
