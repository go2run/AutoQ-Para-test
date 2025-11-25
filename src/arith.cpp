#include <cassert>
#include <cmath>
#include <cstdlib>
#include <new>
#include <ostream>
#include <iostream>
#include <vector>

#include "arith.hpp"


/*
 * operator<< (Complex_Number)
 * 用途：
 *   將 Complex_Number 物件格式化輸出到流中。
 * 怎麼算：
 *   輸出實部、虛部符號和虛部絕對值，格式為 "實部 ±虛部絕對值 i"。
 */
std::ostream& operator<<(std::ostream& os, const Complex_Number& number) {
    const char* sign = number.im > 0 ? "+" : "-";

    os << number.real << " "
       << sign
       << std::abs(number.im) << " i";

    return os;
}

/*
 * operator<< (Algebraic_Complex_Number)
 * 用途：
 *   將 Algebraic_Complex_Number 物件格式化輸出到流中。
 * 怎麼算：
 *   以 5 元組形式輸出 (a, b, c, d, k)，每個分量轉換為 64 位整數。
 */
std::ostream& operator<<(std::ostream& os, const Algebraic_Complex_Number& number) {
    os << "("
       << mpz_get_si(number.a) << ", "
       << mpz_get_si(number.b) << ", "
       << mpz_get_si(number.c) << ", "
       << mpz_get_si(number.d) << ", "
       << mpz_get_si(number.k) << ")";
    return os;
}

/*
 * operator<< (Fixed_Precision_ACN)
 * 用途：
 *   將 Fixed_Precision_ACN 物件格式化輸出到流中。
 * 怎麼算：
 *   以 5 元組形式輸出 (a, b, c, d, k)，每個分量為 64 位整數。
 */
std::ostream& operator<<(std::ostream& os, const Fixed_Precision_ACN& number) {
    os << "("
       << number.a << ", "
       << number.b << ", "
       << number.c << ", "
       << number.d << ", "
       << number.k << ")";
    return os;
}

/*
 * add_row_to_row_echelon_matrix_no_copy
 * 用途：
 *   將一列添加到行簡化列階梯形矩陣中，對輸入列進行原位修改（無複製）。
 *   這是高斯消元法的核心步驟。
 * 怎麼算：
 *   若輸入列全為 0 則返回 -1。否則，與矩陣的每一列進行比較，根據主元位置決定插入位置。
 *   若在某列找到相同主元位置，執行行消元：用矩陣中該列與輸入列的主元係數，
 *   對輸入列進行線性組合以消去該位置的元素。迭代直到找到插入位置或化為全 0。
 *   返回插入行的位置，或若被化為 0 則返回 -1。
 */
s64 add_row_to_row_echelon_matrix_no_copy(ACN_Matrix& matrix, ACN_Matrix& row) {
    assert(matrix.width == row.width);

    if (row.contains_only_zeros()) return -1;

    u64 inserted_row_pivot_idx = row.find_nonzero_elem_in_row(0);

    s64 row_inserted_at = -1; // Not inserted

    for (u64 row_idx = 0; row_idx < matrix.height; row_idx++) {
        u64 row_pivot_idx = matrix.find_nonzero_elem_in_row(row_idx);

        if (row_pivot_idx < inserted_row_pivot_idx) {
            continue;
        }

        if (row_pivot_idx > inserted_row_pivot_idx) {
            matrix.insert_row_at(row, row_idx);
            row_inserted_at = row_idx;
            return row_inserted_at;
        }

        /*
         * 從矩陣的當前列開始減法，並繼續迴圈
         */
        auto& matrix_pivot = matrix.at(row_idx, row_pivot_idx);
        auto& row_pivot    = row.at(0, inserted_row_pivot_idx);

        row.subtract_from_ith_row(
            0,             // Which row of the `row` matrix to subtract from
            matrix_pivot,  // What coefficient should multiply the row before multiplication
            matrix,        // Matrix containing rows that we can subtract from
            row_idx,       // Which of the matrix rows are we subtracting
            row_pivot);    // What coeffcient should multiply the row

        inserted_row_pivot_idx = row.find_nonzero_elem_in_row(0);

        if (inserted_row_pivot_idx >= row.width) {
            return -1; /*
                        * 未插入，高斯消元將該列化簡為 0
                        */
        }
    }

    assert(false);  /*
                     * 不可達 - 該列應為獨立列而被插入，或矩陣滿秩而不需插入
                     */
}

/*
 * add_row_to_row_echelon_matrix
 * 用途：
 *   將一列添加到行簡化列階梯形矩陣中。先對輸入列進行複製，然後調用 no_copy 版本。
 * 怎麼算：
 *   建立輸入列的本地複製以保護原始數據，然後調用 add_row_to_row_echelon_matrix_no_copy 進行處理。
 */
s64 add_row_to_row_echelon_matrix(ACN_Matrix& matrix, const ACN_Matrix& row) {
    /*
     * 建立本地複製，因為將對其進行修改
     */
    ACN_Matrix row_copy (1, row.width);
    for (auto row_elem_idx = 0; row_elem_idx < row.width; row_elem_idx++) {
        auto& elem_value = row.at(0, row_elem_idx);
        row_copy.set(0, row_elem_idx, elem_value);
    }
    return add_row_to_row_echelon_matrix_no_copy(matrix, row_copy);
}

/*
 * compute_square_matrix_dim_from_1d_repr
 * 用途：
 *   從一維表示中計算方陣的維度。
 * 怎麼算：
 *   計算向量大小的平方根，驗證其平方等於原大小（確保是完全平方數），然後返回根值。
 */
u64 compute_square_matrix_dim_from_1d_repr(const std::vector<s64>& repr) {
    u64 dimension = static_cast<s64>(std::sqrt(repr.size()));

    assert(dimension*dimension == repr.size());

    return dimension;
}

/*
 * square_acn_matrix_from_ints
 * 用途：
 *   從一維 64 位整數向量建立方陣，將整數值放入代數複數的 a 分量。
 * 怎麼算：
 *   計算方陣維度，為每個整數建立 ACN 物件，只設置 a 分量。
 */
ACN_Matrix square_acn_matrix_from_ints(const std::vector<s64>& ints) {
    u64 dim = compute_square_matrix_dim_from_1d_repr(ints);

    ACN_Matrix result(dim, dim);
    Algebraic_Complex_Number* matrix_slots = new Algebraic_Complex_Number[dim*dim];

    for (u64 elem_idx = 0; elem_idx < ints.size(); elem_idx++) {
        s64 elem_int_value = ints[elem_idx];
        mpz_set_si(matrix_slots[elem_idx].a, elem_int_value);
    }

    result.data = matrix_slots;
    return result;
}

/*
 * row_from_ints
 * 用途：
 *   從 64 位整數向量建立 1 列矩陣，將整數值放入代數複數的 a 分量。
 * 怎麼算：
 *   為每個整數分量建立 ACN 物件，只設置 a 值，建立高度為向量大小、寬度為 1 的矩陣。
 */
ACN_Matrix row_from_ints(const std::vector<s64>& row_data) {
    Algebraic_Complex_Number* matrix_slots = new Algebraic_Complex_Number[row_data.size()];

    for (u64 elem_idx = 0; elem_idx < row_data.size(); elem_idx++) {
        s64 elem_int_value = row_data[elem_idx];
        mpz_set_si(matrix_slots[elem_idx].a, elem_int_value);
    }

    ACN_Matrix result(1, row_data.size(), matrix_slots);
    return result;
}

/*
 * column_from_ints
 * 用途：
 *   從 64 位整數向量建立 1 行矩陣，將整數值放入代數複數的 a 分量。
 * 怎麼算：
 *   為每個整數分量建立 ACN 物件，只設置 a 值，建立高度為 1、寬度為向量大小的矩陣。
 */
ACN_Matrix column_from_ints(const std::vector<s64>& column_data) {
    Algebraic_Complex_Number* matrix_slots = new Algebraic_Complex_Number[column_data.size()];

    for (u64 elem_idx = 0; elem_idx < column_data.size(); elem_idx++) {
        s64 elem_int_value = column_data[elem_idx];
        mpz_set_si(matrix_slots[elem_idx].a, elem_int_value);
    }

    ACN_Matrix result(column_data.size(), 1, matrix_slots);
    return result;
}


/*
 * convert_acn_into_direct_repr
 * 用途：
 *   將 Algebraic_Complex_Number 轉換為 Direct_ACN 表示。Direct_ACN 使用不同的基準，
 *   使計算更直接但精度略低於 ACN。
 * 怎麼算：
 *   提取 a, b, c, d 分量，進行基準轉換計算：b' = b - d, c 不變, d' = b + d。
 *   根據 k 值的奇偶性，若為奇數需乘以 sqrt(2) 的等效轉換。
 *   最後計算所有分量的 OR 以求公共尾零位數進行正規化，更新 k 值。
 */
Direct_ACN convert_acn_into_direct_repr(const Algebraic_Complex_Number& number) {
    if (number.is_zero()) {
        return {}; /*
                    * 全為 0
                    */
    }

    s64 a = mpz_get_si(number.a);
    s64 b = mpz_get_si(number.b) - mpz_get_si(number.d);
    s64 c = mpz_get_si(number.c);
    s64 d = mpz_get_si(number.b) + mpz_get_si(number.d);

    s64 input_k = mpz_get_si(number.k);
    s64 k = input_k / 2;

    if (input_k % 2) {
        /*
         * 為避免整數除法，當數值較小時縮小縮放因子，使其為偶數
         */
        k += (input_k > 0);

        /*
         * 乘以 sqrt(2) 的等效轉換
         */
        s64 new_a = b;
        s64 new_b = 2*a;
        s64 new_c = d;
        s64 new_d = 2*c;

        a = new_a;
        b = new_b;
        c = new_c;
        d = new_d;
    }

    /*
     * 計算尾零位數以進行正規化 K 值
     */
    s64 product = a | b | c | d;
    u64 trailing_zeros = 0;
    while (!(product % 2)) { /*
                              * 直到最後一位為 1
                              */
        trailing_zeros += 1;
        product = product >> 1;
    }

    a = a >> trailing_zeros;
    b = b >> trailing_zeros;
    c = c >> trailing_zeros;
    d = d >> trailing_zeros;
    k = k - trailing_zeros;

    return {.a = a, .b = b, .c = c, .d = d, .k = k};
}

/*
 * operator<< (Direct_ACN)
 * 用途：
 *   將 Direct_ACN 物件格式化輸出到流中。
 * 怎麼算：
 *   輸出格式為 (a + b/sqrt(2) + i*c + i*d/sqrt(2)) * (1/2)^k，
 *   其中 b, c, d 的符號與絕對值分開輸出。
 */
std::ostream& operator<<(std::ostream& os, const Direct_ACN& number) {
    const char* b_sign = number.b < 0 ? " - " : " + ";
    const char* c_sign = number.c < 0 ? " - " : " + ";
    const char* d_sign = number.d < 0 ? " - " : " + ";

    os << "("
       << number.a
       << b_sign << std::abs(number.b) << "*/sqrt(2)"
       << c_sign << std::abs(number.c) << "*i"
       << d_sign << std::abs(number.d) << "*i/sqrt(2)"
       << ") * (1/2)^(" << number.k << ")";

    return os;
}


/*
 * operator<< (ACN_Matrix)
 * 用途：
 *   將 ACN_Matrix 物件格式化輸出到流中。
 * 怎麼算：
 *   以方括號括住，若矩陣高度 > 1 則每列另起一行。逐個輸出所有元素，以逗號分隔。
 */
std::ostream& operator<<(std::ostream& os, const ACN_Matrix& matrix) {
    os << "[";
    if (matrix.height > 1) os << "\n";
    for (u64 row_idx = 0; row_idx < matrix.height; row_idx++) {
        for (u64 col_idx = 0; col_idx < matrix.width; col_idx++) {
            os << matrix.at(row_idx, col_idx) << ", ";
        }

        if (matrix.height > 1) os << "\n";
    }
    os << "]";

    return os;
}

