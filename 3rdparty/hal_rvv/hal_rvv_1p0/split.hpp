// This file is part of OpenCV project.
// It is subject to the license terms in the LICENSE file found in the top-level directory
// of this distribution and at http://opencv.org/license.html.
#ifndef OPENCV_HAL_RVV_SPLIT_HPP_INCLUDED
#define OPENCV_HAL_RVV_SPLIT_HPP_INCLUDED 

#include <riscv_vector.h>

namespace cv { namespace cv_hal_rvv {

#undef cv_hal_split8u
#define cv_hal_split8u cv::cv_hal_rvv::split8u
#if defined __GNUC__
__attribute__((optimize("no-tree-vectorize")))
#endif
inline int split8u(const uchar* src, uchar** dst, int len, int cn) {
    int k = cn % 4 ? cn % 4 : 4;
    int i = 0;
    int vl = __riscv_vsetvlmax_e8m8();
    if (k == 1) {
        uchar* dst0 = dst[0];
        if (cn == 1) {
            memcpy(dst0, src, len * sizeof(uchar));
        } else {
            for (; i <= len - vl; i += vl) {
                auto vec = __riscv_vlse8_v_u8m4(src + i * cn, cn, vl);
                __riscv_vse8_v_u8m4(dst0 + i, vec, vl);
            }
            #if defined(__clang__)
            #pragma clang loop vectorize(disable)
            #endif
            for (; i < len; i++)
                dst0[i] = src[i * cn];
        }
    } else if (k == 2) {
        uchar *dst0 = dst[0], *dst1 = dst[1];
        for (; i <= len - vl * 4; i += vl * 4) {
            auto a1 = __riscv_vlse8_v_u8m1(src + 0 + i * cn, cn, vl);
            auto b1 = __riscv_vlse8_v_u8m1(src + 1 + i * cn, cn, vl);
            auto a2 = __riscv_vlse8_v_u8m1(src + 0 + (i + vl) * cn, cn, vl);
            auto b2 = __riscv_vlse8_v_u8m1(src + 1 + (i + vl) * cn, cn, vl);
            auto a3 = __riscv_vlse8_v_u8m1(src + 0 + (i + vl + vl) * cn, cn, vl);
            auto b3 = __riscv_vlse8_v_u8m1(src + 1 + (i + vl + vl) * cn, cn, vl);
            auto a4 = __riscv_vlse8_v_u8m1(src + 0 + (i + vl + vl + vl) * cn, cn, vl);
            auto b4 = __riscv_vlse8_v_u8m1(src + 1 + (i + vl + vl + vl) * cn, cn, vl);
            __riscv_vse8_v_u8m1(dst0 + i, a1, vl);
            __riscv_vse8_v_u8m1(dst1 + i, b1, vl);
            __riscv_vse8_v_u8m1(dst0 + i + vl, a2, vl);
            __riscv_vse8_v_u8m1(dst1 + i + vl, b2, vl);
            __riscv_vse8_v_u8m1(dst0 + i + vl + vl, a3, vl);
            __riscv_vse8_v_u8m1(dst1 + i + vl + vl, b3, vl);
            __riscv_vse8_v_u8m1(dst0 + i + vl + vl + vl, a4, vl);
            __riscv_vse8_v_u8m1(dst1 + i + vl + vl + vl, b4, vl);
        }
        #if defined(__clang__)
        #pragma clang loop vectorize(disable)
        #endif
        for (; i < len; i++) {
            dst0[i] = src[i * cn];
            dst1[i] = src[i * cn + 1];
        }
    } else if (k == 3) {
        uchar *dst0 = dst[0], *dst1 = dst[1], *dst2 = dst[2];
        vl = __riscv_vsetvlmax_e8m8();
        for (; i <= len - vl * 4; i += vl * 4) {
            auto a = __riscv_vlse8_v_u8m1(src + 0 + i * cn, cn, vl);
            auto b = __riscv_vlse8_v_u8m1(src + 1 + i * cn, cn, vl);
            auto c = __riscv_vlse8_v_u8m1(src + 2 + i * cn, cn, vl);
            __riscv_vse8_v_u8m1(dst0 + i, a, vl);
            __riscv_vse8_v_u8m1(dst1 + i, b, vl);
            __riscv_vse8_v_u8m1(dst2 + i, c, vl);
            auto a2 = __riscv_vlse8_v_u8m1(src + 0 + (i + vl) * cn, cn, vl);
            auto b2 = __riscv_vlse8_v_u8m1(src + 1 + (i + vl) * cn, cn, vl);
            auto c2 = __riscv_vlse8_v_u8m1(src + 2 + (i + vl) * cn, cn, vl);
            __riscv_vse8_v_u8m1(dst0 + i + vl, a2, vl);
            __riscv_vse8_v_u8m1(dst1 + i + vl, b2, vl);
            __riscv_vse8_v_u8m1(dst2 + i + vl, c2, vl);
            auto a3 = __riscv_vlse8_v_u8m1(src + 0 + (i + vl + vl) * cn, cn, vl);
            auto b3 = __riscv_vlse8_v_u8m1(src + 1 + (i + vl + vl) * cn, cn, vl);
            auto c3 = __riscv_vlse8_v_u8m1(src + 2 + (i + vl + vl) * cn, cn, vl);
            __riscv_vse8_v_u8m1(dst0 + i + vl + vl, a3, vl);
            __riscv_vse8_v_u8m1(dst1 + i + vl + vl, b3, vl);
            __riscv_vse8_v_u8m1(dst2 + i + vl + vl, c3, vl);
            auto a4 = __riscv_vlse8_v_u8m1(src + 0 + (i + vl + vl + vl) * cn, cn, vl);
            auto b4 = __riscv_vlse8_v_u8m1(src + 1 + (i + vl + vl + vl) * cn, cn, vl);
            auto c4 = __riscv_vlse8_v_u8m1(src + 2 + (i + vl + vl + vl) * cn, cn, vl);
            __riscv_vse8_v_u8m1(dst0 + i + vl + vl + vl, a4, vl);
            __riscv_vse8_v_u8m1(dst1 + i + vl + vl + vl, b4, vl);
            __riscv_vse8_v_u8m1(dst2 + i + vl + vl + vl, c4, vl);
        }
        #if defined(__clang__)
        #pragma clang loop vectorize(disable)
        #endif
        for (; i < len; i++) {
            dst0[i] = src[i * cn];
            dst1[i] = src[i * cn + 1];
            dst2[i] = src[i * cn + 2];
        }
    } else {
        uchar *dst0 = dst[0], *dst1 = dst[1], *dst2 = dst[2], *dst3 = dst[3];
        vl = __riscv_vsetvlmax_e8m8();
        for (; i <= len - vl * 4; i += vl * 4) {
            auto a = __riscv_vlse8_v_u8m1(src + 0 + i * cn, cn, vl);
            auto b = __riscv_vlse8_v_u8m1(src + 1 + i * cn, cn, vl);
            auto c = __riscv_vlse8_v_u8m1(src + 2 + i * cn, cn, vl);
            auto d = __riscv_vlse8_v_u8m1(src + 3 + i * cn, cn, vl);
            __riscv_vse8_v_u8m1(dst0 + i, a, vl);
            __riscv_vse8_v_u8m1(dst1 + i, b, vl);
            __riscv_vse8_v_u8m1(dst2 + i, c, vl);
            __riscv_vse8_v_u8m1(dst3 + i, d, vl);

            auto a2 = __riscv_vlse8_v_u8m1(src + 0 + (i + vl) * cn, cn, vl);
            auto b2 = __riscv_vlse8_v_u8m1(src + 1 + (i + vl) * cn, cn, vl);
            auto c2 = __riscv_vlse8_v_u8m1(src + 2 + (i + vl) * cn, cn, vl);
            auto d2 = __riscv_vlse8_v_u8m1(src + 3 + (i + vl) * cn, cn, vl);
            __riscv_vse8_v_u8m1(dst0 + i + vl, a2, vl);
            __riscv_vse8_v_u8m1(dst1 + i + vl, b2, vl);
            __riscv_vse8_v_u8m1(dst2 + i + vl, c2, vl);
            __riscv_vse8_v_u8m1(dst3 + i + vl, d2, vl);

            auto a3 = __riscv_vlse8_v_u8m1(src + 0 + (i + vl + vl) * cn, cn, vl);
            auto b3 = __riscv_vlse8_v_u8m1(src + 1 + (i + vl + vl) * cn, cn, vl);
            auto c3 = __riscv_vlse8_v_u8m1(src + 2 + (i + vl + vl) * cn, cn, vl);
            auto d3 = __riscv_vlse8_v_u8m1(src + 3 + (i + vl + vl) * cn, cn, vl);
            __riscv_vse8_v_u8m1(dst0 + i + vl + vl, a3, vl);
            __riscv_vse8_v_u8m1(dst1 + i + vl + vl, b3, vl);
            __riscv_vse8_v_u8m1(dst2 + i + vl + vl, c3, vl);
            __riscv_vse8_v_u8m1(dst3 + i + vl + vl, d3, vl);

            auto a4 = __riscv_vlse8_v_u8m1(src + 0 + (i + vl + vl + vl) * cn, cn, vl);
            auto b4 = __riscv_vlse8_v_u8m1(src + 1 + (i + vl + vl + vl) * cn, cn, vl);
            auto c4 = __riscv_vlse8_v_u8m1(src + 2 + (i + vl + vl + vl) * cn, cn, vl);
            auto d4 = __riscv_vlse8_v_u8m1(src + 3 + (i + vl + vl + vl) * cn, cn, vl);
            __riscv_vse8_v_u8m1(dst0 + i + vl + vl + vl, a4, vl);
            __riscv_vse8_v_u8m1(dst1 + i + vl + vl + vl, b4, vl);
            __riscv_vse8_v_u8m1(dst2 + i + vl + vl + vl, c4, vl);
            __riscv_vse8_v_u8m1(dst3 + i + vl + vl + vl, d4, vl);
        }
        #if defined(__clang__)
        #pragma clang loop vectorize(disable)
        #endif
        for (; i < len; i++) {
            dst0[i] = src[i * cn];
            dst1[i] = src[i * cn + 1];
            dst2[i] = src[i * cn + 2];
            dst3[i] = src[i * cn + 3];
        }
    }
    for (; k < cn; k += 4) {
        uchar *dst0 = dst[k], *dst1 = dst[k+1], *dst2 = dst[k+2], *dst3 = dst[k+3];
        i = 0;
        vl = __riscv_vsetvlmax_e8m2();
        for (; i <= len - vl * 2; i += vl * 2) {
            auto a = __riscv_vlse8_v_u8m2(src + k + 0 + i * cn, cn, vl);
            auto b = __riscv_vlse8_v_u8m2(src + k + 1 + i * cn, cn, vl);
            auto c = __riscv_vlse8_v_u8m2(src + k + 2 + i * cn, cn, vl);
            auto d = __riscv_vlse8_v_u8m2(src + k + 3 + i * cn, cn, vl);
            __riscv_vse8_v_u8m2(dst0 + i, a, vl);
            __riscv_vse8_v_u8m2(dst1 + i, b, vl);
            __riscv_vse8_v_u8m2(dst2 + i, c, vl);
            __riscv_vse8_v_u8m2(dst3 + i, d, vl);

            auto a2 = __riscv_vlse8_v_u8m2(src + k + 0 + (i + vl) * cn, cn, vl);
            auto b2 = __riscv_vlse8_v_u8m2(src + k + 1 + (i + vl) * cn, cn, vl);
            auto c2 = __riscv_vlse8_v_u8m2(src + k + 2 + (i + vl) * cn, cn, vl);
            auto d2 = __riscv_vlse8_v_u8m2(src + k + 3 + (i + vl) * cn, cn, vl);
            __riscv_vse8_v_u8m2(dst0 + i + vl, a2, vl);
            __riscv_vse8_v_u8m2(dst1 + i + vl, b2, vl);
            __riscv_vse8_v_u8m2(dst2 + i + vl, c2, vl);
            __riscv_vse8_v_u8m2(dst3 + i + vl, d2, vl);
        }
        #if defined(__clang__)
        #pragma clang loop vectorize(disable)
        #endif
        for (; i < len; i++) {
            dst0[i] = src[k + i * cn];
            dst1[i] = src[k + i * cn + 1];
            dst2[i] = src[k + i * cn + 2];
            dst3[i] = src[k + i * cn + 3];
        }
    }
    return CV_HAL_ERROR_OK;
    }
}}
#endif
