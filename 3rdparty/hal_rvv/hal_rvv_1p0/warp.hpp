// This file is part of OpenCV project.
// It is subject to the license terms in the LICENSE file found in the top-level directory
// of this distribution and at http://opencv.org/license.html.

// Copyright (C) 2025, Institute of Software, Chinese Academy of Sciences.

#ifndef OPENCV_HAL_RVV_WARP_HPP_INCLUDED
#define OPENCV_HAL_RVV_WARP_HPP_INCLUDED

#include <riscv_vector.h>

namespace cv { namespace cv_hal_rvv {

namespace remap {
#undef cv_hal_remap32f
#define cv_hal_remap32f cv::cv_hal_rvv::remap::remap32f

class RemapInvoker : public ParallelLoopBody
{
public:
    template<typename... Args>
    RemapInvoker(std::function<int(int, int, Args...)> _func, Args&&... args)
    {
        func = std::bind(_func, std::placeholders::_1, std::placeholders::_2, std::forward<Args>(args)...);
    }

    virtual void operator()(const Range& range) const override
    {
        func(range.start, range.end);
    }

private:
    std::function<int(int, int)> func;
};

template<typename... Args>
static inline int invoke(int width, int height, std::function<int(int, int, Args...)> func, Args&&... args)
{
    cv::parallel_for_(Range(1, height), RemapInvoker(func, std::forward<Args>(args)...), static_cast<double>((width - 1) * height) / (1 << 15));
    return func(0, 1, std::forward<Args>(args)...);
}

template<typename T> struct rvv;
template<> struct rvv<uchar>
{
    static inline vfloat32m8_t vcvt0(vuint8m2_t a, size_t b) { return __riscv_vfcvt_f(__riscv_vzext_vf4(a, b), b); }
    static inline vuint8m2_t vcvt1(vfloat32m8_t a, size_t b) { return __riscv_vnclipu(__riscv_vfncvt_xu(a, b), 0, __RISCV_VXRM_RNU, b); }
    static inline vuint8m2_t vloxei(const uchar* a, vuint32m8_t b, size_t c) { return __riscv_vloxei32_v_u8m2(a, b, c); }
};
template<> struct rvv<ushort>
{
    static inline vfloat32m8_t vcvt0(vuint16m4_t a, size_t b) { return __riscv_vfwcvt_f(a, b); }
    static inline vuint16m4_t vcvt1(vfloat32m8_t a, size_t b) { return __riscv_vfncvt_xu(a, b); }
    static inline vuint16m4_t vloxei(const ushort* a, vuint32m8_t b, size_t c) { return __riscv_vloxei32_v_u16m4(a, b, c); }
};
template<> struct rvv<short>
{
    static inline vfloat32m8_t vcvt0(vint16m4_t a, size_t b) { return __riscv_vfwcvt_f(a, b); }
    static inline vint16m4_t vcvt1(vfloat32m8_t a, size_t b) { return __riscv_vfncvt_x(a, b); }
    static inline vint16m4_t vloxei(const short* a, vuint32m8_t b, size_t c) { return __riscv_vloxei32_v_i16m4(a, b, c); }
};
template<> struct rvv<float>
{
    static inline vfloat32m8_t vcvt0(vfloat32m8_t a, size_t) { return a; }
    static inline vfloat32m8_t vcvt1(vfloat32m8_t a, size_t) { return a; }
    static inline vfloat32m8_t vloxei(const float* a, vuint32m8_t b, size_t c) { return __riscv_vloxei32_v_f32m8(a, b, c); }
};

template<typename helper>
static inline int remap32fC1(int start, int end, const uchar *src_data, size_t src_step, int src_width, int src_height,
                             uchar *dst_data, size_t dst_step, int dst_width,
                             const float* mapx, size_t mapx_step, const float* mapy, size_t mapy_step,
                             int interpolation, int border_type, const double* border_value)
{
    using T = typename helper::ElemType;

    for (int i = start; i < end; i++)
    {
        int vl;
        for (int j = 0; j < dst_width; j += vl)
        {
            vl = helper::setvl(dst_width - j);
            auto mx = RVV_SameLen<float, helper>::vload(mapx + i * mapx_step + j, vl);
            auto my = RVV_SameLen<float, helper>::vload(mapy + i * mapy_step + j, vl);
            if (interpolation & CV_HAL_WARP_RELATIVE_MAP)
            {
                mx = __riscv_vfadd(mx, __riscv_vfcvt_f(__riscv_vadd(RVV_SameLen<uint, helper>::vid(vl), j, vl), vl), vl);
                my = __riscv_vfadd(my, i, vl);
            }

            auto access = [&](typename RVV_SameLen<int, helper>::VecType ix, typename RVV_SameLen<int, helper>::VecType iy) {
                auto ux = RVV_SameLen<uint, helper>::reinterpret(__riscv_vmin(__riscv_vmax(ix, 0, vl), src_width  - 1, vl));
                auto uy = RVV_SameLen<uint, helper>::reinterpret(__riscv_vmin(__riscv_vmax(iy, 0, vl), src_height - 1, vl));
                auto src = rvv<T>::vloxei(reinterpret_cast<const T*>(src_data), __riscv_vmadd(uy, src_step, __riscv_vmul(ux, sizeof(T), vl), vl), vl);
                if (border_type == CV_HAL_BORDER_CONSTANT)
                {
                    auto mask = __riscv_vmor(__riscv_vmsne(ix, RVV_SameLen<int, helper>::reinterpret(ux), vl), __riscv_vmsne(iy, RVV_SameLen<int, helper>::reinterpret(uy), vl), vl);
                    src = __riscv_vmerge(src, helper::vmv(border_value[0], vl), mask, vl);
                }
                return src;
            };
            if ((interpolation & ~CV_HAL_WARP_RELATIVE_MAP) == CV_HAL_INTER_NEAREST)
            {
                auto ix = __riscv_vfcvt_x(mx, vl), iy = __riscv_vfcvt_x(my, vl);
                helper::vstore(reinterpret_cast<T*>(dst_data + i * dst_step) + j, access(ix, iy), vl);
            }
            else
            {
                vint32m8_t ix0, iy0;
                int rd;
                asm volatile("fsrmi %0, 2 \n\t vsetvli zero,%3,e32,m8,ta,ma \n\t vfcvt.x.f.v %1,%4 \n\t vfcvt.x.f.v %2,%5 \n\t fsrm %0"
                             : "=&r"(rd), "=vr"(ix0), "=vr"(iy0)
                             : "r"(vl), "vr"(mx), "vr"(my)); // Rounding Mode: RDN
                auto ix1 = __riscv_vadd(ix0, 1, vl), iy1 = __riscv_vadd(iy0, 1, vl);
                auto v0 = rvv<T>::vcvt0(access(ix0, iy0), vl);
                auto v1 = rvv<T>::vcvt0(access(ix1, iy0), vl);
                auto v2 = rvv<T>::vcvt0(access(ix0, iy1), vl);
                auto v3 = rvv<T>::vcvt0(access(ix1, iy1), vl);

                mx = __riscv_vfsub(mx, __riscv_vfcvt_f(ix0, vl), vl);
                my = __riscv_vfsub(my, __riscv_vfcvt_f(iy0, vl), vl);
                v0 = __riscv_vfmacc(v0, mx, __riscv_vfsub(v1, v0, vl), vl);
                v2 = __riscv_vfmacc(v2, mx, __riscv_vfsub(v3, v2, vl), vl);
                v0 = __riscv_vfmacc(v0, my, __riscv_vfsub(v2, v0, vl), vl);
                helper::vstore(reinterpret_cast<T*>(dst_data + i * dst_step) + j, rvv<T>::vcvt1(v0, vl), vl);
            }
        }
    }

    return CV_HAL_ERROR_OK;
}

static inline int remap32fC3(int start, int end, const uchar *src_data, size_t src_step, int src_width, int src_height,
                             uchar *dst_data, size_t dst_step, int dst_width,
                             const float* mapx, size_t mapx_step, const float* mapy, size_t mapy_step,
                             int interpolation, int border_type, const double* border_value)
{
    for (int i = start; i < end; i++)
    {
        int vl;
        for (int j = 0; j < dst_width; j += vl)
        {
            vl = __riscv_vsetvl_e8mf2(dst_width - j);
            auto mx = __riscv_vle32_v_f32m2(mapx + i * mapx_step + j, vl);
            auto my = __riscv_vle32_v_f32m2(mapy + i * mapy_step + j, vl);
            if (interpolation & CV_HAL_WARP_RELATIVE_MAP)
            {
                mx = __riscv_vfadd(mx, __riscv_vfcvt_f(__riscv_vadd(__riscv_vid_v_u32m2(vl), j, vl), vl), vl);
                my = __riscv_vfadd(my, i, vl);
            }

            auto access = [&](vint32m2_t ix, vint32m2_t iy, vuint8mf2_t& src0, vuint8mf2_t& src1, vuint8mf2_t& src2) {
                auto ux = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vmin(__riscv_vmax(ix, 0, vl), src_width  - 1, vl));
                auto uy = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vmin(__riscv_vmax(iy, 0, vl), src_height - 1, vl));
                auto src = __riscv_vloxseg3ei32_v_u8mf2x3(src_data, __riscv_vmadd(uy, src_step, __riscv_vmul(ux, 3, vl), vl), vl);
                src0 = __riscv_vget_v_u8mf2x3_u8mf2(src, 0);
                src1 = __riscv_vget_v_u8mf2x3_u8mf2(src, 1);
                src2 = __riscv_vget_v_u8mf2x3_u8mf2(src, 2);
                if (border_type == CV_HAL_BORDER_CONSTANT)
                {
                    auto mask = __riscv_vmor(__riscv_vmsne(ix, __riscv_vreinterpret_v_u32m2_i32m2(ux), vl), __riscv_vmsne(iy, __riscv_vreinterpret_v_u32m2_i32m2(uy), vl), vl);
                    src0 = __riscv_vmerge(src0, border_value[0], mask, vl);
                    src1 = __riscv_vmerge(src1, border_value[1], mask, vl);
                    src2 = __riscv_vmerge(src2, border_value[2], mask, vl);
                }
            };
            if ((interpolation & ~CV_HAL_WARP_RELATIVE_MAP) == CV_HAL_INTER_NEAREST)
            {
                auto ix = __riscv_vfcvt_x(mx, vl), iy = __riscv_vfcvt_x(my, vl);
                vuint8mf2_t src0, src1, src2;
                access(ix, iy, src0, src1, src2);
                vuint8mf2x3_t dst{};
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 0, src0);
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 1, src1);
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 2, src2);
                __riscv_vsseg3e8(dst_data + i * dst_step + j * 3, dst, vl);
            }
            else
            {
                vint32m2_t ix0, iy0;
                int rd;
                asm volatile("fsrmi %0, 2 \n\t vsetvli zero,%3,e32,m2,ta,ma \n\t vfcvt.x.f.v %1,%4 \n\t vfcvt.x.f.v %2,%5 \n\t fsrm %0"
                             : "=&r"(rd), "=vr"(ix0), "=vr"(iy0)
                             : "r"(vl), "vr"(mx), "vr"(my)); // Rounding Mode: RDN
                auto ix1 = __riscv_vadd(ix0, 1, vl), iy1 = __riscv_vadd(iy0, 1, vl);

                vfloat32m2_t v00, v10, v20;
                vfloat32m2_t v01, v11, v21;
                vfloat32m2_t v02, v12, v22;
                vfloat32m2_t v03, v13, v23;
                vuint8mf2_t src0, src1, src2;
                access(ix0, iy0, src0, src1, src2);
                v00 = __riscv_vfcvt_f(__riscv_vzext_vf4(src0, vl), vl);
                v10 = __riscv_vfcvt_f(__riscv_vzext_vf4(src1, vl), vl);
                v20 = __riscv_vfcvt_f(__riscv_vzext_vf4(src2, vl), vl);
                access(ix1, iy0, src0, src1, src2);
                v01 = __riscv_vfcvt_f(__riscv_vzext_vf4(src0, vl), vl);
                v11 = __riscv_vfcvt_f(__riscv_vzext_vf4(src1, vl), vl);
                v21 = __riscv_vfcvt_f(__riscv_vzext_vf4(src2, vl), vl);
                access(ix0, iy1, src0, src1, src2);
                v02 = __riscv_vfcvt_f(__riscv_vzext_vf4(src0, vl), vl);
                v12 = __riscv_vfcvt_f(__riscv_vzext_vf4(src1, vl), vl);
                v22 = __riscv_vfcvt_f(__riscv_vzext_vf4(src2, vl), vl);
                access(ix1, iy1, src0, src1, src2);
                v03 = __riscv_vfcvt_f(__riscv_vzext_vf4(src0, vl), vl);
                v13 = __riscv_vfcvt_f(__riscv_vzext_vf4(src1, vl), vl);
                v23 = __riscv_vfcvt_f(__riscv_vzext_vf4(src2, vl), vl);

                mx = __riscv_vfsub(mx, __riscv_vfcvt_f(ix0, vl), vl);
                my = __riscv_vfsub(my, __riscv_vfcvt_f(iy0, vl), vl);
                v00 = __riscv_vfmacc(v00, mx, __riscv_vfsub(v01, v00, vl), vl);
                v02 = __riscv_vfmacc(v02, mx, __riscv_vfsub(v03, v02, vl), vl);
                v00 = __riscv_vfmacc(v00, my, __riscv_vfsub(v02, v00, vl), vl);
                v10 = __riscv_vfmacc(v10, mx, __riscv_vfsub(v11, v10, vl), vl);
                v12 = __riscv_vfmacc(v12, mx, __riscv_vfsub(v13, v12, vl), vl);
                v10 = __riscv_vfmacc(v10, my, __riscv_vfsub(v12, v10, vl), vl);
                v20 = __riscv_vfmacc(v20, mx, __riscv_vfsub(v21, v20, vl), vl);
                v22 = __riscv_vfmacc(v22, mx, __riscv_vfsub(v23, v22, vl), vl);
                v20 = __riscv_vfmacc(v20, my, __riscv_vfsub(v22, v20, vl), vl);
                vuint8mf2x3_t dst{};
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 0, __riscv_vnclipu(__riscv_vfncvt_xu(v00, vl), 0, __RISCV_VXRM_RNU, vl));
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 1, __riscv_vnclipu(__riscv_vfncvt_xu(v10, vl), 0, __RISCV_VXRM_RNU, vl));
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 2, __riscv_vnclipu(__riscv_vfncvt_xu(v20, vl), 0, __RISCV_VXRM_RNU, vl));
                __riscv_vsseg3e8(dst_data + i * dst_step + j * 3, dst, vl);
            }
        }
    }

    return CV_HAL_ERROR_OK;
}

static inline int remap32fC4(int start, int end, const uchar *src_data, size_t src_step, int src_width, int src_height,
                             uchar *dst_data, size_t dst_step, int dst_width,
                             const float* mapx, size_t mapx_step, const float* mapy, size_t mapy_step,
                             int interpolation, int border_type, const double* border_value)
{
    for (int i = start; i < end; i++)
    {
        int vl;
        for (int j = 0; j < dst_width; j += vl)
        {
            vl = __riscv_vsetvl_e8mf2(dst_width - j);
            auto mx = __riscv_vle32_v_f32m2(mapx + i * mapx_step + j, vl);
            auto my = __riscv_vle32_v_f32m2(mapy + i * mapy_step + j, vl);
            if (interpolation & CV_HAL_WARP_RELATIVE_MAP)
            {
                mx = __riscv_vfadd(mx, __riscv_vfcvt_f(__riscv_vadd(__riscv_vid_v_u32m2(vl), j, vl), vl), vl);
                my = __riscv_vfadd(my, i, vl);
            }

            auto access = [&](vint32m2_t ix, vint32m2_t iy, vuint8mf2_t& src0, vuint8mf2_t& src1, vuint8mf2_t& src2, vuint8mf2_t& src3) {
                auto ux = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vmin(__riscv_vmax(ix, 0, vl), src_width  - 1, vl));
                auto uy = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vmin(__riscv_vmax(iy, 0, vl), src_height - 1, vl));
                auto src = __riscv_vloxseg4ei32_v_u8mf2x4(src_data, __riscv_vmadd(uy, src_step, __riscv_vmul(ux, 4, vl), vl), vl);
                src0 = __riscv_vget_v_u8mf2x4_u8mf2(src, 0);
                src1 = __riscv_vget_v_u8mf2x4_u8mf2(src, 1);
                src2 = __riscv_vget_v_u8mf2x4_u8mf2(src, 2);
                src3 = __riscv_vget_v_u8mf2x4_u8mf2(src, 3);
                if (border_type == CV_HAL_BORDER_CONSTANT)
                {
                    auto mask = __riscv_vmor(__riscv_vmsne(ix, __riscv_vreinterpret_v_u32m2_i32m2(ux), vl), __riscv_vmsne(iy, __riscv_vreinterpret_v_u32m2_i32m2(uy), vl), vl);
                    src0 = __riscv_vmerge(src0, border_value[0], mask, vl);
                    src1 = __riscv_vmerge(src1, border_value[1], mask, vl);
                    src2 = __riscv_vmerge(src2, border_value[2], mask, vl);
                    src3 = __riscv_vmerge(src3, border_value[3], mask, vl);
                }
            };
            if ((interpolation & ~CV_HAL_WARP_RELATIVE_MAP) == CV_HAL_INTER_NEAREST)
            {
                auto ix = __riscv_vfcvt_x(mx, vl), iy = __riscv_vfcvt_x(my, vl);
                vuint8mf2_t src0, src1, src2, src3;
                access(ix, iy, src0, src1, src2, src3);
                vuint8mf2x4_t dst{};
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 0, src0);
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 1, src1);
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 2, src2);
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 3, src3);
                __riscv_vsseg4e8(dst_data + i * dst_step + j * 4, dst, vl);
            }
            else
            {
                vint32m2_t ix0, iy0;
                int rd;
                asm volatile("fsrmi %0, 2 \n\t vsetvli zero,%3,e32,m2,ta,ma \n\t vfcvt.x.f.v %1,%4 \n\t vfcvt.x.f.v %2,%5 \n\t fsrm %0"
                             : "=&r"(rd), "=vr"(ix0), "=vr"(iy0)
                             : "r"(vl), "vr"(mx), "vr"(my)); // Rounding Mode: RDN
                auto ix1 = __riscv_vadd(ix0, 1, vl), iy1 = __riscv_vadd(iy0, 1, vl);

                vfloat32m2_t v00, v10, v20, v30;
                vfloat32m2_t v01, v11, v21, v31;
                vfloat32m2_t v02, v12, v22, v32;
                vfloat32m2_t v03, v13, v23, v33;
                vuint8mf2_t src0, src1, src2, src3;
                access(ix0, iy0, src0, src1, src2, src3);
                v00 = __riscv_vfcvt_f(__riscv_vzext_vf4(src0, vl), vl);
                v10 = __riscv_vfcvt_f(__riscv_vzext_vf4(src1, vl), vl);
                v20 = __riscv_vfcvt_f(__riscv_vzext_vf4(src2, vl), vl);
                v30 = __riscv_vfcvt_f(__riscv_vzext_vf4(src3, vl), vl);
                access(ix1, iy0, src0, src1, src2, src3);
                v01 = __riscv_vfcvt_f(__riscv_vzext_vf4(src0, vl), vl);
                v11 = __riscv_vfcvt_f(__riscv_vzext_vf4(src1, vl), vl);
                v21 = __riscv_vfcvt_f(__riscv_vzext_vf4(src2, vl), vl);
                v31 = __riscv_vfcvt_f(__riscv_vzext_vf4(src3, vl), vl);
                access(ix0, iy1, src0, src1, src2, src3);
                v02 = __riscv_vfcvt_f(__riscv_vzext_vf4(src0, vl), vl);
                v12 = __riscv_vfcvt_f(__riscv_vzext_vf4(src1, vl), vl);
                v22 = __riscv_vfcvt_f(__riscv_vzext_vf4(src2, vl), vl);
                v32 = __riscv_vfcvt_f(__riscv_vzext_vf4(src3, vl), vl);
                access(ix1, iy1, src0, src1, src2, src3);
                v03 = __riscv_vfcvt_f(__riscv_vzext_vf4(src0, vl), vl);
                v13 = __riscv_vfcvt_f(__riscv_vzext_vf4(src1, vl), vl);
                v23 = __riscv_vfcvt_f(__riscv_vzext_vf4(src2, vl), vl);
                v33 = __riscv_vfcvt_f(__riscv_vzext_vf4(src3, vl), vl);

                mx = __riscv_vfsub(mx, __riscv_vfcvt_f(ix0, vl), vl);
                my = __riscv_vfsub(my, __riscv_vfcvt_f(iy0, vl), vl);
                v00 = __riscv_vfmacc(v00, mx, __riscv_vfsub(v01, v00, vl), vl);
                v02 = __riscv_vfmacc(v02, mx, __riscv_vfsub(v03, v02, vl), vl);
                v00 = __riscv_vfmacc(v00, my, __riscv_vfsub(v02, v00, vl), vl);
                v10 = __riscv_vfmacc(v10, mx, __riscv_vfsub(v11, v10, vl), vl);
                v12 = __riscv_vfmacc(v12, mx, __riscv_vfsub(v13, v12, vl), vl);
                v10 = __riscv_vfmacc(v10, my, __riscv_vfsub(v12, v10, vl), vl);
                v20 = __riscv_vfmacc(v20, mx, __riscv_vfsub(v21, v20, vl), vl);
                v22 = __riscv_vfmacc(v22, mx, __riscv_vfsub(v23, v22, vl), vl);
                v20 = __riscv_vfmacc(v20, my, __riscv_vfsub(v22, v20, vl), vl);
                v30 = __riscv_vfmacc(v30, mx, __riscv_vfsub(v31, v30, vl), vl);
                v32 = __riscv_vfmacc(v32, mx, __riscv_vfsub(v33, v32, vl), vl);
                v30 = __riscv_vfmacc(v30, my, __riscv_vfsub(v32, v30, vl), vl);
                vuint8mf2x4_t dst{};
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 0, __riscv_vnclipu(__riscv_vfncvt_xu(v00, vl), 0, __RISCV_VXRM_RNU, vl));
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 1, __riscv_vnclipu(__riscv_vfncvt_xu(v10, vl), 0, __RISCV_VXRM_RNU, vl));
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 2, __riscv_vnclipu(__riscv_vfncvt_xu(v20, vl), 0, __RISCV_VXRM_RNU, vl));
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 3, __riscv_vnclipu(__riscv_vfncvt_xu(v30, vl), 0, __RISCV_VXRM_RNU, vl));
                __riscv_vsseg4e8(dst_data + i * dst_step + j * 4, dst, vl);
            }
        }
    }

    return CV_HAL_ERROR_OK;
}

// the algorithm is copied from 3rdparty/carotene/src/remap.cpp,
// in the function void CAROTENE_NS::remapNearestNeighbor and void CAROTENE_NS::remapLinear
inline int remap32f(int src_type, const uchar *src_data, size_t src_step, int src_width, int src_height,
                    uchar *dst_data, size_t dst_step, int dst_width, int dst_height,
                    float* mapx, size_t mapx_step, float* mapy, size_t mapy_step,
                    int interpolation, int border_type, const double border_value[4])
{
    if (src_type != CV_8UC1 && src_type != CV_8UC3 && src_type != CV_8UC4 && src_type != CV_16UC1 && src_type != CV_16SC1 && src_type != CV_32FC1)
        return CV_HAL_ERROR_NOT_IMPLEMENTED;
    if (border_type != CV_HAL_BORDER_CONSTANT && border_type != CV_HAL_BORDER_REPLICATE)
        return CV_HAL_ERROR_NOT_IMPLEMENTED;
    if ((interpolation & ~CV_HAL_WARP_RELATIVE_MAP) != CV_HAL_INTER_NEAREST && (interpolation & ~CV_HAL_WARP_RELATIVE_MAP) != CV_HAL_INTER_LINEAR)
        return CV_HAL_ERROR_NOT_IMPLEMENTED;

    mapx_step /= sizeof(float);
    mapy_step /= sizeof(float);
    switch (src_type)
    {
    case CV_8UC1:
        return invoke(dst_width, dst_height, {remap32fC1<RVV_U8M2>}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, mapx, mapx_step, mapy, mapy_step, interpolation, border_type, border_value);
    case CV_16UC1:
        return invoke(dst_width, dst_height, {remap32fC1<RVV_U16M4>}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, mapx, mapx_step, mapy, mapy_step, interpolation, border_type, border_value);
    case CV_16SC1:
        return invoke(dst_width, dst_height, {remap32fC1<RVV_I16M4>}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, mapx, mapx_step, mapy, mapy_step, interpolation, border_type, border_value);
    case CV_32FC1:
        return invoke(dst_width, dst_height, {remap32fC1<RVV_F32M8>}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, mapx, mapx_step, mapy, mapy_step, interpolation, border_type, border_value);
    case CV_8UC3:
        return invoke(dst_width, dst_height, {remap32fC3}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, mapx, mapx_step, mapy, mapy_step, interpolation, border_type, border_value);
    case CV_8UC4:
        return invoke(dst_width, dst_height, {remap32fC4}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, mapx, mapx_step, mapy, mapy_step, interpolation, border_type, border_value);
    }

    return CV_HAL_ERROR_NOT_IMPLEMENTED;
}
} // cv::cv_hal_rvv::remap

namespace warp {
#undef cv_hal_warpAffine
#define cv_hal_warpAffine cv::cv_hal_rvv::warp::warpAffine
#undef cv_hal_warpPerspective
#define cv_hal_warpPerspective cv::cv_hal_rvv::warp::warpPerspective

template<bool perspective>
static inline int warpC1(int start, int end, const uchar *src_data, size_t src_step, int src_width, int src_height, uchar *dst_data, size_t dst_step, int dst_width, const double* M, int interpolation, int borderType, const double* borderValue)
{
    for (int i = start; i < end; i++)
    {
        int vl;
        for (int j = 0; j < dst_width; j += vl)
        {
            vl = __riscv_vsetvl_e8m1(dst_width - j);
            auto access = [&](vint32m4_t ix, vint32m4_t iy) {
                auto ux = __riscv_vreinterpret_v_i32m4_u32m4(__riscv_vmin(__riscv_vmax(ix, 0, vl), src_width  - 1, vl));
                auto uy = __riscv_vreinterpret_v_i32m4_u32m4(__riscv_vmin(__riscv_vmax(iy, 0, vl), src_height - 1, vl));
                auto src = __riscv_vloxei32_v_u8m1(src_data, __riscv_vmadd(uy, src_step, ux, vl), vl);
                if (borderType == CV_HAL_BORDER_CONSTANT)
                {
                    auto mask = __riscv_vmor(__riscv_vmsne(ix, __riscv_vreinterpret_v_u32m4_i32m4(ux), vl), __riscv_vmsne(iy, __riscv_vreinterpret_v_u32m4_i32m4(uy), vl), vl);
                    src = __riscv_vmerge(src, borderValue[0], mask, vl);
                }
                return src;
            };

            auto id = __riscv_vfcvt_f(__riscv_vadd(__riscv_vid_v_u32m4(vl), j, vl), vl);
            auto mx = __riscv_vfmadd(id, M[0], __riscv_vfmadd(__riscv_vfmv_v_f_f32m4(i, vl), M[1], __riscv_vfmv_v_f_f32m4(M[2], vl), vl), vl);
            auto my = __riscv_vfmadd(id, M[3], __riscv_vfmadd(__riscv_vfmv_v_f_f32m4(i, vl), M[4], __riscv_vfmv_v_f_f32m4(M[5], vl), vl), vl);
            if (perspective)
            {
                auto md = __riscv_vfrdiv(__riscv_vfmadd(id, M[6], __riscv_vfmadd(__riscv_vfmv_v_f_f32m4(i, vl), M[7], __riscv_vfmv_v_f_f32m4(M[8], vl), vl), vl), 1, vl);
                mx = __riscv_vfmul(mx, md, vl);
                my = __riscv_vfmul(my, md, vl);
            }

            if (interpolation == CV_HAL_INTER_NEAREST)
            {
                auto ix = __riscv_vfcvt_x(mx, vl), iy = __riscv_vfcvt_x(my, vl);
                __riscv_vse8(dst_data + i * dst_step + j, access(ix, iy), vl);
            }
            else
            {
                auto ix = __riscv_vfcvt_x(__riscv_vfmadd(mx, 1 << 10, __riscv_vfmv_v_f_f32m4(1 << 4, vl), vl), vl);
                auto iy = __riscv_vfcvt_x(__riscv_vfmadd(my, 1 << 10, __riscv_vfmv_v_f_f32m4(1 << 4, vl), vl), vl);
                auto ix0 = __riscv_vsra(ix, 10, vl), iy0 = __riscv_vsra(iy, 10, vl);
                auto ix1 = __riscv_vadd(ix0, 1, vl), iy1 = __riscv_vadd(iy0, 1, vl);

                auto v0 = __riscv_vzext_vf4(access(ix0, iy0), vl);
                auto v1 = __riscv_vzext_vf4(access(ix1, iy0), vl);
                auto v2 = __riscv_vzext_vf4(access(ix0, iy1), vl);
                auto v3 = __riscv_vzext_vf4(access(ix1, iy1), vl);

                auto rx = __riscv_vreinterpret_v_i32m4_u32m4(__riscv_vand(__riscv_vsra(ix, 5, vl), (1 << 5) - 1, vl));
                auto ry = __riscv_vreinterpret_v_i32m4_u32m4(__riscv_vand(__riscv_vsra(iy, 5, vl), (1 << 5) - 1, vl));
                v0 = __riscv_vmacc(__riscv_vmul(v0, 1 << 5, vl), rx, __riscv_vsub(v1, v0, vl), vl);
                v2 = __riscv_vmacc(__riscv_vmul(v2, 1 << 5, vl), rx, __riscv_vsub(v3, v2, vl), vl);
                v0 = __riscv_vmacc(__riscv_vmul(v0, 1 << 5, vl), ry, __riscv_vsub(v2, v0, vl), vl);
                __riscv_vse8(dst_data + i * dst_step + j, __riscv_vnclipu(__riscv_vnclipu(v0, 10, __RISCV_VXRM_RNU, vl), 0, __RISCV_VXRM_RNU, vl), vl);
            }
        }
    }

    return CV_HAL_ERROR_OK;
}

template<bool perspective>
static inline int warpC3(int start, int end, const uchar *src_data, size_t src_step, int src_width, int src_height, uchar *dst_data, size_t dst_step, int dst_width, const double* M, int interpolation, int borderType, const double* borderValue)
{
    for (int i = start; i < end; i++)
    {
        int vl;
        for (int j = 0; j < dst_width; j += vl)
        {
            vl = __riscv_vsetvl_e8mf2(dst_width - j);
            auto access = [&](vint32m2_t ix, vint32m2_t iy, vuint8mf2_t& src0, vuint8mf2_t& src1, vuint8mf2_t& src2) {
                auto ux = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vmin(__riscv_vmax(ix, 0, vl), src_width  - 1, vl));
                auto uy = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vmin(__riscv_vmax(iy, 0, vl), src_height - 1, vl));
                auto src = __riscv_vloxseg3ei32_v_u8mf2x3(src_data, __riscv_vmadd(uy, src_step, __riscv_vmul(ux, 3, vl), vl), vl);
                src0 = __riscv_vget_v_u8mf2x3_u8mf2(src, 0);
                src1 = __riscv_vget_v_u8mf2x3_u8mf2(src, 1);
                src2 = __riscv_vget_v_u8mf2x3_u8mf2(src, 2);
                if (borderType == CV_HAL_BORDER_CONSTANT)
                {
                    auto mask = __riscv_vmor(__riscv_vmsne(ix, __riscv_vreinterpret_v_u32m2_i32m2(ux), vl), __riscv_vmsne(iy, __riscv_vreinterpret_v_u32m2_i32m2(uy), vl), vl);
                    src0 = __riscv_vmerge(src0, borderValue[0], mask, vl);
                    src1 = __riscv_vmerge(src1, borderValue[1], mask, vl);
                    src2 = __riscv_vmerge(src2, borderValue[2], mask, vl);
                }
            };

            auto id = __riscv_vfcvt_f(__riscv_vadd(__riscv_vid_v_u32m2(vl), j, vl), vl);
            auto mx = __riscv_vfmadd(id, M[0], __riscv_vfmadd(__riscv_vfmv_v_f_f32m2(i, vl), M[1], __riscv_vfmv_v_f_f32m2(M[2], vl), vl), vl);
            auto my = __riscv_vfmadd(id, M[3], __riscv_vfmadd(__riscv_vfmv_v_f_f32m2(i, vl), M[4], __riscv_vfmv_v_f_f32m2(M[5], vl), vl), vl);
            if (perspective)
            {
                auto md = __riscv_vfrdiv(__riscv_vfmadd(id, M[6], __riscv_vfmadd(__riscv_vfmv_v_f_f32m2(i, vl), M[7], __riscv_vfmv_v_f_f32m2(M[8], vl), vl), vl), 1, vl);
                mx = __riscv_vfmul(mx, md, vl);
                my = __riscv_vfmul(my, md, vl);
            }

            if (interpolation == CV_HAL_INTER_NEAREST)
            {
                auto ix = __riscv_vfcvt_x(mx, vl), iy = __riscv_vfcvt_x(my, vl);
                vuint8mf2_t src0, src1, src2;
                access(ix, iy, src0, src1, src2);
                vuint8mf2x3_t dst{};
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 0, src0);
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 1, src1);
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 2, src2);
                __riscv_vsseg3e8(dst_data + i * dst_step + j * 3, dst, vl);
            }
            else
            {
                auto ix = __riscv_vfcvt_x(__riscv_vfmadd(mx, 1 << 10, __riscv_vfmv_v_f_f32m2(1 << 4, vl), vl), vl);
                auto iy = __riscv_vfcvt_x(__riscv_vfmadd(my, 1 << 10, __riscv_vfmv_v_f_f32m2(1 << 4, vl), vl), vl);
                auto ix0 = __riscv_vsra(ix, 10, vl), iy0 = __riscv_vsra(iy, 10, vl);
                auto ix1 = __riscv_vadd(ix0, 1, vl), iy1 = __riscv_vadd(iy0, 1, vl);

                vuint32m2_t v00, v10, v20;
                vuint32m2_t v01, v11, v21;
                vuint32m2_t v02, v12, v22;
                vuint32m2_t v03, v13, v23;
                vuint8mf2_t src0, src1, src2;
                access(ix0, iy0, src0, src1, src2);
                v00 = __riscv_vzext_vf4(src0, vl);
                v10 = __riscv_vzext_vf4(src1, vl);
                v20 = __riscv_vzext_vf4(src2, vl);
                access(ix1, iy0, src0, src1, src2);
                v01 = __riscv_vzext_vf4(src0, vl);
                v11 = __riscv_vzext_vf4(src1, vl);
                v21 = __riscv_vzext_vf4(src2, vl);
                access(ix0, iy1, src0, src1, src2);
                v02 = __riscv_vzext_vf4(src0, vl);
                v12 = __riscv_vzext_vf4(src1, vl);
                v22 = __riscv_vzext_vf4(src2, vl);
                access(ix1, iy1, src0, src1, src2);
                v03 = __riscv_vzext_vf4(src0, vl);
                v13 = __riscv_vzext_vf4(src1, vl);
                v23 = __riscv_vzext_vf4(src2, vl);

                auto rx = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vand(__riscv_vsra(ix, 5, vl), (1 << 5) - 1, vl));
                auto ry = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vand(__riscv_vsra(iy, 5, vl), (1 << 5) - 1, vl));
                v00 = __riscv_vmacc(__riscv_vmul(v00, 1 << 5, vl), rx, __riscv_vsub(v01, v00, vl), vl);
                v02 = __riscv_vmacc(__riscv_vmul(v02, 1 << 5, vl), rx, __riscv_vsub(v03, v02, vl), vl);
                v00 = __riscv_vmacc(__riscv_vmul(v00, 1 << 5, vl), ry, __riscv_vsub(v02, v00, vl), vl);
                v10 = __riscv_vmacc(__riscv_vmul(v10, 1 << 5, vl), rx, __riscv_vsub(v11, v10, vl), vl);
                v12 = __riscv_vmacc(__riscv_vmul(v12, 1 << 5, vl), rx, __riscv_vsub(v13, v12, vl), vl);
                v10 = __riscv_vmacc(__riscv_vmul(v10, 1 << 5, vl), ry, __riscv_vsub(v12, v10, vl), vl);
                v20 = __riscv_vmacc(__riscv_vmul(v20, 1 << 5, vl), rx, __riscv_vsub(v21, v20, vl), vl);
                v22 = __riscv_vmacc(__riscv_vmul(v22, 1 << 5, vl), rx, __riscv_vsub(v23, v22, vl), vl);
                v20 = __riscv_vmacc(__riscv_vmul(v20, 1 << 5, vl), ry, __riscv_vsub(v22, v20, vl), vl);
                vuint8mf2x3_t dst{};
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 0, __riscv_vnclipu(__riscv_vnclipu(v00, 10, __RISCV_VXRM_RNU, vl), 0, __RISCV_VXRM_RNU, vl));
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 1, __riscv_vnclipu(__riscv_vnclipu(v10, 10, __RISCV_VXRM_RNU, vl), 0, __RISCV_VXRM_RNU, vl));
                dst = __riscv_vset_v_u8mf2_u8mf2x3(dst, 2, __riscv_vnclipu(__riscv_vnclipu(v20, 10, __RISCV_VXRM_RNU, vl), 0, __RISCV_VXRM_RNU, vl));
                __riscv_vsseg3e8(dst_data + i * dst_step + j * 3, dst, vl);
            }
        }
    }

    return CV_HAL_ERROR_OK;
}

template<bool perspective>
static inline int warpC4(int start, int end, const uchar *src_data, size_t src_step, int src_width, int src_height, uchar *dst_data, size_t dst_step, int dst_width, const double* M, int interpolation, int borderType, const double* borderValue)
{
    for (int i = start; i < end; i++)
    {
        int vl;
        for (int j = 0; j < dst_width; j += vl)
        {
            vl = __riscv_vsetvl_e8mf2(dst_width - j);
            auto access = [&](vint32m2_t ix, vint32m2_t iy, vuint8mf2_t& src0, vuint8mf2_t& src1, vuint8mf2_t& src2, vuint8mf2_t& src3) {
                auto ux = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vmin(__riscv_vmax(ix, 0, vl), src_width  - 1, vl));
                auto uy = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vmin(__riscv_vmax(iy, 0, vl), src_height - 1, vl));
                auto src = __riscv_vloxseg4ei32_v_u8mf2x4(src_data, __riscv_vmadd(uy, src_step, __riscv_vmul(ux, 4, vl), vl), vl);
                src0 = __riscv_vget_v_u8mf2x4_u8mf2(src, 0);
                src1 = __riscv_vget_v_u8mf2x4_u8mf2(src, 1);
                src2 = __riscv_vget_v_u8mf2x4_u8mf2(src, 2);
                src3 = __riscv_vget_v_u8mf2x4_u8mf2(src, 3);
                if (borderType == CV_HAL_BORDER_CONSTANT)
                {
                    auto mask = __riscv_vmor(__riscv_vmsne(ix, __riscv_vreinterpret_v_u32m2_i32m2(ux), vl), __riscv_vmsne(iy, __riscv_vreinterpret_v_u32m2_i32m2(uy), vl), vl);
                    src0 = __riscv_vmerge(src0, borderValue[0], mask, vl);
                    src1 = __riscv_vmerge(src1, borderValue[1], mask, vl);
                    src2 = __riscv_vmerge(src2, borderValue[2], mask, vl);
                    src3 = __riscv_vmerge(src3, borderValue[3], mask, vl);
                }
            };

            auto id = __riscv_vfcvt_f(__riscv_vadd(__riscv_vid_v_u32m2(vl), j, vl), vl);
            auto mx = __riscv_vfmadd(id, M[0], __riscv_vfmadd(__riscv_vfmv_v_f_f32m2(i, vl), M[1], __riscv_vfmv_v_f_f32m2(M[2], vl), vl), vl);
            auto my = __riscv_vfmadd(id, M[3], __riscv_vfmadd(__riscv_vfmv_v_f_f32m2(i, vl), M[4], __riscv_vfmv_v_f_f32m2(M[5], vl), vl), vl);
            if (perspective)
            {
                auto md = __riscv_vfrdiv(__riscv_vfmadd(id, M[6], __riscv_vfmadd(__riscv_vfmv_v_f_f32m2(i, vl), M[7], __riscv_vfmv_v_f_f32m2(M[8], vl), vl), vl), 1, vl);
                mx = __riscv_vfmul(mx, md, vl);
                my = __riscv_vfmul(my, md, vl);
            }

            if (interpolation == CV_HAL_INTER_NEAREST)
            {
                auto ix = __riscv_vfcvt_x(mx, vl), iy = __riscv_vfcvt_x(my, vl);
                vuint8mf2_t src0, src1, src2, src3;
                access(ix, iy, src0, src1, src2, src3);
                vuint8mf2x4_t dst{};
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 0, src0);
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 1, src1);
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 2, src2);
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 3, src3);
                __riscv_vsseg4e8(dst_data + i * dst_step + j * 4, dst, vl);
            }
            else
            {
                auto ix = __riscv_vfcvt_x(__riscv_vfmadd(mx, 1 << 10, __riscv_vfmv_v_f_f32m2(1 << 4, vl), vl), vl);
                auto iy = __riscv_vfcvt_x(__riscv_vfmadd(my, 1 << 10, __riscv_vfmv_v_f_f32m2(1 << 4, vl), vl), vl);
                auto ix0 = __riscv_vsra(ix, 10, vl), iy0 = __riscv_vsra(iy, 10, vl);
                auto ix1 = __riscv_vadd(ix0, 1, vl), iy1 = __riscv_vadd(iy0, 1, vl);

                vuint32m2_t v00, v10, v20, v30;
                vuint32m2_t v01, v11, v21, v31;
                vuint32m2_t v02, v12, v22, v32;
                vuint32m2_t v03, v13, v23, v33;
                vuint8mf2_t src0, src1, src2, src3;
                access(ix0, iy0, src0, src1, src2, src3);
                v00 = __riscv_vzext_vf4(src0, vl);
                v10 = __riscv_vzext_vf4(src1, vl);
                v20 = __riscv_vzext_vf4(src2, vl);
                v30 = __riscv_vzext_vf4(src3, vl);
                access(ix1, iy0, src0, src1, src2, src3);
                v01 = __riscv_vzext_vf4(src0, vl);
                v11 = __riscv_vzext_vf4(src1, vl);
                v21 = __riscv_vzext_vf4(src2, vl);
                v31 = __riscv_vzext_vf4(src3, vl);
                access(ix0, iy1, src0, src1, src2, src3);
                v02 = __riscv_vzext_vf4(src0, vl);
                v12 = __riscv_vzext_vf4(src1, vl);
                v22 = __riscv_vzext_vf4(src2, vl);
                v32 = __riscv_vzext_vf4(src3, vl);
                access(ix1, iy1, src0, src1, src2, src3);
                v03 = __riscv_vzext_vf4(src0, vl);
                v13 = __riscv_vzext_vf4(src1, vl);
                v23 = __riscv_vzext_vf4(src2, vl);
                v33 = __riscv_vzext_vf4(src3, vl);

                auto rx = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vand(__riscv_vsra(ix, 5, vl), (1 << 5) - 1, vl));
                auto ry = __riscv_vreinterpret_v_i32m2_u32m2(__riscv_vand(__riscv_vsra(iy, 5, vl), (1 << 5) - 1, vl));
                v00 = __riscv_vmacc(__riscv_vmul(v00, 1 << 5, vl), rx, __riscv_vsub(v01, v00, vl), vl);
                v02 = __riscv_vmacc(__riscv_vmul(v02, 1 << 5, vl), rx, __riscv_vsub(v03, v02, vl), vl);
                v00 = __riscv_vmacc(__riscv_vmul(v00, 1 << 5, vl), ry, __riscv_vsub(v02, v00, vl), vl);
                v10 = __riscv_vmacc(__riscv_vmul(v10, 1 << 5, vl), rx, __riscv_vsub(v11, v10, vl), vl);
                v12 = __riscv_vmacc(__riscv_vmul(v12, 1 << 5, vl), rx, __riscv_vsub(v13, v12, vl), vl);
                v10 = __riscv_vmacc(__riscv_vmul(v10, 1 << 5, vl), ry, __riscv_vsub(v12, v10, vl), vl);
                v20 = __riscv_vmacc(__riscv_vmul(v20, 1 << 5, vl), rx, __riscv_vsub(v21, v20, vl), vl);
                v22 = __riscv_vmacc(__riscv_vmul(v22, 1 << 5, vl), rx, __riscv_vsub(v23, v22, vl), vl);
                v20 = __riscv_vmacc(__riscv_vmul(v20, 1 << 5, vl), ry, __riscv_vsub(v22, v20, vl), vl);
                v30 = __riscv_vmacc(__riscv_vmul(v30, 1 << 5, vl), rx, __riscv_vsub(v31, v30, vl), vl);
                v32 = __riscv_vmacc(__riscv_vmul(v32, 1 << 5, vl), rx, __riscv_vsub(v33, v32, vl), vl);
                v30 = __riscv_vmacc(__riscv_vmul(v30, 1 << 5, vl), ry, __riscv_vsub(v32, v30, vl), vl);
                vuint8mf2x4_t dst{};
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 0, __riscv_vnclipu(__riscv_vnclipu(v00, 10, __RISCV_VXRM_RNU, vl), 0, __RISCV_VXRM_RNU, vl));
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 1, __riscv_vnclipu(__riscv_vnclipu(v10, 10, __RISCV_VXRM_RNU, vl), 0, __RISCV_VXRM_RNU, vl));
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 2, __riscv_vnclipu(__riscv_vnclipu(v20, 10, __RISCV_VXRM_RNU, vl), 0, __RISCV_VXRM_RNU, vl));
                dst = __riscv_vset_v_u8mf2_u8mf2x4(dst, 3, __riscv_vnclipu(__riscv_vnclipu(v30, 10, __RISCV_VXRM_RNU, vl), 0, __RISCV_VXRM_RNU, vl));
                __riscv_vsseg4e8(dst_data + i * dst_step + j * 4, dst, vl);
            }
        }
    }

    return CV_HAL_ERROR_OK;
}

// the algorithm is copied from 3rdparty/carotene/src/warp_affine.cpp,
// in the function void CAROTENE_NS::warpAffineNearestNeighbor and void CAROTENE_NS::warpAffineLinear
inline int warpAffine(int src_type, const uchar *src_data, size_t src_step, int src_width, int src_height, uchar *dst_data, size_t dst_step, int dst_width, int dst_height, const double M[6], int interpolation, int borderType, const double borderValue[4])
{
    if (src_type != CV_8UC1 && src_type != CV_8UC3 && src_type != CV_8UC4)
        return CV_HAL_ERROR_NOT_IMPLEMENTED;
    if (borderType != CV_HAL_BORDER_CONSTANT && borderType != CV_HAL_BORDER_REPLICATE)
        return CV_HAL_ERROR_NOT_IMPLEMENTED;
    if (interpolation != CV_HAL_INTER_NEAREST && interpolation != CV_HAL_INTER_LINEAR)
        return CV_HAL_ERROR_NOT_IMPLEMENTED;

    switch (src_type)
    {
    case CV_8UC1:
        return remap::invoke(dst_width, dst_height, {warpC1<false>}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, M, interpolation, borderType, borderValue);
    case CV_8UC3:
        return remap::invoke(dst_width, dst_height, {warpC3<false>}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, M, interpolation, borderType, borderValue);
    case CV_8UC4:
        return remap::invoke(dst_width, dst_height, {warpC4<false>}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, M, interpolation, borderType, borderValue);
    }

    return CV_HAL_ERROR_NOT_IMPLEMENTED;
}

// the algorithm is copied from 3rdparty/carotene/src/warp_perspective.cpp,
// in the function void CAROTENE_NS::warpPerspectiveNearestNeighbor and void CAROTENE_NS::warpPerspectiveLinear
inline int warpPerspective(int src_type, const uchar *src_data, size_t src_step, int src_width, int src_height, uchar *dst_data, size_t dst_step, int dst_width, int dst_height, const double M[9], int interpolation, int borderType, const double borderValue[4])
{
    if (src_type != CV_8UC1 && src_type != CV_8UC3 && src_type != CV_8UC4)
        return CV_HAL_ERROR_NOT_IMPLEMENTED;
    if (borderType != CV_HAL_BORDER_CONSTANT && borderType != CV_HAL_BORDER_REPLICATE)
        return CV_HAL_ERROR_NOT_IMPLEMENTED;
    if (interpolation != CV_HAL_INTER_NEAREST && interpolation != CV_HAL_INTER_LINEAR)
        return CV_HAL_ERROR_NOT_IMPLEMENTED;

    switch (src_type)
    {
    case CV_8UC1:
        return remap::invoke(dst_width, dst_height, {warpC1<true>}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, M, interpolation, borderType, borderValue);
    case CV_8UC3:
        return remap::invoke(dst_width, dst_height, {warpC3<true>}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, M, interpolation, borderType, borderValue);
    case CV_8UC4:
        return remap::invoke(dst_width, dst_height, {warpC4<true>}, src_data, src_step, src_width, src_height, dst_data, dst_step, dst_width, M, interpolation, borderType, borderValue);
    }

    return CV_HAL_ERROR_NOT_IMPLEMENTED;
}
} // cv::cv_hal_rvv::warp

}}

#endif
