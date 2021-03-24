/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                          License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
// Copyright (C) 2013, OpenCV Foundation, all rights reserved.
// Copyright (C) 2015, Itseez Inc., all rights reserved.
// Copyright (C) 2021, PingTouGe Semiconductor Co., Ltd all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

#ifndef OPENCV_HAL_INTRIN_RISCVV_HPP
#define OPENCV_HAL_INTRIN_RISCVV_HPP

#include <float.h>
#include <algorithm>
#include "opencv2/core/utility.hpp"


namespace cv
{

//! @cond IGNORED

CV_CPU_OPTIMIZATION_HAL_NAMESPACE_BEGIN

#define CV_SIMD128 1
#define CV_SIMD128_64F 1
//////////// Types ////////////
struct v_uint8x16
{
    typedef uchar lane_type;
    enum { nlanes = 16 };

    v_uint8x16() {}
    explicit v_uint8x16(vuint8m1_t v) : val(v) {}
    v_uint8x16(uchar v0, uchar v1, uchar v2, uchar v3, uchar v4, uchar v5, uchar v6, uchar v7,
               uchar v8, uchar v9, uchar v10, uchar v11, uchar v12, uchar v13, uchar v14, uchar v15)
    {
        uchar v[] = {v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15};
        val = (vuint8m1_t)vle_v_u8m1((unsigned char*)v, 16);
    }
    uchar get0() const
    {
        return vmv_x_s_u8m1_u8(val, 16);
    }

    vuint8m1_t val;
};

struct v_int8x16
{
    typedef schar lane_type;
    enum { nlanes = 16 };

    v_int8x16() {}
    explicit v_int8x16(vint8m1_t v) : val(v) {}
    v_int8x16(schar v0, schar v1, schar v2, schar v3, schar v4, schar v5, schar v6, schar v7,
               schar v8, schar v9, schar v10, schar v11, schar v12, schar v13, schar v14, schar v15)
    {
        schar v[] = {v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15};
        val = (vint8m1_t)vle_v_i8m1((schar*)v, 16);
    }
    schar get0() const
    {
        return vmv_x_s_i8m1_i8(val, 16);
    }

    vint8m1_t val;
};

struct v_uint16x8
{
    typedef ushort lane_type;
    enum { nlanes = 8 };

    v_uint16x8() {}
    explicit v_uint16x8(vuint16m1_t v) : val(v) {}
    v_uint16x8(ushort v0, ushort v1, ushort v2, ushort v3, ushort v4, ushort v5, ushort v6, ushort v7)
    {
        ushort v[] = {v0, v1, v2, v3, v4, v5, v6, v7};
        val = (vuint16m1_t)vle_v_u16m1((unsigned short*)v, 8);
    }
    ushort get0() const
    {
        return vmv_x_s_u16m1_u16(val, 8);
    }

    vuint16m1_t val;
};

struct v_int16x8
{
    typedef short lane_type;
    enum { nlanes = 8 };

    v_int16x8() {}
    explicit v_int16x8(vint16m1_t v) : val(v) {}
    v_int16x8(short v0, short v1, short v2, short v3, short v4, short v5, short v6, short v7)
    {
        short v[] = {v0, v1, v2, v3, v4, v5, v6, v7};
        val = (vint16m1_t)vle_v_i16m1((signed short*)v, 8);
    }
    short get0() const
    {
        return vmv_x_s_i16m1_i16(val, 8);
    }

    vint16m1_t val;
};

struct v_uint32x4
{
    typedef unsigned lane_type;
    enum { nlanes = 4 };

    v_uint32x4() {}
    explicit v_uint32x4(vuint32m1_t v) : val(v) {}
    v_uint32x4(unsigned v0, unsigned v1, unsigned v2, unsigned v3)
    {
        unsigned v[] = {v0, v1, v2, v3};
        val = (vuint32m1_t)vle_v_u32m1((unsigned int*)v, 4);
    }
    unsigned get0() const
    {
        return vmv_x_s_u32m1_u32(val, 4);
    }

    vuint32m1_t val;
};

struct v_int32x4
{
    typedef int lane_type;
    enum { nlanes = 4 };

    v_int32x4() {}
    explicit v_int32x4(vint32m1_t v) : val(v) {}
    v_int32x4(int v0, int v1, int v2, int v3)
    {
        int v[] = {v0, v1, v2, v3};
        val = (vint32m1_t)vle_v_i32m1((signed int*)v, 4);
    }
    int get0() const
    {
        return vmv_x_s_i32m1_i32(val, 4);
    }
    vint32m1_t val;
};

struct v_float32x4
{
    typedef float lane_type;
    enum { nlanes = 4 };

    v_float32x4() {}
    explicit v_float32x4(vfloat32m1_t v) : val(v) {}
    v_float32x4(float v0, float v1, float v2, float v3)
    {
        float v[] = {v0, v1, v2, v3};
        val = (vfloat32m1_t)vle_v_f32m1((float*)v, 4);
    }
    float get0() const
    {
        return vfmv_f_s_f32m1_f32(val, 4);
    }
    vfloat32m1_t val;
};

struct v_uint64x2
{
    typedef uint64 lane_type;
    enum { nlanes = 2 };

    v_uint64x2() {}
    explicit v_uint64x2(vuint64m1_t v) : val(v) {}
    v_uint64x2(uint64 v0, uint64 v1)
    {
        uint64 v[] = {v0, v1};
        val = (vuint64m1_t)vle_v_u64m1((unsigned long*)v, 2);
    }
    uint64 get0() const
    {
        return vmv_x_s_u64m1_u64(val, 2);
    }
    vuint64m1_t val;
};

struct v_int64x2
{
    typedef int64 lane_type;
    enum { nlanes = 2 };

    v_int64x2() {}
    explicit v_int64x2(vint64m1_t v) : val(v) {}
    v_int64x2(int64 v0, int64 v1)
    {
        int64 v[] = {v0, v1};
        val = (vint64m1_t)vle_v_i64m1((long*)v, 2);
    }
    int64 get0() const
    {
        return vmv_x_s_i64m1_i64(val, 2);
    }
    vint64m1_t val;
};

struct v_float64x2
{
    typedef double lane_type;
    enum { nlanes = 2 };

    v_float64x2() {}
    explicit v_float64x2(vfloat64m1_t v) : val(v) {}
    v_float64x2(double v0, double v1)
    {
        double v[] = {v0, v1};
        val = (vfloat64m1_t)vle_v_f64m1((double*)v, 2);
    }
    double get0() const
    {
        return vfmv_f_s_f64m1_f64(val, 2);
    }
    vfloat64m1_t val;
};

#define OPENCV_HAL_IMPL_RISCVV_INIT(_Tpv, _Tp, suffix) \
inline _Tp##m1_t vreinterpretq_##suffix##_##suffix(_Tp##m1_t v) { return v; } \
inline v_uint8x16 v_reinterpret_as_u8(const v_##_Tpv& v) { return v_uint8x16((vuint8m1_t)(v.val)); } \
inline v_int8x16 v_reinterpret_as_s8(const v_##_Tpv& v) { return v_int8x16((vint8m1_t)(v.val)); } \
inline v_uint16x8 v_reinterpret_as_u16(const v_##_Tpv& v) { return v_uint16x8((vuint16m1_t)(v.val)); } \
inline v_int16x8 v_reinterpret_as_s16(const v_##_Tpv& v) { return v_int16x8((vint16m1_t)(v.val)); } \
inline v_uint32x4 v_reinterpret_as_u32(const v_##_Tpv& v) { return v_uint32x4((vuint32m1_t)(v.val)); } \
inline v_int32x4 v_reinterpret_as_s32(const v_##_Tpv& v) { return v_int32x4((vint32m1_t)(v.val)); } \
inline v_uint64x2 v_reinterpret_as_u64(const v_##_Tpv& v) { return v_uint64x2((vuint64m1_t)(v.val)); } \
inline v_int64x2 v_reinterpret_as_s64(const v_##_Tpv& v) { return v_int64x2((vint64m1_t)(v.val)); } \
inline v_float32x4 v_reinterpret_as_f32(const v_##_Tpv& v) { return v_float32x4((vfloat32m1_t)(v.val)); }\
inline v_float64x2 v_reinterpret_as_f64(const v_##_Tpv& v) { return v_float64x2((vfloat64m1_t)(v.val)); }


OPENCV_HAL_IMPL_RISCVV_INIT(uint8x16, vuint8, u8)
OPENCV_HAL_IMPL_RISCVV_INIT(int8x16, vint8, s8)
OPENCV_HAL_IMPL_RISCVV_INIT(uint16x8, vuint16, u16)
OPENCV_HAL_IMPL_RISCVV_INIT(int16x8, vint16, s16)
OPENCV_HAL_IMPL_RISCVV_INIT(uint32x4, vuint32, u32)
OPENCV_HAL_IMPL_RISCVV_INIT(int32x4, vint32, s32)
OPENCV_HAL_IMPL_RISCVV_INIT(uint64x2, vuint64, u64)
OPENCV_HAL_IMPL_RISCVV_INIT(int64x2, vint64, s64)
OPENCV_HAL_IMPL_RISCVV_INIT(float64x2, vfloat64, f64)
OPENCV_HAL_IMPL_RISCVV_INIT(float32x4, vfloat32, f32)
#define OPENCV_HAL_IMPL_RISCVV_INIT_SET(__Tp, _Tp, suffix, len, num) \
inline v_##_Tp##x##num v_setzero_##suffix() { return v_##_Tp##x##num((v##_Tp##m1_t){0}); }     \
inline v_##_Tp##x##num v_setall_##suffix(__Tp v) { return v_##_Tp##x##num(vmv_v_x_##len##m1(v, num)); }

OPENCV_HAL_IMPL_RISCVV_INIT_SET(uchar, uint8, u8, u8, 16)
OPENCV_HAL_IMPL_RISCVV_INIT_SET(char, int8, s8, i8, 16)
OPENCV_HAL_IMPL_RISCVV_INIT_SET(ushort, uint16, u16, u16, 8)
OPENCV_HAL_IMPL_RISCVV_INIT_SET(short, int16, s16, i16, 8)
OPENCV_HAL_IMPL_RISCVV_INIT_SET(unsigned int, uint32, u32, u32, 4)
OPENCV_HAL_IMPL_RISCVV_INIT_SET(int, int32, s32, i32, 4)
OPENCV_HAL_IMPL_RISCVV_INIT_SET(unsigned long, uint64, u64, u64, 2)
OPENCV_HAL_IMPL_RISCVV_INIT_SET(long, int64, s64, i64, 2)
inline v_float32x4 v_setzero_f32() { return v_float32x4((vfloat32m1_t){0}); }
inline v_float32x4 v_setall_f32(float v) { return v_float32x4(vfmv_v_f_f32m1(v, 4)); }

inline v_float64x2 v_setzero_f64() { return v_float64x2(vfmv_v_f_f64m1(0, 2)); }
inline v_float64x2 v_setall_f64(double v) { return v_float64x2(vfmv_v_f_f64m1(v, 2)); }

inline void v_cleanup() {}

CV_CPU_OPTIMIZATION_HAL_NAMESPACE_END

//! @endcond

}
#endif
