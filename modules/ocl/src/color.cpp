/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2010-2012, Institute Of Software Chinese Academy Of Science, all rights reserved.
// Copyright (C) 2010-2012, Advanced Micro Devices, Inc., all rights reserved.
// Third party copyrights are property of their respective owners.
//
// @Authors
//    Wang Weiyan, wangweiyanster@gmail.com
//    Peng Xiao, pengxiao@multicorewareinc.com
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

#include "precomp.hpp"
#include "opencl_kernels.hpp"

using namespace cv;
using namespace cv::ocl;

static void fromRGB_caller(const oclMat &src, oclMat &dst, int bidx, const std::string & kernelName,
                           const oclMat & data = oclMat())
{
    int src_offset = src.offset / src.elemSize1(), src_step = src.step1();
    int dst_offset = dst.offset / dst.elemSize1(), dst_step = dst.step1();

    std::string build_options = format("-D DEPTH_%d", src.depth());

    vector<pair<size_t , const void *> > args;
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst.cols));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst.rows));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&src_step));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst_step));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&bidx));
    args.push_back( make_pair( sizeof(cl_mem) , (void *)&src.data));
    args.push_back( make_pair( sizeof(cl_mem) , (void *)&dst.data));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&src_offset ));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst_offset ));

    if (!data.empty())
        args.push_back( make_pair( sizeof(cl_mem) , (void *)&data.data ));

    size_t gt[3] = { dst.cols, dst.rows, 1 }, lt[3] = { 16, 16, 1 };
    openCLExecuteKernel(src.clCxt, &cvt_color, kernelName.c_str(), gt, lt, args, -1, -1, build_options.c_str());
}

static void toRGB_caller(const oclMat &src, oclMat &dst, int bidx, const std::string & kernelName,
                         const oclMat & data = oclMat())
{
    std::string build_options = format("-D DEPTH_%d -D dcn=%d", src.depth(), dst.channels());
    int src_offset = src.offset / src.elemSize1(), src_step = src.step1();
    int dst_offset = dst.offset / dst.elemSize1(), dst_step = dst.step1();

    vector<pair<size_t , const void *> > args;
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst.cols));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst.rows));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&src_step));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst_step));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&bidx));
    args.push_back( make_pair( sizeof(cl_mem) , (void *)&src.data));
    args.push_back( make_pair( sizeof(cl_mem) , (void *)&dst.data));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&src_offset ));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst_offset ));

    if (!data.empty())
        args.push_back( make_pair( sizeof(cl_mem) , (void *)&data.data ));

    size_t gt[3] = { dst.cols, dst.rows, 1 }, lt[3] = { 16, 16, 1 };
    openCLExecuteKernel(src.clCxt, &cvt_color, kernelName.c_str(), gt, lt, args, -1, -1, build_options.c_str());
}

static void RGB_caller(const oclMat &src, oclMat &dst, bool reverse)
{
    std::string build_options = format("-D DEPTH_%d -D dcn=%d -D scn=%d -D %s", src.depth(),
                                       dst.channels(), src.channels(), reverse ? "REVERSE" : "ORDER");
    int src_offset = src.offset / src.elemSize1(), src_step = src.step1();
    int dst_offset = dst.offset / dst.elemSize1(), dst_step = dst.step1();

    vector<pair<size_t , const void *> > args;
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst.cols));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst.rows));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&src_step));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst_step));
    args.push_back( make_pair( sizeof(cl_mem) , (void *)&src.data));
    args.push_back( make_pair( sizeof(cl_mem) , (void *)&dst.data));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&src_offset ));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst_offset ));

    size_t gt[3] = { dst.cols, dst.rows, 1 }, lt[3] = { 16, 16, 1 };
    openCLExecuteKernel(src.clCxt, &cvt_color, "RGB", gt, lt, args, -1, -1, build_options.c_str());
}

static void fromRGB5x5_caller(const oclMat &src, oclMat &dst, int bidx, int greenbits, const std::string & kernelName)
{
    std::string build_options = format("-D DEPTH_%d -D greenbits=%d -D dcn=%d",
                                       src.depth(), greenbits, dst.channels());
    int src_offset = src.offset >> 1, src_step = src.step >> 1;
    int dst_offset = dst.offset / dst.elemSize1(), dst_step = dst.step / dst.elemSize1();

    vector<pair<size_t , const void *> > args;
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst.cols));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst.rows));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&src_step));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst_step));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&bidx));
    args.push_back( make_pair( sizeof(cl_mem) , (void *)&src.data));
    args.push_back( make_pair( sizeof(cl_mem) , (void *)&dst.data));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&src_offset ));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst_offset ));

    size_t gt[3] = { dst.cols, dst.rows, 1 }, lt[3] = { 16, 16, 1 };
    openCLExecuteKernel(src.clCxt, &cvt_color, kernelName.c_str(), gt, lt, args, -1, -1, build_options.c_str());
}

static void toRGB5x5_caller(const oclMat &src, oclMat &dst, int bidx, int greenbits, const std::string & kernelName)
{
    std::string build_options = format("-D DEPTH_%d -D greenbits=%d -D scn=%d",
                                       src.depth(), greenbits, src.channels());
    int src_offset = (int)src.offset, src_step = (int)src.step;
    int dst_offset = dst.offset >> 1, dst_step = dst.step >> 1;

    vector<pair<size_t , const void *> > args;
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst.cols));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst.rows));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&src_step));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst_step));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&bidx));
    args.push_back( make_pair( sizeof(cl_mem) , (void *)&src.data));
    args.push_back( make_pair( sizeof(cl_mem) , (void *)&dst.data));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&src_offset ));
    args.push_back( make_pair( sizeof(cl_int) , (void *)&dst_offset ));

    size_t gt[3] = { dst.cols, dst.rows, 1 }, lt[3] = { 16, 16, 1 };
    openCLExecuteKernel(src.clCxt, &cvt_color, kernelName.c_str(), gt, lt, args, -1, -1, build_options.c_str());
}

static void cvtColor_caller(const oclMat &src, oclMat &dst, int code, int dcn)
{
    Size sz = src.size();
    int scn = src.channels(), depth = src.depth(), bidx;

    CV_Assert(depth == CV_8U || depth == CV_16U || depth == CV_32F);

    switch (code)
    {
    case CV_BGR2BGRA: case CV_RGB2BGRA: case CV_BGRA2BGR:
    case CV_RGBA2BGR: case CV_RGB2BGR: case CV_BGRA2RGBA:
    {
        CV_Assert(scn == 3 || scn == 4);
        dcn = code == CV_BGR2BGRA || code == CV_RGB2BGRA || code == CV_BGRA2RGBA ? 4 : 3;
        bool reverse = !(code == CV_BGR2BGRA || code == CV_BGRA2BGR);
        dst.create(sz, CV_MAKE_TYPE(depth, dcn));
        RGB_caller(src, dst, reverse);
        break;
    }
    case CV_BGR2BGR565: case CV_BGR2BGR555: case CV_RGB2BGR565: case CV_RGB2BGR555:
    case CV_BGRA2BGR565: case CV_BGRA2BGR555: case CV_RGBA2BGR565: case CV_RGBA2BGR555:
    {
        CV_Assert((scn == 3 || scn == 4) && depth == CV_8U );
        bidx = code == CV_BGR2BGR565 || code == CV_BGR2BGR555 ||
            code == CV_BGRA2BGR565 || code == CV_BGRA2BGR555 ? 0 : 2;
        int greenbits = code == CV_BGR2BGR565 || code == CV_RGB2BGR565 ||
            code == CV_BGRA2BGR565 || code == CV_RGBA2BGR565 ? 6 : 5;
        dst.create(sz, CV_8UC2);
        toRGB5x5_caller(src, dst, bidx, greenbits, "RGB2RGB5x5");
        break;
    }
    case CV_BGR5652BGR: case CV_BGR5552BGR: case CV_BGR5652RGB: case CV_BGR5552RGB:
    case CV_BGR5652BGRA: case CV_BGR5552BGRA: case CV_BGR5652RGBA: case CV_BGR5552RGBA:
    {
        dcn = code == CV_BGR5652BGRA || code == CV_BGR5552BGRA || code == CV_BGR5652RGBA || code == CV_BGR5552RGBA ? 4 : 3;
        CV_Assert((dcn == 3 || dcn == 4) && scn == 2 && depth == CV_8U);
        bidx = code == CV_BGR5652BGR || code == CV_BGR5552BGR ||
            code == CV_BGR5652BGRA || code == CV_BGR5552BGRA ? 0 : 2;
        int greenbits = code == CV_BGR5652BGR || code == CV_BGR5652RGB ||
            code == CV_BGR5652BGRA || code == CV_BGR5652RGBA ? 6 : 5;
        dst.create(sz, CV_MAKETYPE(depth, dcn));
        fromRGB5x5_caller(src, dst, bidx, greenbits, "RGB5x52RGB");
        break;
    }
    case CV_BGR5652GRAY: case CV_BGR5552GRAY:
    {
        CV_Assert(scn == 2 && depth == CV_8U);
        dst.create(sz, CV_8UC1);
        int greenbits = code == CV_BGR5652GRAY ? 6 : 5;
        fromRGB5x5_caller(src, dst, -1, greenbits, "BGR5x52Gray");
        break;
    }
    case CV_GRAY2BGR565: case CV_GRAY2BGR555:
    {
        CV_Assert(scn == 1 && depth == CV_8U);
        dst.create(sz, CV_8UC2);
        int greenbits = code == CV_GRAY2BGR565 ? 6 : 5;
        toRGB5x5_caller(src, dst, -1, greenbits, "Gray2BGR5x5");
        break;
    }
    case CV_RGB2GRAY: case CV_BGR2GRAY:
    case CV_RGBA2GRAY: case CV_BGRA2GRAY:
    {
        CV_Assert(scn == 3 || scn == 4);
        bidx = code == CV_BGR2GRAY || code == CV_BGRA2GRAY ? 0 : 2;
        dst.create(sz, CV_MAKETYPE(depth, 1));
        fromRGB_caller(src, dst, bidx, "RGB2Gray");
        break;
    }
    case CV_GRAY2BGR:
    case CV_GRAY2BGRA:
    {
        CV_Assert(scn == 1);
        dcn  = code == CV_GRAY2BGRA ? 4 : 3;
        dst.create(sz, CV_MAKETYPE(depth, dcn));
        toRGB_caller(src, dst, 0, "Gray2RGB");
        break;
    }
    case CV_BGR2YUV:
    case CV_RGB2YUV:
    {
        CV_Assert(scn == 3 || scn == 4);
        bidx = code == CV_BGR2YUV ? 0 : 2;
        dst.create(sz, CV_MAKETYPE(depth, 3));
        fromRGB_caller(src, dst, bidx, "RGB2YUV");
        break;
    }
    case CV_YUV2BGR:
    case CV_YUV2RGB:
    {
        if( dcn <= 0 )
            dcn = 3;
        CV_Assert(scn == 3 && (dcn == 3 || dcn == 4));
        bidx = code == CV_YUV2BGR ? 0 : 2;
        dst.create(sz, CV_MAKETYPE(depth, dcn));
        toRGB_caller(src, dst, bidx, "YUV2RGB");
        break;
    }
    case CV_YUV2RGB_NV12: case CV_YUV2BGR_NV12:
    case CV_YUV2RGBA_NV12: case CV_YUV2BGRA_NV12:
    {
        CV_Assert(scn == 1);
        CV_Assert( sz.width % 2 == 0 && sz.height % 3 == 0 && depth == CV_8U );
        dcn = code == CV_YUV2BGRA_NV12 || code == CV_YUV2RGBA_NV12 ? 4 : 3;
        bidx = code == CV_YUV2BGRA_NV12 || code == CV_YUV2BGR_NV12 ? 0 : 2;

        Size dstSz(sz.width, sz.height * 2 / 3);
        dst.create(dstSz, CV_MAKETYPE(depth, dcn));
        toRGB_caller(src, dst, bidx, "YUV2RGBA_NV12");
        break;
    }
    case CV_BGR2YCrCb:
    case CV_RGB2YCrCb:
    {
        CV_Assert(scn == 3 || scn == 4);
        bidx = code == CV_BGR2YCrCb ? 0 : 2;
        dst.create(sz, CV_MAKETYPE(depth, 3));
        fromRGB_caller(src, dst, bidx, "RGB2YCrCb");
        break;
    }
    case CV_YCrCb2BGR:
    case CV_YCrCb2RGB:
    {
        if( dcn <= 0 )
            dcn = 3;
        CV_Assert(scn == 3 && (dcn == 3 || dcn == 4));
        bidx = code == CV_YCrCb2BGR ? 0 : 2;
        dst.create(sz, CV_MAKETYPE(depth, dcn));
        toRGB_caller(src, dst, bidx, "YCrCb2RGB");
        break;
    }
    /*
    case CV_BGR5652GRAY: case CV_BGR5552GRAY:
    case CV_GRAY2BGR565: case CV_GRAY2BGR555:
    */
    case CV_BGR2XYZ:
    case CV_RGB2XYZ:
    {
        CV_Assert(scn == 3 || scn == 4);
        bidx = code == CV_BGR2XYZ ? 0 : 2;
        dst.create(sz, CV_MAKE_TYPE(depth, 3));

        void * pdata = NULL;
        if (depth == CV_32F)
        {
            float coeffs[] =
            {
                0.412453f, 0.357580f, 0.180423f,
                0.212671f, 0.715160f, 0.072169f,
                0.019334f, 0.119193f, 0.950227f
            };
            if (bidx == 0)
            {
                std::swap(coeffs[0], coeffs[2]);
                std::swap(coeffs[3], coeffs[5]);
                std::swap(coeffs[6], coeffs[8]);
            }
            pdata = coeffs;
        }
        else
        {
            int coeffs[] =
            {
                1689,    1465,    739,
                871,     2929,    296,
                79,      488,     3892
            };
            if (bidx == 0)
            {
                std::swap(coeffs[0], coeffs[2]);
                std::swap(coeffs[3], coeffs[5]);
                std::swap(coeffs[6], coeffs[8]);
            }
            pdata = coeffs;
        }
        oclMat oclCoeffs(1, 9, depth == CV_32F ? CV_32FC1 : CV_32SC1, pdata);

        fromRGB_caller(src, dst, bidx, "RGB2XYZ", oclCoeffs);
        break;
    }
    case CV_XYZ2BGR:
    case CV_XYZ2RGB:
    {
        if (dcn <= 0)
            dcn = 3;
        CV_Assert(scn == 3 && (dcn == 3 || dcn == 4));
        bidx = code == CV_XYZ2BGR ? 0 : 2;
        dst.create(sz, CV_MAKE_TYPE(depth, dcn));

        void * pdata = NULL;
        if (depth == CV_32F)
        {
            float coeffs[] =
            {
                3.240479f, -1.53715f, -0.498535f,
                -0.969256f, 1.875991f, 0.041556f,
                0.055648f, -0.204043f, 1.057311f
            };
            if (bidx == 0)
            {
                std::swap(coeffs[0], coeffs[6]);
                std::swap(coeffs[1], coeffs[7]);
                std::swap(coeffs[2], coeffs[8]);
            }
            pdata = coeffs;
        }
        else
        {
            int coeffs[] =
            {
                13273,  -6296,  -2042,
                -3970,   7684,    170,
                  228,   -836,   4331
            };
            if (bidx == 0)
            {
                std::swap(coeffs[0], coeffs[6]);
                std::swap(coeffs[1], coeffs[7]);
                std::swap(coeffs[2], coeffs[8]);
            }
            pdata = coeffs;
        }
        oclMat oclCoeffs(1, 9, depth == CV_32F ? CV_32FC1 : CV_32SC1, pdata);

        toRGB_caller(src, dst, bidx, "XYZ2RGB", oclCoeffs);
        break;
    }
    /*
    case CV_BGR2HSV: case CV_RGB2HSV: case CV_BGR2HSV_FULL: case CV_RGB2HSV_FULL:
    case CV_BGR2HLS: case CV_RGB2HLS: case CV_BGR2HLS_FULL: case CV_RGB2HLS_FULL:
    case CV_HSV2BGR: case CV_HSV2RGB: case CV_HSV2BGR_FULL: case CV_HSV2RGB_FULL:
    case CV_HLS2BGR: case CV_HLS2RGB: case CV_HLS2BGR_FULL: case CV_HLS2RGB_FULL:
    */
    default:
        CV_Error( CV_StsBadFlag, "Unknown/unsupported color conversion code" );
    }
}

void cv::ocl::cvtColor(const oclMat &src, oclMat &dst, int code, int dcn)
{
    cvtColor_caller(src, dst, code, dcn);
}
