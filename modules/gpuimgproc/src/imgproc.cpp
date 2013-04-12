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
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
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

#include "precomp.hpp"

using namespace cv;
using namespace cv::gpu;

#if !defined (HAVE_CUDA) || defined (CUDA_DISABLER)

void cv::gpu::meanShiftFiltering(const GpuMat&, GpuMat&, int, int, TermCriteria, Stream&) { throw_no_cuda(); }
void cv::gpu::meanShiftProc(const GpuMat&, GpuMat&, GpuMat&, int, int, TermCriteria, Stream&) { throw_no_cuda(); }
void cv::gpu::drawColorDisp(const GpuMat&, GpuMat&, int, Stream&) { throw_no_cuda(); }
void cv::gpu::reprojectImageTo3D(const GpuMat&, GpuMat&, const Mat&, int, Stream&) { throw_no_cuda(); }
void cv::gpu::buildWarpPlaneMaps(Size, Rect, const Mat&, const Mat&, const Mat&, float, GpuMat&, GpuMat&, Stream&) { throw_no_cuda(); }
void cv::gpu::buildWarpCylindricalMaps(Size, Rect, const Mat&, const Mat&, float, GpuMat&, GpuMat&, Stream&) { throw_no_cuda(); }
void cv::gpu::buildWarpSphericalMaps(Size, Rect, const Mat&, const Mat&, float, GpuMat&, GpuMat&, Stream&) { throw_no_cuda(); }
void cv::gpu::rotate(const GpuMat&, GpuMat&, Size, double, double, double, int, Stream&) { throw_no_cuda(); }
void cv::gpu::evenLevels(GpuMat&, int, int, int) { throw_no_cuda(); }
void cv::gpu::histEven(const GpuMat&, GpuMat&, int, int, int, Stream&) { throw_no_cuda(); }
void cv::gpu::histEven(const GpuMat&, GpuMat&, GpuMat&, int, int, int, Stream&) { throw_no_cuda(); }
void cv::gpu::histEven(const GpuMat&, GpuMat*, int*, int*, int*, Stream&) { throw_no_cuda(); }
void cv::gpu::histEven(const GpuMat&, GpuMat*, GpuMat&, int*, int*, int*, Stream&) { throw_no_cuda(); }
void cv::gpu::histRange(const GpuMat&, GpuMat&, const GpuMat&, Stream&) { throw_no_cuda(); }
void cv::gpu::histRange(const GpuMat&, GpuMat&, const GpuMat&, GpuMat&, Stream&) { throw_no_cuda(); }
void cv::gpu::histRange(const GpuMat&, GpuMat*, const GpuMat*, Stream&) { throw_no_cuda(); }
void cv::gpu::histRange(const GpuMat&, GpuMat*, const GpuMat*, GpuMat&, Stream&) { throw_no_cuda(); }
void cv::gpu::calcHist(const GpuMat&, GpuMat&, Stream&) { throw_no_cuda(); }
void cv::gpu::equalizeHist(const GpuMat&, GpuMat&, Stream&) { throw_no_cuda(); }
void cv::gpu::equalizeHist(const GpuMat&, GpuMat&, GpuMat&, GpuMat&, Stream&) { throw_no_cuda(); }
void cv::gpu::cornerHarris(const GpuMat&, GpuMat&, int, int, double, int) { throw_no_cuda(); }
void cv::gpu::cornerHarris(const GpuMat&, GpuMat&, GpuMat&, GpuMat&, int, int, double, int) { throw_no_cuda(); }
void cv::gpu::cornerHarris(const GpuMat&, GpuMat&, GpuMat&, GpuMat&, GpuMat&, int, int, double, int, Stream&) { throw_no_cuda(); }
void cv::gpu::cornerMinEigenVal(const GpuMat&, GpuMat&, int, int, int) { throw_no_cuda(); }
void cv::gpu::cornerMinEigenVal(const GpuMat&, GpuMat&, GpuMat&, GpuMat&, int, int, int) { throw_no_cuda(); }
void cv::gpu::cornerMinEigenVal(const GpuMat&, GpuMat&, GpuMat&, GpuMat&, GpuMat&, int, int, int, Stream&) { throw_no_cuda(); }
void cv::gpu::Canny(const GpuMat&, GpuMat&, double, double, int, bool) { throw_no_cuda(); }
void cv::gpu::Canny(const GpuMat&, CannyBuf&, GpuMat&, double, double, int, bool) { throw_no_cuda(); }
void cv::gpu::Canny(const GpuMat&, const GpuMat&, GpuMat&, double, double, bool) { throw_no_cuda(); }
void cv::gpu::Canny(const GpuMat&, const GpuMat&, CannyBuf&, GpuMat&, double, double, bool) { throw_no_cuda(); }
void cv::gpu::CannyBuf::create(const Size&, int) { throw_no_cuda(); }
void cv::gpu::CannyBuf::release() { throw_no_cuda(); }
cv::Ptr<cv::gpu::CLAHE> cv::gpu::createCLAHE(double, cv::Size) { throw_no_cuda(); return cv::Ptr<cv::gpu::CLAHE>(); }
void cv::gpu::alphaComp(const GpuMat&, const GpuMat&, GpuMat&, int, Stream&) { throw_no_cuda(); }

#else /* !defined (HAVE_CUDA) */

////////////////////////////////////////////////////////////////////////
// meanShiftFiltering_GPU

namespace cv { namespace gpu { namespace cudev
{
    namespace imgproc
    {
        void meanShiftFiltering_gpu(const PtrStepSzb& src, PtrStepSzb dst, int sp, int sr, int maxIter, float eps, cudaStream_t stream);
    }
}}}

void cv::gpu::meanShiftFiltering(const GpuMat& src, GpuMat& dst, int sp, int sr, TermCriteria criteria, Stream& stream)
{
    using namespace ::cv::gpu::cudev::imgproc;

    if( src.empty() )
        CV_Error( cv::Error::StsBadArg, "The input image is empty" );

    if( src.depth() != CV_8U || src.channels() != 4 )
        CV_Error( cv::Error::StsUnsupportedFormat, "Only 8-bit, 4-channel images are supported" );

    dst.create( src.size(), CV_8UC4 );

    if( !(criteria.type & TermCriteria::MAX_ITER) )
        criteria.maxCount = 5;

    int maxIter = std::min(std::max(criteria.maxCount, 1), 100);

    float eps;
    if( !(criteria.type & TermCriteria::EPS) )
        eps = 1.f;
    eps = (float)std::max(criteria.epsilon, 0.0);

    meanShiftFiltering_gpu(src, dst, sp, sr, maxIter, eps, StreamAccessor::getStream(stream));
}

////////////////////////////////////////////////////////////////////////
// meanShiftProc_GPU

namespace cv { namespace gpu { namespace cudev
{
    namespace imgproc
    {
        void meanShiftProc_gpu(const PtrStepSzb& src, PtrStepSzb dstr, PtrStepSzb dstsp, int sp, int sr, int maxIter, float eps, cudaStream_t stream);
    }
}}}

void cv::gpu::meanShiftProc(const GpuMat& src, GpuMat& dstr, GpuMat& dstsp, int sp, int sr, TermCriteria criteria, Stream& stream)
{
    using namespace ::cv::gpu::cudev::imgproc;

    if( src.empty() )
        CV_Error( cv::Error::StsBadArg, "The input image is empty" );

    if( src.depth() != CV_8U || src.channels() != 4 )
        CV_Error( cv::Error::StsUnsupportedFormat, "Only 8-bit, 4-channel images are supported" );

    dstr.create( src.size(), CV_8UC4 );
    dstsp.create( src.size(), CV_16SC2 );

    if( !(criteria.type & TermCriteria::MAX_ITER) )
        criteria.maxCount = 5;

    int maxIter = std::min(std::max(criteria.maxCount, 1), 100);

    float eps;
    if( !(criteria.type & TermCriteria::EPS) )
        eps = 1.f;
    eps = (float)std::max(criteria.epsilon, 0.0);

    meanShiftProc_gpu(src, dstr, dstsp, sp, sr, maxIter, eps, StreamAccessor::getStream(stream));
}

////////////////////////////////////////////////////////////////////////
// drawColorDisp

namespace cv { namespace gpu { namespace cudev
{
    namespace imgproc
    {
        void drawColorDisp_gpu(const PtrStepSzb& src, const PtrStepSzb& dst, int ndisp, const cudaStream_t& stream);
        void drawColorDisp_gpu(const PtrStepSz<short>& src, const PtrStepSzb& dst, int ndisp, const cudaStream_t& stream);
    }
}}}

namespace
{
    template <typename T>
    void drawColorDisp_caller(const GpuMat& src, GpuMat& dst, int ndisp, const cudaStream_t& stream)
    {
        using namespace ::cv::gpu::cudev::imgproc;

        dst.create(src.size(), CV_8UC4);

        drawColorDisp_gpu((PtrStepSz<T>)src, dst, ndisp, stream);
    }

    typedef void (*drawColorDisp_caller_t)(const GpuMat& src, GpuMat& dst, int ndisp, const cudaStream_t& stream);

    const drawColorDisp_caller_t drawColorDisp_callers[] = {drawColorDisp_caller<unsigned char>, 0, 0, drawColorDisp_caller<short>, 0, 0, 0, 0};
}

void cv::gpu::drawColorDisp(const GpuMat& src, GpuMat& dst, int ndisp, Stream& stream)
{
    CV_Assert(src.type() == CV_8U || src.type() == CV_16S);

    drawColorDisp_callers[src.type()](src, dst, ndisp, StreamAccessor::getStream(stream));
}

////////////////////////////////////////////////////////////////////////
// reprojectImageTo3D

namespace cv { namespace gpu { namespace cudev
{
    namespace imgproc
    {
        template <typename T, typename D>
        void reprojectImageTo3D_gpu(const PtrStepSzb disp, PtrStepSzb xyz, const float* q, cudaStream_t stream);
    }
}}}

void cv::gpu::reprojectImageTo3D(const GpuMat& disp, GpuMat& xyz, const Mat& Q, int dst_cn, Stream& stream)
{
    using namespace cv::gpu::cudev::imgproc;

    typedef void (*func_t)(const PtrStepSzb disp, PtrStepSzb xyz, const float* q, cudaStream_t stream);
    static const func_t funcs[2][4] =
    {
        {reprojectImageTo3D_gpu<uchar, float3>, 0, 0, reprojectImageTo3D_gpu<short, float3>},
        {reprojectImageTo3D_gpu<uchar, float4>, 0, 0, reprojectImageTo3D_gpu<short, float4>}
    };

    CV_Assert(disp.type() == CV_8U || disp.type() == CV_16S);
    CV_Assert(Q.type() == CV_32F && Q.rows == 4 && Q.cols == 4 && Q.isContinuous());
    CV_Assert(dst_cn == 3 || dst_cn == 4);

    xyz.create(disp.size(), CV_MAKE_TYPE(CV_32F, dst_cn));

    funcs[dst_cn == 4][disp.type()](disp, xyz, Q.ptr<float>(), StreamAccessor::getStream(stream));
}

//////////////////////////////////////////////////////////////////////////////
// buildWarpPlaneMaps

namespace cv { namespace gpu { namespace cudev
{
    namespace imgproc
    {
        void buildWarpPlaneMaps(int tl_u, int tl_v, PtrStepSzf map_x, PtrStepSzf map_y,
                                const float k_rinv[9], const float r_kinv[9], const float t[3], float scale,
                                cudaStream_t stream);
    }
}}}

void cv::gpu::buildWarpPlaneMaps(Size src_size, Rect dst_roi, const Mat &K, const Mat& R, const Mat &T,
                                 float scale, GpuMat& map_x, GpuMat& map_y, Stream& stream)
{
    (void)src_size;
    using namespace ::cv::gpu::cudev::imgproc;

    CV_Assert(K.size() == Size(3,3) && K.type() == CV_32F);
    CV_Assert(R.size() == Size(3,3) && R.type() == CV_32F);
    CV_Assert((T.size() == Size(3,1) || T.size() == Size(1,3)) && T.type() == CV_32F && T.isContinuous());

    Mat K_Rinv = K * R.t();
    Mat R_Kinv = R * K.inv();
    CV_Assert(K_Rinv.isContinuous());
    CV_Assert(R_Kinv.isContinuous());

    map_x.create(dst_roi.size(), CV_32F);
    map_y.create(dst_roi.size(), CV_32F);
    cudev::imgproc::buildWarpPlaneMaps(dst_roi.tl().x, dst_roi.tl().y, map_x, map_y, K_Rinv.ptr<float>(), R_Kinv.ptr<float>(),
                       T.ptr<float>(), scale, StreamAccessor::getStream(stream));
}

//////////////////////////////////////////////////////////////////////////////
// buildWarpCylyndricalMaps

namespace cv { namespace gpu { namespace cudev
{
    namespace imgproc
    {
        void buildWarpCylindricalMaps(int tl_u, int tl_v, PtrStepSzf map_x, PtrStepSzf map_y,
                                      const float k_rinv[9], const float r_kinv[9], float scale,
                                      cudaStream_t stream);
    }
}}}

void cv::gpu::buildWarpCylindricalMaps(Size src_size, Rect dst_roi, const Mat &K, const Mat& R, float scale,
                                       GpuMat& map_x, GpuMat& map_y, Stream& stream)
{
    (void)src_size;
    using namespace ::cv::gpu::cudev::imgproc;

    CV_Assert(K.size() == Size(3,3) && K.type() == CV_32F);
    CV_Assert(R.size() == Size(3,3) && R.type() == CV_32F);

    Mat K_Rinv = K * R.t();
    Mat R_Kinv = R * K.inv();
    CV_Assert(K_Rinv.isContinuous());
    CV_Assert(R_Kinv.isContinuous());

    map_x.create(dst_roi.size(), CV_32F);
    map_y.create(dst_roi.size(), CV_32F);
    cudev::imgproc::buildWarpCylindricalMaps(dst_roi.tl().x, dst_roi.tl().y, map_x, map_y, K_Rinv.ptr<float>(), R_Kinv.ptr<float>(), scale, StreamAccessor::getStream(stream));
}


//////////////////////////////////////////////////////////////////////////////
// buildWarpSphericalMaps

namespace cv { namespace gpu { namespace cudev
{
    namespace imgproc
    {
        void buildWarpSphericalMaps(int tl_u, int tl_v, PtrStepSzf map_x, PtrStepSzf map_y,
                                    const float k_rinv[9], const float r_kinv[9], float scale,
                                    cudaStream_t stream);
    }
}}}

void cv::gpu::buildWarpSphericalMaps(Size src_size, Rect dst_roi, const Mat &K, const Mat& R, float scale,
                                     GpuMat& map_x, GpuMat& map_y, Stream& stream)
{
    (void)src_size;
    using namespace ::cv::gpu::cudev::imgproc;

    CV_Assert(K.size() == Size(3,3) && K.type() == CV_32F);
    CV_Assert(R.size() == Size(3,3) && R.type() == CV_32F);

    Mat K_Rinv = K * R.t();
    Mat R_Kinv = R * K.inv();
    CV_Assert(K_Rinv.isContinuous());
    CV_Assert(R_Kinv.isContinuous());

    map_x.create(dst_roi.size(), CV_32F);
    map_y.create(dst_roi.size(), CV_32F);
    cudev::imgproc::buildWarpSphericalMaps(dst_roi.tl().x, dst_roi.tl().y, map_x, map_y, K_Rinv.ptr<float>(), R_Kinv.ptr<float>(), scale, StreamAccessor::getStream(stream));
}

////////////////////////////////////////////////////////////////////////
// rotate

namespace
{
    template<int DEPTH> struct NppTypeTraits;
    template<> struct NppTypeTraits<CV_8U>  { typedef Npp8u npp_t; };
    template<> struct NppTypeTraits<CV_8S>  { typedef Npp8s npp_t; };
    template<> struct NppTypeTraits<CV_16U> { typedef Npp16u npp_t; };
    template<> struct NppTypeTraits<CV_16S> { typedef Npp16s npp_t; };
    template<> struct NppTypeTraits<CV_32S> { typedef Npp32s npp_t; };
    template<> struct NppTypeTraits<CV_32F> { typedef Npp32f npp_t; };
    template<> struct NppTypeTraits<CV_64F> { typedef Npp64f npp_t; };

    template <int DEPTH> struct NppRotateFunc
    {
        typedef typename NppTypeTraits<DEPTH>::npp_t npp_t;

        typedef NppStatus (*func_t)(const npp_t* pSrc, NppiSize oSrcSize, int nSrcStep, NppiRect oSrcROI,
                                    npp_t* pDst, int nDstStep, NppiRect oDstROI,
                                    double nAngle, double nShiftX, double nShiftY, int eInterpolation);
    };

    template <int DEPTH, typename NppRotateFunc<DEPTH>::func_t func> struct NppRotate
    {
        typedef typename NppRotateFunc<DEPTH>::npp_t npp_t;

        static void call(const GpuMat& src, GpuMat& dst, Size dsize, double angle, double xShift, double yShift, int interpolation, cudaStream_t stream)
        {
            (void)dsize;
            static const int npp_inter[] = {NPPI_INTER_NN, NPPI_INTER_LINEAR, NPPI_INTER_CUBIC};

            NppStreamHandler h(stream);

            NppiSize srcsz;
            srcsz.height = src.rows;
            srcsz.width = src.cols;
            NppiRect srcroi;
            srcroi.x = srcroi.y = 0;
            srcroi.height = src.rows;
            srcroi.width = src.cols;
            NppiRect dstroi;
            dstroi.x = dstroi.y = 0;
            dstroi.height = dst.rows;
            dstroi.width = dst.cols;

            nppSafeCall( func(src.ptr<npp_t>(), srcsz, static_cast<int>(src.step), srcroi,
                dst.ptr<npp_t>(), static_cast<int>(dst.step), dstroi, angle, xShift, yShift, npp_inter[interpolation]) );

            if (stream == 0)
                cudaSafeCall( cudaDeviceSynchronize() );
        }
    };
}

void cv::gpu::rotate(const GpuMat& src, GpuMat& dst, Size dsize, double angle, double xShift, double yShift, int interpolation, Stream& stream)
{
    typedef void (*func_t)(const GpuMat& src, GpuMat& dst, Size dsize, double angle, double xShift, double yShift, int interpolation, cudaStream_t stream);

    static const func_t funcs[6][4] =
    {
        {NppRotate<CV_8U, nppiRotate_8u_C1R>::call, 0, NppRotate<CV_8U, nppiRotate_8u_C3R>::call, NppRotate<CV_8U, nppiRotate_8u_C4R>::call},
        {0,0,0,0},
        {NppRotate<CV_16U, nppiRotate_16u_C1R>::call, 0, NppRotate<CV_16U, nppiRotate_16u_C3R>::call, NppRotate<CV_16U, nppiRotate_16u_C4R>::call},
        {0,0,0,0},
        {0,0,0,0},
        {NppRotate<CV_32F, nppiRotate_32f_C1R>::call, 0, NppRotate<CV_32F, nppiRotate_32f_C3R>::call, NppRotate<CV_32F, nppiRotate_32f_C4R>::call}
    };

    CV_Assert(src.depth() == CV_8U || src.depth() == CV_16U || src.depth() == CV_32F);
    CV_Assert(src.channels() == 1 || src.channels() == 3 || src.channels() == 4);
    CV_Assert(interpolation == INTER_NEAREST || interpolation == INTER_LINEAR || interpolation == INTER_CUBIC);

    dst.create(dsize, src.type());
    dst.setTo(Scalar::all(0));

    funcs[src.depth()][src.channels() - 1](src, dst, dsize, angle, xShift, yShift, interpolation, StreamAccessor::getStream(stream));
}


////////////////////////////////////////////////////////////////////////
// Histogram

namespace
{
    typedef NppStatus (*get_buf_size_c1_t)(NppiSize oSizeROI, int nLevels, int* hpBufferSize);
    typedef NppStatus (*get_buf_size_c4_t)(NppiSize oSizeROI, int nLevels[], int* hpBufferSize);

    template<int SDEPTH> struct NppHistogramEvenFuncC1
    {
        typedef typename NppTypeTraits<SDEPTH>::npp_t src_t;

    typedef NppStatus (*func_ptr)(const src_t* pSrc, int nSrcStep, NppiSize oSizeROI, Npp32s * pHist,
            int nLevels, Npp32s nLowerLevel, Npp32s nUpperLevel, Npp8u * pBuffer);
    };
    template<int SDEPTH> struct NppHistogramEvenFuncC4
    {
        typedef typename NppTypeTraits<SDEPTH>::npp_t src_t;

        typedef NppStatus (*func_ptr)(const src_t* pSrc, int nSrcStep, NppiSize oSizeROI,
            Npp32s * pHist[4], int nLevels[4], Npp32s nLowerLevel[4], Npp32s nUpperLevel[4], Npp8u * pBuffer);
    };

    template<int SDEPTH, typename NppHistogramEvenFuncC1<SDEPTH>::func_ptr func, get_buf_size_c1_t get_buf_size>
    struct NppHistogramEvenC1
    {
        typedef typename NppHistogramEvenFuncC1<SDEPTH>::src_t src_t;

        static void hist(const GpuMat& src, GpuMat& hist, GpuMat& buffer, int histSize, int lowerLevel, int upperLevel, cudaStream_t stream)
        {
            int levels = histSize + 1;
            hist.create(1, histSize, CV_32S);

            NppiSize sz;
            sz.width = src.cols;
            sz.height = src.rows;

            int buf_size;
            get_buf_size(sz, levels, &buf_size);

            ensureSizeIsEnough(1, buf_size, CV_8U, buffer);

            NppStreamHandler h(stream);

            nppSafeCall( func(src.ptr<src_t>(), static_cast<int>(src.step), sz, hist.ptr<Npp32s>(), levels,
                lowerLevel, upperLevel, buffer.ptr<Npp8u>()) );

            if (stream == 0)
                cudaSafeCall( cudaDeviceSynchronize() );
        }
    };
    template<int SDEPTH, typename NppHistogramEvenFuncC4<SDEPTH>::func_ptr func, get_buf_size_c4_t get_buf_size>
    struct NppHistogramEvenC4
    {
        typedef typename NppHistogramEvenFuncC4<SDEPTH>::src_t src_t;

        static void hist(const GpuMat& src, GpuMat hist[4], GpuMat& buffer, int histSize[4], int lowerLevel[4], int upperLevel[4], cudaStream_t stream)
        {
            int levels[] = {histSize[0] + 1, histSize[1] + 1, histSize[2] + 1, histSize[3] + 1};
            hist[0].create(1, histSize[0], CV_32S);
            hist[1].create(1, histSize[1], CV_32S);
            hist[2].create(1, histSize[2], CV_32S);
            hist[3].create(1, histSize[3], CV_32S);

            NppiSize sz;
            sz.width = src.cols;
            sz.height = src.rows;

            Npp32s* pHist[] = {hist[0].ptr<Npp32s>(), hist[1].ptr<Npp32s>(), hist[2].ptr<Npp32s>(), hist[3].ptr<Npp32s>()};

            int buf_size;
            get_buf_size(sz, levels, &buf_size);

            ensureSizeIsEnough(1, buf_size, CV_8U, buffer);

            NppStreamHandler h(stream);

            nppSafeCall( func(src.ptr<src_t>(), static_cast<int>(src.step), sz, pHist, levels, lowerLevel, upperLevel, buffer.ptr<Npp8u>()) );

            if (stream == 0)
                cudaSafeCall( cudaDeviceSynchronize() );
        }
    };

    template<int SDEPTH> struct NppHistogramRangeFuncC1
    {
        typedef typename NppTypeTraits<SDEPTH>::npp_t src_t;
        typedef Npp32s level_t;
        enum {LEVEL_TYPE_CODE=CV_32SC1};

        typedef NppStatus (*func_ptr)(const src_t* pSrc, int nSrcStep, NppiSize oSizeROI, Npp32s* pHist,
            const Npp32s* pLevels, int nLevels, Npp8u* pBuffer);
    };
    template<> struct NppHistogramRangeFuncC1<CV_32F>
    {
        typedef Npp32f src_t;
        typedef Npp32f level_t;
        enum {LEVEL_TYPE_CODE=CV_32FC1};

        typedef NppStatus (*func_ptr)(const Npp32f* pSrc, int nSrcStep, NppiSize oSizeROI, Npp32s* pHist,
            const Npp32f* pLevels, int nLevels, Npp8u* pBuffer);
    };
    template<int SDEPTH> struct NppHistogramRangeFuncC4
    {
        typedef typename NppTypeTraits<SDEPTH>::npp_t src_t;
        typedef Npp32s level_t;
        enum {LEVEL_TYPE_CODE=CV_32SC1};

        typedef NppStatus (*func_ptr)(const src_t* pSrc, int nSrcStep, NppiSize oSizeROI, Npp32s* pHist[4],
            const Npp32s* pLevels[4], int nLevels[4], Npp8u* pBuffer);
    };
    template<> struct NppHistogramRangeFuncC4<CV_32F>
    {
        typedef Npp32f src_t;
        typedef Npp32f level_t;
        enum {LEVEL_TYPE_CODE=CV_32FC1};

        typedef NppStatus (*func_ptr)(const Npp32f* pSrc, int nSrcStep, NppiSize oSizeROI, Npp32s* pHist[4],
            const Npp32f* pLevels[4], int nLevels[4], Npp8u* pBuffer);
    };

    template<int SDEPTH, typename NppHistogramRangeFuncC1<SDEPTH>::func_ptr func, get_buf_size_c1_t get_buf_size>
    struct NppHistogramRangeC1
    {
        typedef typename NppHistogramRangeFuncC1<SDEPTH>::src_t src_t;
        typedef typename NppHistogramRangeFuncC1<SDEPTH>::level_t level_t;
        enum {LEVEL_TYPE_CODE=NppHistogramRangeFuncC1<SDEPTH>::LEVEL_TYPE_CODE};

        static void hist(const GpuMat& src, GpuMat& hist, const GpuMat& levels, GpuMat& buffer, cudaStream_t stream)
        {
            CV_Assert(levels.type() == LEVEL_TYPE_CODE && levels.rows == 1);

            hist.create(1, levels.cols - 1, CV_32S);

            NppiSize sz;
            sz.width = src.cols;
            sz.height = src.rows;

            int buf_size;
            get_buf_size(sz, levels.cols, &buf_size);

            ensureSizeIsEnough(1, buf_size, CV_8U, buffer);

            NppStreamHandler h(stream);

            nppSafeCall( func(src.ptr<src_t>(), static_cast<int>(src.step), sz, hist.ptr<Npp32s>(), levels.ptr<level_t>(), levels.cols, buffer.ptr<Npp8u>()) );

            if (stream == 0)
                cudaSafeCall( cudaDeviceSynchronize() );
        }
    };
    template<int SDEPTH, typename NppHistogramRangeFuncC4<SDEPTH>::func_ptr func, get_buf_size_c4_t get_buf_size>
    struct NppHistogramRangeC4
    {
        typedef typename NppHistogramRangeFuncC4<SDEPTH>::src_t src_t;
        typedef typename NppHistogramRangeFuncC1<SDEPTH>::level_t level_t;
        enum {LEVEL_TYPE_CODE=NppHistogramRangeFuncC1<SDEPTH>::LEVEL_TYPE_CODE};

        static void hist(const GpuMat& src, GpuMat hist[4], const GpuMat levels[4], GpuMat& buffer, cudaStream_t stream)
        {
            CV_Assert(levels[0].type() == LEVEL_TYPE_CODE && levels[0].rows == 1);
            CV_Assert(levels[1].type() == LEVEL_TYPE_CODE && levels[1].rows == 1);
            CV_Assert(levels[2].type() == LEVEL_TYPE_CODE && levels[2].rows == 1);
            CV_Assert(levels[3].type() == LEVEL_TYPE_CODE && levels[3].rows == 1);

            hist[0].create(1, levels[0].cols - 1, CV_32S);
            hist[1].create(1, levels[1].cols - 1, CV_32S);
            hist[2].create(1, levels[2].cols - 1, CV_32S);
            hist[3].create(1, levels[3].cols - 1, CV_32S);

            Npp32s* pHist[] = {hist[0].ptr<Npp32s>(), hist[1].ptr<Npp32s>(), hist[2].ptr<Npp32s>(), hist[3].ptr<Npp32s>()};
            int nLevels[] = {levels[0].cols, levels[1].cols, levels[2].cols, levels[3].cols};
            const level_t* pLevels[] = {levels[0].ptr<level_t>(), levels[1].ptr<level_t>(), levels[2].ptr<level_t>(), levels[3].ptr<level_t>()};

            NppiSize sz;
            sz.width = src.cols;
            sz.height = src.rows;

            int buf_size;
            get_buf_size(sz, nLevels, &buf_size);

            ensureSizeIsEnough(1, buf_size, CV_8U, buffer);

            NppStreamHandler h(stream);

            nppSafeCall( func(src.ptr<src_t>(), static_cast<int>(src.step), sz, pHist, pLevels, nLevels, buffer.ptr<Npp8u>()) );

            if (stream == 0)
                cudaSafeCall( cudaDeviceSynchronize() );
        }
    };
}

void cv::gpu::evenLevels(GpuMat& levels, int nLevels, int lowerLevel, int upperLevel)
{
    Mat host_levels(1, nLevels, CV_32SC1);
    nppSafeCall( nppiEvenLevelsHost_32s(host_levels.ptr<Npp32s>(), nLevels, lowerLevel, upperLevel) );
    levels.upload(host_levels);
}

void cv::gpu::histEven(const GpuMat& src, GpuMat& hist, int histSize, int lowerLevel, int upperLevel, Stream& stream)
{
    GpuMat buf;
    histEven(src, hist, buf, histSize, lowerLevel, upperLevel, stream);
}

void cv::gpu::histEven(const GpuMat& src, GpuMat& hist, GpuMat& buf, int histSize, int lowerLevel, int upperLevel, Stream& stream)
{
    CV_Assert(src.type() == CV_8UC1 || src.type() == CV_16UC1 || src.type() == CV_16SC1 );

    typedef void (*hist_t)(const GpuMat& src, GpuMat& hist, GpuMat& buf, int levels, int lowerLevel, int upperLevel, cudaStream_t stream);
    static const hist_t hist_callers[] =
    {
        NppHistogramEvenC1<CV_8U , nppiHistogramEven_8u_C1R , nppiHistogramEvenGetBufferSize_8u_C1R >::hist,
        0,
        NppHistogramEvenC1<CV_16U, nppiHistogramEven_16u_C1R, nppiHistogramEvenGetBufferSize_16u_C1R>::hist,
        NppHistogramEvenC1<CV_16S, nppiHistogramEven_16s_C1R, nppiHistogramEvenGetBufferSize_16s_C1R>::hist
    };

    hist_callers[src.depth()](src, hist, buf, histSize, lowerLevel, upperLevel, StreamAccessor::getStream(stream));
}

void cv::gpu::histEven(const GpuMat& src, GpuMat hist[4], int histSize[4], int lowerLevel[4], int upperLevel[4], Stream& stream)
{
    GpuMat buf;
    histEven(src, hist, buf, histSize, lowerLevel, upperLevel, stream);
}

void cv::gpu::histEven(const GpuMat& src, GpuMat hist[4], GpuMat& buf, int histSize[4], int lowerLevel[4], int upperLevel[4], Stream& stream)
{
    CV_Assert(src.type() == CV_8UC4 || src.type() == CV_16UC4 || src.type() == CV_16SC4 );

    typedef void (*hist_t)(const GpuMat& src, GpuMat hist[4], GpuMat& buf, int levels[4], int lowerLevel[4], int upperLevel[4], cudaStream_t stream);
    static const hist_t hist_callers[] =
    {
        NppHistogramEvenC4<CV_8U , nppiHistogramEven_8u_C4R , nppiHistogramEvenGetBufferSize_8u_C4R >::hist,
        0,
        NppHistogramEvenC4<CV_16U, nppiHistogramEven_16u_C4R, nppiHistogramEvenGetBufferSize_16u_C4R>::hist,
        NppHistogramEvenC4<CV_16S, nppiHistogramEven_16s_C4R, nppiHistogramEvenGetBufferSize_16s_C4R>::hist
    };

    hist_callers[src.depth()](src, hist, buf, histSize, lowerLevel, upperLevel, StreamAccessor::getStream(stream));
}

void cv::gpu::histRange(const GpuMat& src, GpuMat& hist, const GpuMat& levels, Stream& stream)
{
    GpuMat buf;
    histRange(src, hist, levels, buf, stream);
}

void cv::gpu::histRange(const GpuMat& src, GpuMat& hist, const GpuMat& levels, GpuMat& buf, Stream& stream)
{
    CV_Assert(src.type() == CV_8UC1 || src.type() == CV_16UC1 || src.type() == CV_16SC1 || src.type() == CV_32FC1);

    typedef void (*hist_t)(const GpuMat& src, GpuMat& hist, const GpuMat& levels, GpuMat& buf, cudaStream_t stream);
    static const hist_t hist_callers[] =
    {
        NppHistogramRangeC1<CV_8U , nppiHistogramRange_8u_C1R , nppiHistogramRangeGetBufferSize_8u_C1R >::hist,
        0,
        NppHistogramRangeC1<CV_16U, nppiHistogramRange_16u_C1R, nppiHistogramRangeGetBufferSize_16u_C1R>::hist,
        NppHistogramRangeC1<CV_16S, nppiHistogramRange_16s_C1R, nppiHistogramRangeGetBufferSize_16s_C1R>::hist,
        0,
        NppHistogramRangeC1<CV_32F, nppiHistogramRange_32f_C1R, nppiHistogramRangeGetBufferSize_32f_C1R>::hist
    };

    hist_callers[src.depth()](src, hist, levels, buf, StreamAccessor::getStream(stream));
}

void cv::gpu::histRange(const GpuMat& src, GpuMat hist[4], const GpuMat levels[4], Stream& stream)
{
    GpuMat buf;
    histRange(src, hist, levels, buf, stream);
}

void cv::gpu::histRange(const GpuMat& src, GpuMat hist[4], const GpuMat levels[4], GpuMat& buf, Stream& stream)
{
    CV_Assert(src.type() == CV_8UC4 || src.type() == CV_16UC4 || src.type() == CV_16SC4 || src.type() == CV_32FC4);

    typedef void (*hist_t)(const GpuMat& src, GpuMat hist[4], const GpuMat levels[4], GpuMat& buf, cudaStream_t stream);
    static const hist_t hist_callers[] =
    {
        NppHistogramRangeC4<CV_8U , nppiHistogramRange_8u_C4R , nppiHistogramRangeGetBufferSize_8u_C4R >::hist,
        0,
        NppHistogramRangeC4<CV_16U, nppiHistogramRange_16u_C4R, nppiHistogramRangeGetBufferSize_16u_C4R>::hist,
        NppHistogramRangeC4<CV_16S, nppiHistogramRange_16s_C4R, nppiHistogramRangeGetBufferSize_16s_C4R>::hist,
        0,
        NppHistogramRangeC4<CV_32F, nppiHistogramRange_32f_C4R, nppiHistogramRangeGetBufferSize_32f_C4R>::hist
    };

    hist_callers[src.depth()](src, hist, levels, buf, StreamAccessor::getStream(stream));
}

namespace hist
{
    void histogram256(PtrStepSzb src, int* hist, cudaStream_t stream);
    void equalizeHist(PtrStepSzb src, PtrStepSzb dst, const int* lut, cudaStream_t stream);
}

void cv::gpu::calcHist(const GpuMat& src, GpuMat& hist, Stream& stream)
{
    CV_Assert(src.type() == CV_8UC1);

    hist.create(1, 256, CV_32SC1);
    hist.setTo(Scalar::all(0));

    hist::histogram256(src, hist.ptr<int>(), StreamAccessor::getStream(stream));
}

void cv::gpu::equalizeHist(const GpuMat& src, GpuMat& dst, Stream& stream)
{
    GpuMat hist;
    GpuMat buf;
    equalizeHist(src, dst, hist, buf, stream);
}

void cv::gpu::equalizeHist(const GpuMat& src, GpuMat& dst, GpuMat& hist, GpuMat& buf, Stream& s)
{
    CV_Assert(src.type() == CV_8UC1);

    dst.create(src.size(), src.type());

    int intBufSize;
    nppSafeCall( nppsIntegralGetBufferSize_32s(256, &intBufSize) );

    ensureSizeIsEnough(1, intBufSize + 256 * sizeof(int), CV_8UC1, buf);

    GpuMat intBuf(1, intBufSize, CV_8UC1, buf.ptr());
    GpuMat lut(1, 256, CV_32S, buf.ptr() + intBufSize);

    calcHist(src, hist, s);

    cudaStream_t stream = StreamAccessor::getStream(s);

    NppStreamHandler h(stream);

    nppSafeCall( nppsIntegral_32s(hist.ptr<Npp32s>(), lut.ptr<Npp32s>(), 256, intBuf.ptr<Npp8u>()) );

    hist::equalizeHist(src, dst, lut.ptr<int>(), stream);
}

////////////////////////////////////////////////////////////////////////
// cornerHarris & minEgenVal

namespace cv { namespace gpu { namespace cudev
{
    namespace imgproc
    {
        void cornerHarris_gpu(int block_size, float k, PtrStepSzf Dx, PtrStepSzf Dy, PtrStepSzf dst, int border_type, cudaStream_t stream);
        void cornerMinEigenVal_gpu(int block_size, PtrStepSzf Dx, PtrStepSzf Dy, PtrStepSzf dst, int border_type, cudaStream_t stream);
    }
}}}

namespace
{
    void extractCovData(const GpuMat& src, GpuMat& Dx, GpuMat& Dy, GpuMat& buf, int blockSize, int ksize, int borderType, Stream& stream)
    {
        double scale = static_cast<double>(1 << ((ksize > 0 ? ksize : 3) - 1)) * blockSize;

        if (ksize < 0)
            scale *= 2.;

        if (src.depth() == CV_8U)
            scale *= 255.;

        scale = 1./scale;

        Dx.create(src.size(), CV_32F);
        Dy.create(src.size(), CV_32F);

        if (ksize > 0)
        {
            Sobel(src, Dx, CV_32F, 1, 0, buf, ksize, scale, borderType, -1, stream);
            Sobel(src, Dy, CV_32F, 0, 1, buf, ksize, scale, borderType, -1, stream);
        }
        else
        {
            Scharr(src, Dx, CV_32F, 1, 0, buf, scale, borderType, -1, stream);
            Scharr(src, Dy, CV_32F, 0, 1, buf, scale, borderType, -1, stream);
        }
    }
}

void cv::gpu::cornerHarris(const GpuMat& src, GpuMat& dst, int blockSize, int ksize, double k, int borderType)
{
    GpuMat Dx, Dy;
    cornerHarris(src, dst, Dx, Dy, blockSize, ksize, k, borderType);
}

void cv::gpu::cornerHarris(const GpuMat& src, GpuMat& dst, GpuMat& Dx, GpuMat& Dy, int blockSize, int ksize, double k, int borderType)
{
    GpuMat buf;
    cornerHarris(src, dst, Dx, Dy, buf, blockSize, ksize, k, borderType);
}

void cv::gpu::cornerHarris(const GpuMat& src, GpuMat& dst, GpuMat& Dx, GpuMat& Dy, GpuMat& buf, int blockSize, int ksize, double k, int borderType, Stream& stream)
{
    using namespace cv::gpu::cudev::imgproc;

    CV_Assert(borderType == cv::BORDER_REFLECT101 || borderType == cv::BORDER_REPLICATE || borderType == cv::BORDER_REFLECT);

    int gpuBorderType;
    CV_Assert(tryConvertToGpuBorderType(borderType, gpuBorderType));

    extractCovData(src, Dx, Dy, buf, blockSize, ksize, borderType, stream);

    dst.create(src.size(), CV_32F);

    cornerHarris_gpu(blockSize, static_cast<float>(k), Dx, Dy, dst, gpuBorderType, StreamAccessor::getStream(stream));
}

void cv::gpu::cornerMinEigenVal(const GpuMat& src, GpuMat& dst, int blockSize, int ksize, int borderType)
{
    GpuMat Dx, Dy;
    cornerMinEigenVal(src, dst, Dx, Dy, blockSize, ksize, borderType);
}

void cv::gpu::cornerMinEigenVal(const GpuMat& src, GpuMat& dst, GpuMat& Dx, GpuMat& Dy, int blockSize, int ksize, int borderType)
{
    GpuMat buf;
    cornerMinEigenVal(src, dst, Dx, Dy, buf, blockSize, ksize, borderType);
}

void cv::gpu::cornerMinEigenVal(const GpuMat& src, GpuMat& dst, GpuMat& Dx, GpuMat& Dy, GpuMat& buf, int blockSize, int ksize, int borderType, Stream& stream)
{
    using namespace ::cv::gpu::cudev::imgproc;

    CV_Assert(borderType == cv::BORDER_REFLECT101 || borderType == cv::BORDER_REPLICATE || borderType == cv::BORDER_REFLECT);

    int gpuBorderType;
    CV_Assert(tryConvertToGpuBorderType(borderType, gpuBorderType));

    extractCovData(src, Dx, Dy, buf, blockSize, ksize, borderType, stream);

    dst.create(src.size(), CV_32F);

    cornerMinEigenVal_gpu(blockSize, Dx, Dy, dst, gpuBorderType, StreamAccessor::getStream(stream));
}


//////////////////////////////////////////////////////////////////////////////
// Canny

void cv::gpu::CannyBuf::create(const Size& image_size, int apperture_size)
{
    if (apperture_size > 0)
    {
        ensureSizeIsEnough(image_size, CV_32SC1, dx);
        ensureSizeIsEnough(image_size, CV_32SC1, dy);

        if (apperture_size != 3)
        {
            filterDX = createDerivFilter_GPU(CV_8UC1, CV_32S, 1, 0, apperture_size, BORDER_REPLICATE);
            filterDY = createDerivFilter_GPU(CV_8UC1, CV_32S, 0, 1, apperture_size, BORDER_REPLICATE);
        }
    }

    ensureSizeIsEnough(image_size, CV_32FC1, mag);
    ensureSizeIsEnough(image_size, CV_32SC1, map);

    ensureSizeIsEnough(1, image_size.area(), CV_16UC2, st1);
    ensureSizeIsEnough(1, image_size.area(), CV_16UC2, st2);
}

void cv::gpu::CannyBuf::release()
{
    dx.release();
    dy.release();
    mag.release();
    map.release();
    st1.release();
    st2.release();
}

namespace canny
{
    void calcMagnitude(PtrStepSzb srcWhole, int xoff, int yoff, PtrStepSzi dx, PtrStepSzi dy, PtrStepSzf mag, bool L2Grad);
    void calcMagnitude(PtrStepSzi dx, PtrStepSzi dy, PtrStepSzf mag, bool L2Grad);

    void calcMap(PtrStepSzi dx, PtrStepSzi dy, PtrStepSzf mag, PtrStepSzi map, float low_thresh, float high_thresh);

    void edgesHysteresisLocal(PtrStepSzi map, ushort2* st1);

    void edgesHysteresisGlobal(PtrStepSzi map, ushort2* st1, ushort2* st2);

    void getEdges(PtrStepSzi map, PtrStepSzb dst);
}

namespace
{
    void CannyCaller(const GpuMat& dx, const GpuMat& dy, CannyBuf& buf, GpuMat& dst, float low_thresh, float high_thresh)
    {
        using namespace canny;

        buf.map.setTo(Scalar::all(0));
        calcMap(dx, dy, buf.mag, buf.map, low_thresh, high_thresh);

        edgesHysteresisLocal(buf.map, buf.st1.ptr<ushort2>());

        edgesHysteresisGlobal(buf.map, buf.st1.ptr<ushort2>(), buf.st2.ptr<ushort2>());

        getEdges(buf.map, dst);
    }
}

void cv::gpu::Canny(const GpuMat& src, GpuMat& dst, double low_thresh, double high_thresh, int apperture_size, bool L2gradient)
{
    CannyBuf buf;
    Canny(src, buf, dst, low_thresh, high_thresh, apperture_size, L2gradient);
}

void cv::gpu::Canny(const GpuMat& src, CannyBuf& buf, GpuMat& dst, double low_thresh, double high_thresh, int apperture_size, bool L2gradient)
{
    using namespace canny;

    CV_Assert(src.type() == CV_8UC1);

    if (!deviceSupports(SHARED_ATOMICS))
        CV_Error(cv::Error::StsNotImplemented, "The device doesn't support shared atomics");

    if( low_thresh > high_thresh )
        std::swap( low_thresh, high_thresh);

    dst.create(src.size(), CV_8U);
    buf.create(src.size(), apperture_size);

    if (apperture_size == 3)
    {
        Size wholeSize;
        Point ofs;
        src.locateROI(wholeSize, ofs);
        GpuMat srcWhole(wholeSize, src.type(), src.datastart, src.step);

        calcMagnitude(srcWhole, ofs.x, ofs.y, buf.dx, buf.dy, buf.mag, L2gradient);
    }
    else
    {
        buf.filterDX->apply(src, buf.dx, Rect(0, 0, src.cols, src.rows));
        buf.filterDY->apply(src, buf.dy, Rect(0, 0, src.cols, src.rows));

        calcMagnitude(buf.dx, buf.dy, buf.mag, L2gradient);
    }

    CannyCaller(buf.dx, buf.dy, buf, dst, static_cast<float>(low_thresh), static_cast<float>(high_thresh));
}

void cv::gpu::Canny(const GpuMat& dx, const GpuMat& dy, GpuMat& dst, double low_thresh, double high_thresh, bool L2gradient)
{
    CannyBuf buf;
    Canny(dx, dy, buf, dst, low_thresh, high_thresh, L2gradient);
}

void cv::gpu::Canny(const GpuMat& dx, const GpuMat& dy, CannyBuf& buf, GpuMat& dst, double low_thresh, double high_thresh, bool L2gradient)
{
    using namespace canny;

    CV_Assert(TargetArchs::builtWith(SHARED_ATOMICS) && DeviceInfo().supports(SHARED_ATOMICS));
    CV_Assert(dx.type() == CV_32SC1 && dy.type() == CV_32SC1 && dx.size() == dy.size());

    if( low_thresh > high_thresh )
        std::swap( low_thresh, high_thresh);

    dst.create(dx.size(), CV_8U);
    buf.create(dx.size(), -1);

    calcMagnitude(dx, dy, buf.mag, L2gradient);

    CannyCaller(dx, dy, buf, dst, static_cast<float>(low_thresh), static_cast<float>(high_thresh));
}

////////////////////////////////////////////////////////////////////////
// CLAHE

namespace clahe
{
    void calcLut(PtrStepSzb src, PtrStepb lut, int tilesX, int tilesY, int2 tileSize, int clipLimit, float lutScale, cudaStream_t stream);
    void transform(PtrStepSzb src, PtrStepSzb dst, PtrStepb lut, int tilesX, int tilesY, int2 tileSize, cudaStream_t stream);
}

namespace
{
    class CLAHE_Impl : public cv::gpu::CLAHE
    {
    public:
        CLAHE_Impl(double clipLimit = 40.0, int tilesX = 8, int tilesY = 8);

        cv::AlgorithmInfo* info() const;

        void apply(cv::InputArray src, cv::OutputArray dst);
        void apply(InputArray src, OutputArray dst, Stream& stream);

        void setClipLimit(double clipLimit);
        double getClipLimit() const;

        void setTilesGridSize(cv::Size tileGridSize);
        cv::Size getTilesGridSize() const;

        void collectGarbage();

    private:
        double clipLimit_;
        int tilesX_;
        int tilesY_;

        GpuMat srcExt_;
        GpuMat lut_;
    };

    CLAHE_Impl::CLAHE_Impl(double clipLimit, int tilesX, int tilesY) :
        clipLimit_(clipLimit), tilesX_(tilesX), tilesY_(tilesY)
    {
    }

    CV_INIT_ALGORITHM(CLAHE_Impl, "CLAHE_GPU",
        obj.info()->addParam(obj, "clipLimit", obj.clipLimit_);
        obj.info()->addParam(obj, "tilesX", obj.tilesX_);
        obj.info()->addParam(obj, "tilesY", obj.tilesY_))

    void CLAHE_Impl::apply(cv::InputArray _src, cv::OutputArray _dst)
    {
        apply(_src, _dst, Stream::Null());
    }

    void CLAHE_Impl::apply(InputArray _src, OutputArray _dst, Stream& s)
    {
        GpuMat src = _src.getGpuMat();

        CV_Assert( src.type() == CV_8UC1 );

        _dst.create( src.size(), src.type() );
        GpuMat dst = _dst.getGpuMat();

        const int histSize = 256;

        ensureSizeIsEnough(tilesX_ * tilesY_, histSize, CV_8UC1, lut_);

        cudaStream_t stream = StreamAccessor::getStream(s);

        cv::Size tileSize;
        GpuMat srcForLut;

        if (src.cols % tilesX_ == 0 && src.rows % tilesY_ == 0)
        {
            tileSize = cv::Size(src.cols / tilesX_, src.rows / tilesY_);
            srcForLut = src;
        }
        else
        {
            cv::gpu::copyMakeBorder(src, srcExt_, 0, tilesY_ - (src.rows % tilesY_), 0, tilesX_ - (src.cols % tilesX_), cv::BORDER_REFLECT_101, cv::Scalar(), s);

            tileSize = cv::Size(srcExt_.cols / tilesX_, srcExt_.rows / tilesY_);
            srcForLut = srcExt_;
        }

        const int tileSizeTotal = tileSize.area();
        const float lutScale = static_cast<float>(histSize - 1) / tileSizeTotal;

        int clipLimit = 0;
        if (clipLimit_ > 0.0)
        {
            clipLimit = static_cast<int>(clipLimit_ * tileSizeTotal / histSize);
            clipLimit = std::max(clipLimit, 1);
        }

        clahe::calcLut(srcForLut, lut_, tilesX_, tilesY_, make_int2(tileSize.width, tileSize.height), clipLimit, lutScale, stream);

        clahe::transform(src, dst, lut_, tilesX_, tilesY_, make_int2(tileSize.width, tileSize.height), stream);
    }

    void CLAHE_Impl::setClipLimit(double clipLimit)
    {
        clipLimit_ = clipLimit;
    }

    double CLAHE_Impl::getClipLimit() const
    {
        return clipLimit_;
    }

    void CLAHE_Impl::setTilesGridSize(cv::Size tileGridSize)
    {
        tilesX_ = tileGridSize.width;
        tilesY_ = tileGridSize.height;
    }

    cv::Size CLAHE_Impl::getTilesGridSize() const
    {
        return cv::Size(tilesX_, tilesY_);
    }

    void CLAHE_Impl::collectGarbage()
    {
        srcExt_.release();
        lut_.release();
    }
}

cv::Ptr<cv::gpu::CLAHE> cv::gpu::createCLAHE(double clipLimit, cv::Size tileGridSize)
{
    return new CLAHE_Impl(clipLimit, tileGridSize.width, tileGridSize.height);
}

////////////////////////////////////////////////////////////////////////
// alphaComp

namespace
{
    template <int DEPTH> struct NppAlphaCompFunc
    {
        typedef typename NppTypeTraits<DEPTH>::npp_t npp_t;

        typedef NppStatus (*func_t)(const npp_t* pSrc1, int nSrc1Step, const npp_t* pSrc2, int nSrc2Step, npp_t* pDst, int nDstStep, NppiSize oSizeROI, NppiAlphaOp eAlphaOp);
    };

    template <int DEPTH, typename NppAlphaCompFunc<DEPTH>::func_t func> struct NppAlphaComp
    {
        typedef typename NppTypeTraits<DEPTH>::npp_t npp_t;

        static void call(const GpuMat& img1, const GpuMat& img2, GpuMat& dst, NppiAlphaOp eAlphaOp, cudaStream_t stream)
        {
            NppStreamHandler h(stream);

            NppiSize oSizeROI;
            oSizeROI.width = img1.cols;
            oSizeROI.height = img2.rows;

            nppSafeCall( func(img1.ptr<npp_t>(), static_cast<int>(img1.step), img2.ptr<npp_t>(), static_cast<int>(img2.step),
                              dst.ptr<npp_t>(), static_cast<int>(dst.step), oSizeROI, eAlphaOp) );

            if (stream == 0)
                cudaSafeCall( cudaDeviceSynchronize() );
        }
    };
}

void cv::gpu::alphaComp(const GpuMat& img1, const GpuMat& img2, GpuMat& dst, int alpha_op, Stream& stream)
{
    static const NppiAlphaOp npp_alpha_ops[] = {
        NPPI_OP_ALPHA_OVER,
        NPPI_OP_ALPHA_IN,
        NPPI_OP_ALPHA_OUT,
        NPPI_OP_ALPHA_ATOP,
        NPPI_OP_ALPHA_XOR,
        NPPI_OP_ALPHA_PLUS,
        NPPI_OP_ALPHA_OVER_PREMUL,
        NPPI_OP_ALPHA_IN_PREMUL,
        NPPI_OP_ALPHA_OUT_PREMUL,
        NPPI_OP_ALPHA_ATOP_PREMUL,
        NPPI_OP_ALPHA_XOR_PREMUL,
        NPPI_OP_ALPHA_PLUS_PREMUL,
        NPPI_OP_ALPHA_PREMUL
    };

    typedef void (*func_t)(const GpuMat& img1, const GpuMat& img2, GpuMat& dst, NppiAlphaOp eAlphaOp, cudaStream_t stream);

    static const func_t funcs[] =
    {
        NppAlphaComp<CV_8U, nppiAlphaComp_8u_AC4R>::call,
        0,
        NppAlphaComp<CV_16U, nppiAlphaComp_16u_AC4R>::call,
        0,
        NppAlphaComp<CV_32S, nppiAlphaComp_32s_AC4R>::call,
        NppAlphaComp<CV_32F, nppiAlphaComp_32f_AC4R>::call
    };

    CV_Assert( img1.type() == CV_8UC4 || img1.type() == CV_16UC4 || img1.type() == CV_32SC4 || img1.type() == CV_32FC4 );
    CV_Assert( img1.size() == img2.size() && img1.type() == img2.type() );

    dst.create(img1.size(), img1.type());

    const func_t func = funcs[img1.depth()];

    func(img1, img2, dst, npp_alpha_ops[alpha_op], StreamAccessor::getStream(stream));
}

#endif /* !defined (HAVE_CUDA) */
