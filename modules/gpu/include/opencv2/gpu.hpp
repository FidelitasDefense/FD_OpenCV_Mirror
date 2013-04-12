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

#ifndef __OPENCV_GPU_HPP__
#define __OPENCV_GPU_HPP__

#ifndef SKIP_INCLUDES
#include <vector>
#include <memory>
#include <iosfwd>
#endif

#include "opencv2/core/gpumat.hpp"
#include "opencv2/gpuarithm.hpp"
#include "opencv2/gpufilters.hpp"
#include "opencv2/gpuimgproc.hpp"
#include "opencv2/gpufeatures2d.hpp"
#include "opencv2/gpuvideo.hpp"
#include "opencv2/gpucalib3d.hpp"
#include "opencv2/gpuobjdetect.hpp"

#include "opencv2/imgproc.hpp"
#include "opencv2/objdetect.hpp"
#include "opencv2/features2d.hpp"

namespace cv { namespace gpu {
////////////////////////////// Image processing //////////////////////////////




///////////////////////////// Calibration 3D //////////////////////////////////

//////////////////////////////// Image Labeling ////////////////////////////////



////////////////////////////////// Histograms //////////////////////////////////



//////////////////////////////// StereoBM_GPU ////////////////////////////////



////////////////////////// StereoBeliefPropagation ///////////////////////////


/////////////////////////// StereoConstantSpaceBP ///////////////////////////



/////////////////////////// DisparityBilateralFilter ///////////////////////////






////////////////////////////////// BruteForceMatcher //////////////////////////////////



template <class Distance>
class CV_EXPORTS BruteForceMatcher_GPU;

template <typename T>
class CV_EXPORTS BruteForceMatcher_GPU< L1<T> > : public BFMatcher_GPU
{
public:
    explicit BruteForceMatcher_GPU() : BFMatcher_GPU(NORM_L1) {}
    explicit BruteForceMatcher_GPU(L1<T> /*d*/) : BFMatcher_GPU(NORM_L1) {}
};
template <typename T>
class CV_EXPORTS BruteForceMatcher_GPU< L2<T> > : public BFMatcher_GPU
{
public:
    explicit BruteForceMatcher_GPU() : BFMatcher_GPU(NORM_L2) {}
    explicit BruteForceMatcher_GPU(L2<T> /*d*/) : BFMatcher_GPU(NORM_L2) {}
};
template <> class CV_EXPORTS BruteForceMatcher_GPU< Hamming > : public BFMatcher_GPU
{
public:
    explicit BruteForceMatcher_GPU() : BFMatcher_GPU(NORM_HAMMING) {}
    explicit BruteForceMatcher_GPU(Hamming /*d*/) : BFMatcher_GPU(NORM_HAMMING) {}
};

////////////////////////////////// CascadeClassifier_GPU //////////////////////////////////////////


////////////////////////////////// FAST //////////////////////////////////////////



////////////////////////////////// ORB //////////////////////////////////////////





















//! removes points (CV_32FC2, single row matrix) with zero mask value
CV_EXPORTS void compactPoints(GpuMat &points0, GpuMat &points1, const GpuMat &mask);

CV_EXPORTS void calcWobbleSuppressionMaps(
        int left, int idx, int right, Size size, const Mat &ml, const Mat &mr,
        GpuMat &mapx, GpuMat &mapy);

} // namespace gpu

} // namespace cv

#endif /* __OPENCV_GPU_HPP__ */
