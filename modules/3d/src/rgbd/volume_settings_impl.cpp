// This file is part of OpenCV project.
// It is subject to the license terms in the LICENSE file found in the top-level directory
// of this distribution and at http://opencv.org/license.html


#include "../precomp.hpp"

namespace cv
{

Vec4i calcVolumeDimentions(Point3i volumeResolution, bool ZFirstMemOrder);

class VolumeSettings::Impl
{
public:
    Impl() {};
    virtual ~Impl() {};

    virtual void  setWidth(int  val) = 0;
    virtual int   getWidth() const = 0;
    virtual void  setHeight(int  val) = 0;
    virtual int   getHeight() const = 0;
    virtual void  setDepthFactor(float val) = 0;
    virtual float getDepthFactor() const = 0;
    virtual void  setVoxelSize(float  val) = 0;
    virtual float getVoxelSize() const = 0;
    virtual void  setTsdfTruncateDistance(float val) = 0;
    virtual float getTsdfTruncateDistance() const = 0;
    virtual void  setMaxDepth(float val) = 0;
    virtual float getMaxDepth() const = 0;
    virtual void  setMaxWeight(int val) = 0;
    virtual int   getMaxWeight() const = 0;
    virtual void  setRaycastStepFactor(float val) = 0;
    virtual float getRaycastStepFactor() const = 0;

    virtual void setVolumePose(InputArray val) = 0;
    virtual void getVolumePose(OutputArray val) const = 0;
    virtual void setVolumeResolution(InputArray val) = 0;
    virtual void getVolumeResolution(OutputArray val) const = 0;
    virtual void getVolumeDimentions(OutputArray val) const = 0;
    virtual void setCameraIntegrateIntrinsics(InputArray val) = 0;
    virtual void getCameraIntegrateIntrinsics(OutputArray val) const = 0;
    virtual void setRGBCameraIntegrateIntrinsics(InputArray val) = 0;
    virtual void getRGBCameraIntegrateIntrinsics(OutputArray val) const = 0;
    virtual void setCameraRaycastIntrinsics(InputArray val) = 0;
    virtual void getCameraRaycastIntrinsics(OutputArray val) const = 0;
    virtual void setRGBCameraRaycastIntrinsics(InputArray val) = 0;
    virtual void getRGBCameraRaycastIntrinsics(OutputArray val) const = 0;
};

class VolumeSettingsImpl : public VolumeSettings::Impl
{
public:
    VolumeSettingsImpl();
    VolumeSettingsImpl(VolumeType volumeType);
    ~VolumeSettingsImpl();

    virtual void  setWidth(int  val) override;
    virtual int   getWidth() const override;
    virtual void  setHeight(int  val) override;
    virtual int   getHeight() const override;
    virtual void  setDepthFactor(float val) override;
    virtual float getDepthFactor() const override;
    virtual void  setVoxelSize(float  val) override;
    virtual float getVoxelSize() const override;
    virtual void  setTsdfTruncateDistance(float val) override;
    virtual float getTsdfTruncateDistance() const override;
    virtual void  setMaxDepth(float val) override;
    virtual float getMaxDepth() const override;
    virtual void  setMaxWeight(int val) override;
    virtual int   getMaxWeight() const override;
    virtual void  setRaycastStepFactor(float val) override;
    virtual float getRaycastStepFactor() const override;

    virtual void setVolumePose(InputArray val) override;
    virtual void getVolumePose(OutputArray val) const override;
    virtual void setVolumeResolution(InputArray val) override;
    virtual void getVolumeResolution(OutputArray val) const override;
    virtual void getVolumeDimentions(OutputArray val) const override;
    virtual void setCameraIntegrateIntrinsics(InputArray val) override;
    virtual void getCameraIntegrateIntrinsics(OutputArray val) const override;
    virtual void setRGBCameraIntegrateIntrinsics(InputArray val) override;
    virtual void getRGBCameraIntegrateIntrinsics(OutputArray val) const override;
    virtual void setCameraRaycastIntrinsics(InputArray val) override;
    virtual void getCameraRaycastIntrinsics(OutputArray val) const override;
    virtual void setRGBCameraRaycastIntrinsics(InputArray val) override;
    virtual void getRGBCameraRaycastIntrinsics(OutputArray val) const override;

private:
    VolumeType volumeType;

    int   width;
    int   height;
    float depthFactor;
    float voxelSize;
    float tsdfTruncateDistance;
    float maxDepth;
    int   maxWeight;
    float raycastStepFactor;
    bool  zFirstMemOrder;

    Matx44f volumePose;
    Point3i volumeResolution;
    Vec4i volumeDimentions;
    Matx33f cameraIntegrateIntrinsics;
    Matx33f rgb_cameraIntegrateIntrinsics;
    Matx33f cameraRaycastIntrinsics;
    Matx33f rgb_cameraRaycastIntrinsics;

public:
    // duplicate classes for all volumes

    class DefaultTsdfSets {
    public:
        static const int width  = 640;
        static const int height = 480;
        static constexpr float fx = 525.f; // focus point x axis
        static constexpr float fy = 525.f; // focus point y axis
        static constexpr float cx = float(width) / 2.f - 0.5f;  // central point x axis
        static constexpr float cy = float(height) / 2.f - 0.5f; // central point y axis
        static constexpr float depthFactor = 5000.f; // 5000 for the 16-bit PNG files, 1 for the 32-bit float images in the ROS bag files
        static constexpr float volumeSize = 3.f; // meters
        static constexpr float voxelSize = volumeSize / 128.f; //meters
        static constexpr float tsdfTruncateDistance = 2 * voxelSize;
        static constexpr float maxDepth = 0.f;
        static const int maxWeight = 64; // number of frames
        static constexpr float raycastStepFactor = 0.75f;
        static const bool zFirstMemOrder = true; // order of voxels in volume

        const Matx33f  cameraIntegrateIntrinsics = Matx33f(fx, 0, cx, 0, fy, cy, 0, 0, 1); // camera settings
        const Matx33f  rgb_cameraIntegrateIntrinsics = Matx33f(fx, 0, cx, 0, fy, cy, 0, 0, 1); // camera settings
        const Matx33f  cameraRaycastIntrinsics = Matx33f(fx, 0, cx, 0, fy, cy, 0, 0, 1); // camera settings
        const Matx33f  rgb_cameraRaycastIntrinsics = Matx33f(fx, 0, cx, 0, fy, cy, 0, 0, 1); // camera settings
        const Affine3f volumePose = Affine3f().translate(Vec3f(-volumeSize / 2.f, -volumeSize / 2.f, 0.5f));
        const Matx44f  volumePoseMatrix = volumePose.matrix;
        // Unlike original code, this should work with any volume size
        // Not only when (x,y,z % 32) == 0
        const Point3i  volumeResolution = Vec3i::all(128); //number of voxels
    };

    class DefaultHashTsdfSets {
    public:
        static const int width = 640;
        static const int height = 480;
        static constexpr float fx = 525.f; // focus point x axis
        static constexpr float fy = 525.f; // focus point y axis
        static constexpr float cx = float(width) / 2.f - 0.5f;  // central point x axis
        static constexpr float cy = float(height) / 2.f - 0.5f; // central point y axis
        static constexpr float depthFactor = 5000.f; // 5000 for the 16-bit PNG files, 1 for the 32-bit float images in the ROS bag files
        static constexpr float volumeSize = 3.f; // meters
        static constexpr float voxelSize = volumeSize / 512.f; //meters
        static constexpr float tsdfTruncateDistance = 7 * voxelSize;
        static constexpr float maxDepth = 4.f;
        static const int maxWeight = 64; // number of frames
        static constexpr float raycastStepFactor = 0.25f;
        static const bool zFirstMemOrder = true; // order of voxels in volume

        const Matx33f  cameraIntegrateIntrinsics = Matx33f(fx, 0, cx, 0, fy, cy, 0, 0, 1); // camera settings
        const Matx33f  rgb_cameraIntegrateIntrinsics = Matx33f(fx, 0, cx, 0, fy, cy, 0, 0, 1); // camera settings
        const Matx33f  cameraRaycastIntrinsics = Matx33f(fx, 0, cx, 0, fy, cy, 0, 0, 1); // camera settings
        const Matx33f  rgb_cameraRaycastIntrinsics = Matx33f(fx, 0, cx, 0, fy, cy, 0, 0, 1); // camera settings
        const Affine3f volumePose = Affine3f().translate(Vec3f(-volumeSize / 2.f, -volumeSize / 2.f, 0.5f));
        const Matx44f  volumePoseMatrix = volumePose.matrix;
        // Unlike original code, this should work with any volume size
        // Not only when (x,y,z % 32) == 0
        const Point3i  volumeResolution = Vec3i::all(16); //number of voxels
    };

    class DefaultColorTsdfSets {
    public:
        static const int width = 640;
        static const int height = 480;
        static constexpr float fx = 525.f; // focus point x axis
        static constexpr float fy = 525.f; // focus point y axis
        static constexpr float cx = float(width) / 2.f - 0.5f;  // central point x axis
        static constexpr float cy = float(height) / 2.f - 0.5f; // central point y axis
        static constexpr float rgb_fx = 525.f; // focus point x axis
        static constexpr float rgb_fy = 525.f; // focus point y axis
        static constexpr float rgb_cx = float(width) / 2.f - 0.5f;  // central point x axis
        static constexpr float rgb_cy = float(height) / 2.f - 0.5f; // central point y axis
        static constexpr float depthFactor = 5000.f; // 5000 for the 16-bit PNG files, 1 for the 32-bit float images in the ROS bag files
        static constexpr float volumeSize = 3.f; // meters
        static constexpr float voxelSize = volumeSize / 128.f; //meters
        static constexpr float tsdfTruncateDistance = 2 * voxelSize;
        static constexpr float maxDepth = 0.f;
        static const int maxWeight = 64; // number of frames
        static constexpr float raycastStepFactor = 0.75f;
        static const bool zFirstMemOrder = true; // order of voxels in volume

        const Matx33f  cameraIntegrateIntrinsics = Matx33f(fx, 0, cx, 0, fy, cy, 0, 0, 1); // camera settings
        const Matx33f  rgb_cameraIntegrateIntrinsics = Matx33f(rgb_fx, 0, rgb_cx, 0, rgb_fy, rgb_cy, 0, 0, 1); // camera settings
        const Matx33f  cameraRaycastIntrinsics = Matx33f(fx, 0, cx, 0, fy, cy, 0, 0, 1); // camera settings
        const Matx33f  rgb_cameraRaycastIntrinsics = Matx33f(rgb_fx, 0, rgb_cx, 0, rgb_fy, rgb_cy, 0, 0, 1); // camera settings
        const Affine3f volumePose = Affine3f().translate(Vec3f(-volumeSize / 2.f, -volumeSize / 2.f, 0.5f));
        const Matx44f  volumePoseMatrix = volumePose.matrix;
        // Unlike original code, this should work with any volume size
        // Not only when (x,y,z % 32) == 0
        const Point3i  volumeResolution = Vec3i::all(128); //number of voxels
    };

};


VolumeSettings::VolumeSettings()
{
    this->impl = makePtr<VolumeSettingsImpl>();
}

VolumeSettings::VolumeSettings(VolumeType volumeType)
{
    this->impl = makePtr<VolumeSettingsImpl>(volumeType);
}

VolumeSettings::~VolumeSettings() {}

void  VolumeSettings::setWidth(int val) { this->impl->setWidth(val); };
int   VolumeSettings::getWidth() const { return this->impl->getWidth(); };
void  VolumeSettings::setHeight(int val) { this->impl->setHeight(val); };
int   VolumeSettings::getHeight() const { return this->impl->getHeight(); };
void  VolumeSettings::setVoxelSize(float val) { this->impl->setVoxelSize(val); };
float VolumeSettings::getVoxelSize() const { return this->impl->getVoxelSize(); };
void  VolumeSettings::setRaycastStepFactor(float val) { this->impl->setRaycastStepFactor(val); };
float VolumeSettings::getRaycastStepFactor() const { return this->impl->getRaycastStepFactor(); };
void  VolumeSettings::setTsdfTruncateDistance(float val) { this->impl->setTsdfTruncateDistance(val); };
float VolumeSettings::getTsdfTruncateDistance() const { return this->impl->getTsdfTruncateDistance(); };
void  VolumeSettings::setMaxDepth(float val) { this->impl->setMaxDepth(val); };
float VolumeSettings::getMaxDepth() const { return this->impl->getMaxDepth(); };
void  VolumeSettings::setDepthFactor(float val) { this->impl->setDepthFactor(val); };
float VolumeSettings::getDepthFactor() const { return this->impl->getDepthFactor(); };
void  VolumeSettings::setMaxWeight(int val) { this->impl->setMaxWeight(val); };
int   VolumeSettings::getMaxWeight() const { return this->impl->getMaxWeight(); };

void VolumeSettings::setVolumePose(InputArray val) { this->impl->setVolumePose(val); };
void VolumeSettings::getVolumePose(OutputArray val) const { this->impl->getVolumePose(val); };
void VolumeSettings::setVolumeResolution(InputArray val) { this->impl->setVolumeResolution(val); };
void VolumeSettings::getVolumeResolution(OutputArray val) const { this->impl->getVolumeResolution(val); };
void VolumeSettings::getVolumeDimentions(OutputArray val) const { this->impl->getVolumeDimentions(val); };
void VolumeSettings::setCameraIntegrateIntrinsics(InputArray val) { this->impl->setCameraIntegrateIntrinsics(val); };
void VolumeSettings::getCameraIntegrateIntrinsics(OutputArray val) const { this->impl->getCameraIntegrateIntrinsics(val); };
void VolumeSettings::setRGBCameraIntegrateIntrinsics(InputArray val) { this->impl->setRGBCameraIntegrateIntrinsics(val); };
void VolumeSettings::getRGBCameraIntegrateIntrinsics(OutputArray val) const { this->impl->getRGBCameraIntegrateIntrinsics(val); };
void VolumeSettings::setCameraRaycastIntrinsics(InputArray val) { this->impl->setCameraRaycastIntrinsics(val); };
void VolumeSettings::getCameraRaycastIntrinsics(OutputArray val) const { this->impl->getCameraRaycastIntrinsics(val); };
void VolumeSettings::setRGBCameraRaycastIntrinsics(InputArray val) { this->impl->setRGBCameraRaycastIntrinsics(val); };
void VolumeSettings::getRGBCameraRaycastIntrinsics(OutputArray val) const { this->impl->getRGBCameraRaycastIntrinsics(val); };


VolumeSettingsImpl::VolumeSettingsImpl()
    : VolumeSettingsImpl(VolumeType::TSDF)
{
}

VolumeSettingsImpl::VolumeSettingsImpl(VolumeType _volumeType)
{
    volumeType = _volumeType;
    if (volumeType == VolumeType::TSDF)
    {
        DefaultTsdfSets ds = DefaultTsdfSets();

        this->width = ds.width;
        this->height = ds.height;
        this->depthFactor = ds.depthFactor;
        this->voxelSize = ds.voxelSize;
        this->tsdfTruncateDistance = ds.tsdfTruncateDistance;
        this->maxDepth = ds.maxDepth;
        this->maxWeight = ds.maxWeight;
        this->raycastStepFactor = ds.raycastStepFactor;
        this->zFirstMemOrder = ds.zFirstMemOrder;

        this->volumePose = ds.volumePoseMatrix;
        this->volumeResolution = ds.volumeResolution;
        this->volumeDimentions = calcVolumeDimentions(ds.volumeResolution, ds.zFirstMemOrder);
        this->cameraIntegrateIntrinsics = ds.cameraIntegrateIntrinsics;
        this->rgb_cameraIntegrateIntrinsics = ds.rgb_cameraIntegrateIntrinsics;
        this->cameraRaycastIntrinsics = ds.cameraRaycastIntrinsics;
        this->rgb_cameraRaycastIntrinsics = ds.rgb_cameraRaycastIntrinsics;
    }
    else if (volumeType == VolumeType::HashTSDF)
    {
        DefaultHashTsdfSets ds = DefaultHashTsdfSets();

        this->width = ds.width;
        this->height = ds.height;
        this->depthFactor = ds.depthFactor;
        this->voxelSize = ds.voxelSize;
        this->tsdfTruncateDistance = ds.tsdfTruncateDistance;
        this->maxDepth = ds.maxDepth;
        this->maxWeight = ds.maxWeight;
        this->raycastStepFactor = ds.raycastStepFactor;
        this->zFirstMemOrder = ds.zFirstMemOrder;

        this->volumePose = ds.volumePoseMatrix;
        this->volumeResolution = ds.volumeResolution;
        this->volumeDimentions = calcVolumeDimentions(ds.volumeResolution, ds.zFirstMemOrder);
        this->cameraIntegrateIntrinsics = ds.cameraIntegrateIntrinsics;
        this->rgb_cameraIntegrateIntrinsics = ds.rgb_cameraIntegrateIntrinsics;
        this->cameraRaycastIntrinsics = ds.cameraRaycastIntrinsics;
        this->rgb_cameraRaycastIntrinsics = ds.rgb_cameraRaycastIntrinsics;
    }
    else if (volumeType == VolumeType::ColorTSDF)
    {
        DefaultColorTsdfSets ds = DefaultColorTsdfSets();

        this->width = ds.width;
        this->height = ds.height;
        this->depthFactor = ds.depthFactor;
        this->voxelSize = ds.voxelSize;
        this->tsdfTruncateDistance = ds.tsdfTruncateDistance;
        this->maxDepth = ds.maxDepth;
        this->maxWeight = ds.maxWeight;
        this->raycastStepFactor = ds.raycastStepFactor;
        this->zFirstMemOrder = ds.zFirstMemOrder;

        this->volumePose = ds.volumePoseMatrix;
        this->volumeResolution = ds.volumeResolution;
        this->volumeDimentions = calcVolumeDimentions(ds.volumeResolution, ds.zFirstMemOrder);
        this->cameraIntegrateIntrinsics = ds.cameraIntegrateIntrinsics;
        this->rgb_cameraIntegrateIntrinsics = ds.rgb_cameraIntegrateIntrinsics;
        this->cameraRaycastIntrinsics = ds.cameraRaycastIntrinsics;
        this->rgb_cameraRaycastIntrinsics = ds.rgb_cameraRaycastIntrinsics;
    }
}


VolumeSettingsImpl::~VolumeSettingsImpl() {}


void VolumeSettingsImpl::setWidth(int val)
{
    this->width = val;
}

int VolumeSettingsImpl::getWidth() const
{
    return this->width;
}

void VolumeSettingsImpl::setHeight(int val)
{
    this->height = val;
}

int VolumeSettingsImpl::getHeight() const
{
    return this->height;
}

void VolumeSettingsImpl::setDepthFactor(float val)
{
    this->depthFactor = val;
}

float VolumeSettingsImpl::getDepthFactor() const
{
    return this->depthFactor;
}

void VolumeSettingsImpl::setVoxelSize(float  val)
{
    this->voxelSize = val;
}

float VolumeSettingsImpl::getVoxelSize() const
{
    return this->voxelSize;
}

void VolumeSettingsImpl::setTsdfTruncateDistance(float val)
{
    this->tsdfTruncateDistance = val;
}

float VolumeSettingsImpl::getTsdfTruncateDistance() const
{
    return this->tsdfTruncateDistance;
}

void VolumeSettingsImpl::setMaxDepth(float val)
{
    this->maxDepth = val;
}

float VolumeSettingsImpl::getMaxDepth() const
{
    return this->maxDepth;
}

void VolumeSettingsImpl::setMaxWeight(int val)
{
    this->maxWeight = val;
}

int VolumeSettingsImpl::getMaxWeight() const
{
    return this->maxWeight;
}

void VolumeSettingsImpl::setRaycastStepFactor(float val)
{
    this->raycastStepFactor = val;
}

float VolumeSettingsImpl::getRaycastStepFactor() const
{
    return this->raycastStepFactor;
}

void VolumeSettingsImpl::setVolumePose(InputArray val)
{
    if (!val.empty())
    {
        val.copyTo(this->volumePose);
    }
}

void VolumeSettingsImpl::getVolumePose(OutputArray val) const
{
    Mat(this->volumePose).copyTo(val);
}

void VolumeSettingsImpl::setVolumeResolution(InputArray val)
{
    if (!val.empty())
    {
        this->volumeResolution = Point3i(val.getMat());
    }
}

void VolumeSettingsImpl::getVolumeResolution(OutputArray val) const
{
    Mat(this->volumeResolution).copyTo(val);
}

void VolumeSettingsImpl::getVolumeDimentions(OutputArray val) const
{
    Mat(this->volumeDimentions).copyTo(val);
}

void VolumeSettingsImpl::setCameraIntegrateIntrinsics(InputArray val)
{
    if (!val.empty())
    {
        this->cameraIntegrateIntrinsics = Matx33f(val.getMat());
    }
}

void VolumeSettingsImpl::getCameraIntegrateIntrinsics(OutputArray val) const
{
    Mat(this->cameraIntegrateIntrinsics).copyTo(val);
}

void VolumeSettingsImpl::setRGBCameraIntegrateIntrinsics(InputArray val)
{
    if (!val.empty())
    {
        this->rgb_cameraIntegrateIntrinsics = Matx33f(val.getMat());
    }
}

void VolumeSettingsImpl::getRGBCameraIntegrateIntrinsics(OutputArray val) const
{
    Mat(this->rgb_cameraIntegrateIntrinsics).copyTo(val);
}


void VolumeSettingsImpl::setCameraRaycastIntrinsics(InputArray val)
{
    if (!val.empty())
    {
        this->cameraRaycastIntrinsics = Matx33f(val.getMat());
    }
}

void VolumeSettingsImpl::getCameraRaycastIntrinsics(OutputArray val) const
{
    Mat(this->cameraRaycastIntrinsics).copyTo(val);
}

void VolumeSettingsImpl::setRGBCameraRaycastIntrinsics(InputArray val)
{
    if (!val.empty())
    {
        this->rgb_cameraRaycastIntrinsics = Matx33f(val.getMat());
    }
}

void VolumeSettingsImpl::getRGBCameraRaycastIntrinsics(OutputArray val) const
{
    Mat(this->rgb_cameraRaycastIntrinsics).copyTo(val);
}


Vec4i calcVolumeDimentions(Point3i volumeResolution, bool ZFirstMemOrder)
{
    // (xRes*yRes*zRes) array
    // Depending on zFirstMemOrder arg:
    // &elem(x, y, z) = data + x*zRes*yRes + y*zRes + z;
    // &elem(x, y, z) = data + x + y*xRes + z*xRes*yRes;
    int xdim, ydim, zdim;
    if (ZFirstMemOrder)
    {
        xdim = volumeResolution.z * volumeResolution.y;
        ydim = volumeResolution.z;
        zdim = 1;
    }
    else
    {
        xdim = 1;
        ydim = volumeResolution.x;
        zdim = volumeResolution.x * volumeResolution.y;
    }
    return Vec4i(xdim, ydim, zdim);
}

}
