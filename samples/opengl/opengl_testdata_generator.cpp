#include <iostream>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN 1
#define NOMINMAX 1
#include <windows.h>
#endif

#if defined(__APPLE__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "opencv2/core.hpp"
#include "opencv2/core/opengl.hpp"
#include "opencv2/3d.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"

using namespace std;
using namespace cv;
using namespace cv::cuda;

// model data should be identical to the code from tests
enum class ModelType
{
    Empty = 0,
    File = 1,
    Clipping = 2,
    Color = 3,
    Centered = 4
};

static void generateNormals(const std::vector<Vec3f>& points, const std::vector<std::vector<int>>& indices,
                            std::vector<Vec3f>& normals)
{
    std::vector<std::vector<Vec3f>> preNormals(points.size(), std::vector<Vec3f>());

    for (const auto& tri : indices)
    {
        Vec3f p0 = points[tri[0]];
        Vec3f p1 = points[tri[1]];
        Vec3f p2 = points[tri[2]];

        Vec3f cross = cv::normalize((p1 - p0).cross(p2 - p0));
        for (int i = 0; i < 3; i++)
        {
            preNormals[tri[i]].push_back(cross);
        }
    }

    normals.reserve(points.size());
    for (const auto& pn : preNormals)
    {
        Vec3f sum { };
        for (const auto& n : pn)
        {
            sum += n;
        }
        normals.push_back(cv::normalize(sum));
    }
}

class ModelData
{
public:
    ModelData(ModelType type = ModelType::Empty, std::string objPath = { })
    {
        switch (type)
        {
        case ModelType::Empty:
        {
            position = Vec3d(0.0, 0.0, 0.0);
            lookat   = Vec3d(0.0, 0.0, 0.0);
            upVector = Vec3d(0.0, 1.0, 0.0);

            fovy = 45.0;

            vertices = std::vector<Vec3f>(4, {2.0f, 0, -2.0f});
            colors   = std::vector<Vec3f>(4, {0, 0, 1.0f});
            indices = { };
        }
        break;
        case ModelType::File:
        {
            position = Vec3d( 1.9, 0.4, 1.3);
            lookat   = Vec3d( 0.0, 0.0, 0.0);
            upVector = Vec3d( 0.0, 1.0, 0.0);

            fovy = 45.0;

            objectPath = objPath;
            std::vector<vector<int>> indvec;
            loadMesh(objectPath, vertices, indvec);
            // using per-vertex normals as colors
            generateNormals(vertices, indvec, colors);
            if (vertices.size() != colors.size())
            {
                std::runtime_error("Model should contain normals for each vertex");
            }
            for (const auto &vec : indvec)
            {
                indices.push_back({vec[0], vec[1], vec[2]});
            }

            for (auto &color : colors)
            {
                color = Vec3f(abs(color[0]), abs(color[1]), abs(color[2]));
            }
        }
        break;
        case ModelType::Clipping:
        {
            position = Vec3d(0.0, 0.0, 5.0);
            lookat   = Vec3d(0.0, 0.0, 0.0);
            upVector = Vec3d(0.0, 1.0, 0.0);

            fovy = 45.0;

            vertices =
            {
                { 2.0,  0.0, -2.0}, { 0.0, -6.0, -2.0}, {-2.0,  0.0, -2.0},
                { 3.5, -1.0, -5.0}, { 2.5, -2.5, -5.0}, {-1.0,  1.0, -5.0},
                {-6.5, -1.0, -3.0}, {-2.5, -2.0, -3.0}, { 1.0,  1.0, -5.0},
            };

            indices = { {0, 1, 2}, {3, 4, 5}, {6, 7, 8} };

            Vec3f col1(217.0, 238.0, 185.0);
            Vec3f col2(185.0, 217.0, 238.0);
            Vec3f col3(150.0,  10.0, 238.0);

            col1 *= (1.f / 255.f);
            col2 *= (1.f / 255.f);
            col3 *= (1.f / 255.f);

            colors =
            {
                col1, col2, col3,
                col2, col3, col1,
                col3, col1, col2,
            };
        }
        break;
        case ModelType::Centered:
        {
            position = Vec3d(0.0, 0.0, 5.0);
            lookat   = Vec3d(0.0, 0.0, 0.0);
            upVector = Vec3d(0.0, 1.0, 0.0);

            fovy = 45.0;

            vertices =
            {
                { 2.0,  0.0, -2.0}, { 0.0, -2.0, -2.0}, {-2.0,  0.0, -2.0},
                { 3.5, -1.0, -5.0}, { 2.5, -1.5, -5.0}, {-1.0,  0.5, -5.0},
            };

            indices = { {0, 1, 2}, {3, 4, 5} };

            Vec3f col1(217.0, 238.0, 185.0);
            Vec3f col2(185.0, 217.0, 238.0);

            col1 *= (1.f / 255.f);
            col2 *= (1.f / 255.f);

            colors =
            {
                col1, col2, col1,
                col2, col1, col2,
            };
        }
        break;
        case ModelType::Color:
        {
            position = Vec3d(0.0, 0.0, 5.0);
            lookat   = Vec3d(0.0, 0.0, 0.0);
            upVector = Vec3d(0.0, 1.0, 0.0);

            fovy = 60.0;

            vertices =
            {
                { 2.0,  0.0, -2.0},
                { 0.0,  2.0, -3.0},
                {-2.0,  0.0, -2.0},
                { 0.0, -2.0,  1.0},
            };

            indices = { {0, 1, 2}, {0, 2, 3} };

            colors =
            {
                { 0.0f, 0.0f, 1.0f},
                { 0.0f, 1.0f, 0.0f},
                { 1.0f, 0.0f, 0.0f},
                { 0.0f, 1.0f, 0.0f},
            };
        }
        break;

        default:
            CV_Error(Error::StsBadArg, "Unknown model type");
            break;
        }
    }

    Vec3d position;
    Vec3d lookat;
    Vec3d upVector;

    double fovy;

    std::vector<Vec3f> vertices;
    std::vector<Vec3i> indices;
    std::vector<Vec3f> colors;

    string objectPath;
};


struct DrawData
{
    ogl::Arrays arr;
    ogl::Buffer indices;
};

void draw(void* userdata);

void draw(void* userdata)
{
    DrawData* data = static_cast<DrawData*>(userdata);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    ogl::render(data->arr, data->indices, ogl::TRIANGLES);
}

static void generateImage(cv::Size imgSz, TriangleShadingType shadingType, TriangleCullingMode cullingMode,
                          ModelType modelType, std::string modelPath, cv::Mat& colorImage, cv::Mat& depthImage)
{
    namedWindow("OpenGL", WINDOW_OPENGL);
    resizeWindow("OpenGL", imgSz.width, imgSz.height);

    ModelData modelData(modelType, modelPath);

    DrawData data;

    std::vector<Vec3f> vertices;
    std::vector<Vec4f> colors4f;
    std::vector<int> idxLinear;

    if (shadingType == RASTERIZE_SHADING_FLAT)
    {
        // rearrange vertices and colors for flat shading
        int ctr = 0;
        for (const auto& idx : modelData.indices)
        {
            for (int i = 0; i < 3; i++)
            {
                vertices.push_back(modelData.vertices[idx[i]]);
                idxLinear.push_back(ctr++);
            }

            Vec3f ci = modelData.colors[idx[0]];
            for (int i = 0; i < 3; i++)
            {
                colors4f.emplace_back(ci[0], ci[1], ci[2], 1.f);
            }
        }
    }
    else
    {
        vertices = modelData.vertices;
        for (const auto& c : modelData.colors)
        {
            Vec3f ci = (shadingType == RASTERIZE_SHADING_SHADED) ? c: cv::Vec3f::all(1.f);
            colors4f.emplace_back(ci[0], ci[1], ci[2], 1.0);
        }

        for (const auto& idx : modelData.indices)
        {
            for (int i = 0; i < 3; i++)
            {
                idxLinear.push_back(idx[i]);
            }
        }
    }

    data.arr.setVertexArray(vertices);
    data.arr.setColorArray(colors4f);
    data.indices.copyFrom(idxLinear);

    double zNear = 0.1, zFar = 50;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(modelData.fovy, (double)imgSz.width / imgSz.height, zNear, zFar);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    //gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
    gluLookAt(modelData.position[0], modelData.position[1], modelData.position[2],
              modelData.lookat  [0], modelData.lookat  [1], modelData.lookat  [2],
              modelData.upVector[0], modelData.upVector[1], modelData.upVector[2]);

    if (cullingMode == RASTERIZE_CULLING_NONE)
    {
        glDisable(GL_CULL_FACE);
    }
    else
    {
        glEnable(GL_CULL_FACE);
        glCullFace(GL_FRONT);
        if (cullingMode == RASTERIZE_CULLING_CW)
        {
            glFrontFace(GL_CW);
        }
        else
        {
            glFrontFace(GL_CCW);
        }
    }

    glEnable(GL_DEPTH_TEST);

    cv::setOpenGlDrawCallback("OpenGL", draw, &data);

    const int framesToSkip = 10;
    for (int f = 0; f < framesToSkip; f++)
    {
        updateWindow("OpenGL");

        colorImage = cv::Mat(imgSz.height, imgSz.width, CV_8UC3);
        glReadPixels(0, 0, imgSz.width, imgSz.height, GL_RGB, GL_UNSIGNED_BYTE, colorImage.data);
        cv::cvtColor(colorImage, colorImage, cv::COLOR_RGB2BGR);
        cv::flip(colorImage, colorImage, 0);

        depthImage = cv::Mat(imgSz.height, imgSz.width, CV_32F);
        glReadPixels(0, 0, imgSz.width, imgSz.height, GL_DEPTH_COMPONENT, GL_FLOAT, depthImage.data);
        // map from [0, 1] to [zNear, zFar]
        for (auto it = depthImage.begin<float>(); it != depthImage.end<float>(); ++it)
        {
            *it = (float)(zNear * zFar / (double(*it) * (zNear - zFar) + zFar));
        }
        cv::flip(depthImage, depthImage, 0);
        depthImage.convertTo(depthImage, CV_16U, 1000.0);

        char key = (char)waitKey(40);
        if (key == 27)
            break;
    }

    cv::setOpenGlDrawCallback("OpenGL", 0, 0);
    cv::destroyAllWindows();
}


int main(int argc, char* argv[])
{
    cv::CommandLineParser parser(argc, argv,
            "{ help h usage ? |      | show this message }"
            "{ outPath        |      | output path for generated images }"
            "{ modelPath      |      | path to 3d model to render }"
    );
    parser.about("This app is used to generate test data for triangleRasterize() function");

    if (parser.has("help"))
    {
        parser.printMessage();
        return 0;
    }

    std::string modelPath = parser.get<std::string>("modelPath");
    if (modelPath.empty())
    {
        std::cout << "No model path given" << std::endl;
        return -1;
    }

    std::string outPath = parser.get<std::string>("outPath");
    if (outPath.empty())
    {
        std::cout << "No output path given" << std::endl;
        return -1;
    }

    std::array<cv::Size, 4> resolutions = { cv::Size {700, 700}, cv::Size {640, 480}, cv::Size(256, 256), cv::Size(320, 240) };
    for (const auto& res : resolutions)
    {
        for (const auto shadingType : {
                RASTERIZE_SHADING_WHITE,
                RASTERIZE_SHADING_FLAT,
                RASTERIZE_SHADING_SHADED
            })
        {
            std::string shadingName;
            switch (shadingType)
            {
            case RASTERIZE_SHADING_WHITE:  shadingName = "White";  break;
            case RASTERIZE_SHADING_FLAT:   shadingName = "Flat";   break;
            case RASTERIZE_SHADING_SHADED: shadingName = "Shaded"; break;
            default:
                break;
            }

            for (const auto cullingMode : {
                    RASTERIZE_CULLING_NONE,
                    RASTERIZE_CULLING_CW,
                    RASTERIZE_CULLING_CCW
            })
            {
                std::string cullingName;
                switch (cullingMode)
                {
                    case RASTERIZE_CULLING_NONE: cullingName = "None"; break;
                    case RASTERIZE_CULLING_CW:   cullingName = "CW"; break;
                    case RASTERIZE_CULLING_CCW:  cullingName = "CCW"; break;
                    default: break;
                }

                for (const auto modelType : {
                            ModelType::File,
                            ModelType::Clipping,
                            ModelType::Color,
                            ModelType::Centered,
                    })
                {
                    std::string modelName;
                    switch (modelType)
                    {
                    case ModelType::File:     modelName = "File";     break;
                    case ModelType::Clipping: modelName = "Clipping"; break;
                    case ModelType::Color:    modelName = "Color";    break;
                    case ModelType::Centered: modelName = "Centered"; break;
                    default:
                        break;
                    }

                    std::string suffix = cv::format("%s_%dx%d_Cull%s", modelName.c_str(), res.width, res.height, cullingName.c_str());

                    std::cout << suffix + "_" + shadingName << "..." << std::endl;

                    cv::Mat colorImage, depthImage;
                    generateImage(res, shadingType, cullingMode, modelType, modelPath, colorImage, depthImage);

                    std::string gtPathColor = outPath + "/example_image_" + suffix + "_" + shadingName + ".png";
                    std::string gtPathDepth = outPath + "/depth_image_"   + suffix + ".png";

                    cv::imwrite(gtPathColor, colorImage);
                    cv::imwrite(gtPathDepth, depthImage);
                }
            }
        }
    }

    return 0;
}
