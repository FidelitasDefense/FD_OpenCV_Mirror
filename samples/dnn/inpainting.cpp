/*
This sample inpaints the masked area in the given image.

Copyright (C) 20245, Bigvision LLC.

How to use:
    sample command to run:

        ./example_dnn_inpainting
    The system will ask you to draw the mask to be inpainted

    You can download lama inpainting model using:
       `python download_models.py lama`

    Set environment variable OPENCV_DOWNLOAD_CACHE_DIR to point to the directory where models are downloaded. Also, point OPENCV_SAMPLES_DATA_PATH to opencv/samples/data.
*/

#include <iostream>
#include <fstream>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/dnn.hpp>

#include "common.hpp"

using namespace cv;
using namespace dnn;
using namespace std;

const string about = "Use this script for Person Re-identification using OpenCV. \n\n"
        "Firstly, download required models i.e. lama using `download_models.py` (if not already done). Set environment variable OPENCV_DOWNLOAD_CACHE_DIR to point to the directory where models are downloaded. Also, point OPENCV_SAMPLES_DATA_PATH to opencv/samples/data.\n"
        "To run:\n"
        "\t Example: ./example_dnn_inpainting\n\n"
        "Inpainting model path can also be specified using --model argument.\n\n";

const string param_keys =
    "{help    h  |                   | show help message}"
    "{ @alias    |       lama        | An alias name of model to extract preprocessing parameters from models.yml file. }"
    "{ zoo       | ../dnn/models.yml | An optional path to file with preprocessing parameters }"
    "{input   i  |                   | image file path}";

const string backend_keys = format(
    "{ backend | default | Choose one of computation backends: "
    "default: automatically (by default), "
    "openvino: Intel's Deep Learning Inference Engine (https://software.intel.com/openvino-toolkit), "
    "opencv: OpenCV implementation, "
    "vkcom: VKCOM, "
    "cuda: CUDA, "
    "webnn: WebNN }");

const string target_keys = format(
    "{ target | cpu | Choose one of target computation devices: "
    "cpu: CPU target (by default), "
    "opencl: OpenCL, "
    "opencl_fp16: OpenCL fp16 (half-float precision), "
    "vpu: VPU, "
    "vulkan: Vulkan, "
    "cuda: CUDA, "
    "cuda_fp16: CUDA fp16 (half-float preprocess) }");

string keys = param_keys + backend_keys + target_keys;
bool drawing = false;
Mat maskGray;
int brush_size = 25;


static void drawMask(int event, int x, int y, int, void*) {
    if (event == EVENT_LBUTTONDOWN) {
        drawing = true;
    } else if (event == EVENT_MOUSEMOVE) {
        if (drawing) {
            circle(maskGray, Point(x, y), brush_size, Scalar(255), -1);
        }
    } else if (event == EVENT_LBUTTONUP) {
        drawing = false;
    }
}

int main(int argc, char **argv)
{
    CommandLineParser parser(argc, argv, keys);

    if (!parser.has("@alias") || parser.has("help"))
    {
        cout<<about<<endl;
        parser.printMessage();
        return 0;
    }
    string modelName = parser.get<String>("@alias");
    string zooFile = findFile(parser.get<String>("zoo"));
    keys += genPreprocArguments(modelName, zooFile);
    parser = CommandLineParser(argc, argv, keys);
    parser.about("Use this script to run ReID networks using OpenCV.");

    const string sha1 = parser.get<String>("sha1");
    const string modelPath = findModel(parser.get<String>("model"), sha1);
    string imgPath = parser.get<String>("input");
    const string backend = parser.get<String>("backend");
    const string target = parser.get<String>("target");
    int height = parser.get<int>("height");
    int width = parser.get<int>("width");
    float scale = parser.get<float>("scale");
    bool swapRB = parser.get<bool>("rgb");
    Scalar mean_v = parser.get<Scalar>("mean");
    int stdSize = 20;
    int stdWeight = 400;
    int stdImgSize = 512;
    int imgWidth = -1; // Initialization
    int fontSize = 50;
    int fontWeight = 500;

    EngineType engine = ENGINE_AUTO;
    if (backend != "default" || target != "cpu"){
        engine = ENGINE_CLASSIC;
    }

    Net net = readNetFromONNX(modelPath, engine);
    net.setPreferableBackend(getBackendID(backend));
    net.setPreferableTarget(getTargetID(target));

    FontFace fontFace("sans");

    VideoCapture cap;
    Mat image = imread(findFile(imgPath));
    if (image.empty()) {
        cerr << "Error: Input image could not be loaded." << endl;
        return -1;
    }
    double aspectRatio = static_cast<double>(image.rows) / static_cast<double>(image.cols);
    imgWidth = min(image.rows, image.cols);
    fontSize = min(fontSize, (stdSize*imgWidth)/stdImgSize);
    fontWeight = min(fontWeight, (stdWeight*imgWidth)/stdImgSize);

    maskGray = Mat::zeros(image.size(), CV_8U);

    namedWindow("Draw Mask");
    setMouseCallback("Draw Mask", drawMask);

    const string label = "Draw the mask on the image. Press space bar when done ";

    for(;;) {
        Mat displayImage = image.clone();
        Mat overlay = image.clone();

        double alpha = 0.5;
        Rect r = getTextSize(Size(), label, Point(), fontFace, fontSize, fontWeight);
        r.height += 2 * fontSize; // padding
        r.width += 10; // padding
        rectangle(overlay, r, Scalar::all(255), FILLED);
        addWeighted(overlay, alpha, displayImage, 1 - alpha, 0, displayImage);
        putText(displayImage, label, Point(10, fontSize), Scalar(0,0,0), fontFace, fontSize, fontWeight);
        putText(displayImage, "Press 'i' to increase and 'd' to decrease brush size", Point(10, 2*fontSize), Scalar(0,0,0), fontFace, fontSize, fontWeight);

        displayImage.setTo(Scalar(255, 255, 255), maskGray > 0); // Highlight mask area
        imshow("Draw Mask", displayImage);

        char key = waitKey(1);
        if (key == 'i') {
            brush_size += 1;
            cout << "Brush size increased to " << brush_size << endl;
        } else if (key == 'd') {
            brush_size = max(1, brush_size - 1);
            cout << "Brush size decreased to " << brush_size << endl;
        } else if (key == ' ') {
            break;
        } else if (key == 27){
            return -1;
        }
    }
    destroyAllWindows();

    Mat image_blob = blobFromImage(image, scale, Size(width, height), mean_v, swapRB, false);

    Mat mask_blob;
    mask_blob = blobFromImage(maskGray, 1.0, Size(width, height), Scalar(0), false, false);
    mask_blob = (mask_blob > 0);
    mask_blob.convertTo(mask_blob, CV_32F);
    mask_blob = mask_blob/255.0;

    net.setInput(image_blob, "image");
    net.setInput(mask_blob, "mask");

    Mat output = net.forward();

    Mat output_transposed(3, &output.size[1], CV_32F, output.ptr<float>());

    vector<Mat> channels;
    for (int i = 0; i < 3; ++i) {
        channels.push_back(Mat(output_transposed.size[1], output_transposed.size[2], CV_32F,
                                    output_transposed.ptr<float>(i)));
    }
    Mat output_image;
    merge(channels, output_image);
    output_image.convertTo(output_image, CV_8U);

    int h = static_cast<int>(width * aspectRatio);
    resize(output_image, output_image, Size(width, h));

    imshow("Inpainted Output", output_image);
    waitKey(0);

    return 0;
}
