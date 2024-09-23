/*
This sample detects the query person in the given video file.

Authors of samples and Youtu ReID baseline:
        Xing Sun <winfredsun@tencent.com>
        Feng Zheng <zhengf@sustech.edu.cn>
        Xinyang Jiang <sevjiang@tencent.com>
        Fufu Yu <fufuyu@tencent.com>
        Enwei Zhang <miyozhang@tencent.com>

Copyright (C) 2020-2021, Tencent.
Copyright (C) 2020-2021, SUSTech.

How to use:
    sample command to run:

        ./person_reid --video=/path/to/videofile --model=path/to/youtu_reid_baseline_medium.onnx --yolo=path/to/yolov8n.onnx
    The system will ask you to mark the person to be tracked

    You can download a baseline ReID model from:
        https://github.com/ReID-Team/ReID_extra_testdata
    Drive Link: https://drive.google.com/drive/folders/1wFGcuolSzX3_PqNKb4BAV3DNac7tYpc2
*/

#include <iostream>
#include <fstream>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/dnn.hpp>

#include "common.hpp"

using namespace cv;
using namespace cv::dnn;
using namespace std;

const string about = "Use this script for Person Re-identification using OpenCV. \n\n"
        "Firstly, download required models i.e. reid and yolov8 using `download_models.py` (if not already done). Set environment variable OPENCV_DOWNLOAD_CACHE_DIR to point to the directory where models are downloaded. Also, point OPENCV_SAMPLES_DATA_PATH to opencv/samples/data.\n"
        "To run:\n"
        "\t Example: ./example_dnn_person_reid reid\n\n"
        "Re-Identification model path can also be specified using --model argument. Detection model can be set using --yolo_model argument.\n\n";

const string param_keys =
    "{help    h  |                   | show help message}"
    "{ @alias    |                   | An alias name of model to extract preprocessing parameters from models.yml file. }"
    "{ zoo       | ../dnn/models.yml | An optional path to file with preprocessing parameters }"
    "{query   q  |                   | Path to target image. Skip this argument to select target in the video frame.}"
    "{input   i  |                   | video file path}"
    "{yolo_model |                   | Path to yolov8n.onnx}";

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


struct MatComparator
{
    bool operator()(const Mat &a, const Mat &b) const
    {
        return a.data < b.data; // This is a simple pointer comparison, not content!
    }
};

// Global variable for drawing function to select target
map<Mat, Rect, MatComparator> imgDict;
bool drawing = false;
int ix = -1, iy = -1;
Rect rect;

string sha1, yoloSha1, modelPath;
string queryImagePath, videoPath, yoloPath, backend, target;
int height, width, yoloHeight, yoloWidth;
float scale, yoloScale;
bool swapRB, yoloSwapRB;
Scalar mean_v, stnd;

static void extractFrames(Net& reidNet);

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
    keys += genPreprocArguments(modelName, zooFile, "yolo_");
    parser = CommandLineParser(argc, argv, keys);
    parser.about("Use this script to run ReID networks using OpenCV.");

    sha1 = parser.get<String>("sha1");
    yoloSha1 = parser.get<String>("yolo_sha1");
    modelPath = findModel(parser.get<String>("model"), sha1);
    queryImagePath = parser.get<String>("query");
    videoPath = parser.get<String>("input");
    yoloPath = findModel(parser.get<String>("yolo_model"), yoloSha1);
    backend = parser.get<String>("backend");
    target = parser.get<String>("target");
    height = parser.get<int>("height");
    width = parser.get<int>("width");
    yoloHeight = parser.get<int>("yolo_height");
    yoloWidth = parser.get<int>("yolo_width");
    scale = parser.get<float>("scale");
    yoloScale = parser.get<float>("yolo_scale");
    swapRB = parser.get<bool>("rgb");
    yoloSwapRB = parser.get<bool>("yolo_rgb");
    mean_v = parser.get<Scalar>("mean");
    stnd = parser.get<Scalar>("std");

    Net net = readNetFromONNX(modelPath);
    net.setPreferableBackend(getBackendID(backend));
    net.setPreferableTarget(getTargetID(target));

    if(yoloPath.empty()){
        cout<<"[ERROR] Please pass path to yolov8.onnx model file using --yolo."<<endl;
        return -1;
    }

    extractFrames(net);
    return 0;
}

static void extractFeatures(vector<Mat> &imglist, Net &net, vector<Mat> &features)
{
    for (int st = 0; st < (int)imglist.size(); st++)
    {
        Mat blob;
        blobFromImage(imglist[st], blob, scale, Size(width, height), mean_v, swapRB, false, CV_32F);

        // Check if standard deviation values are non-zero
        if (stnd[0] != 0.0 && stnd[1] != 0.0 && stnd[2] != 0.0)
        {
            // Divide blob by std for each channel
            divide(blob, stnd, blob);
        }
        net.setInput(blob);
        Mat out=net.forward();
        vector<int> s {out.size[0], out.size[1]};
        out = out.reshape(1, s);
        for (int i = 0; i < out.rows; i++)
        {
            Mat temp_feature(1, out.cols, CV_32F);
            for (int j = 0; j < out.cols; j++)
            {
                temp_feature.at<float>(0, j) = out.at<float>(i, j);
            }
            Mat norm_feature;
            normalize(temp_feature, norm_feature, 1.0, 0.0, NORM_L2);
            features.push_back(norm_feature);
        }
    }
    return;
}

static int findMatching(const Mat &queryFeatures, const vector<Mat> &galleryFeatures)
{
    if (queryFeatures.empty() || galleryFeatures.empty())
        return -1; // No valid index if either feature list is empty

    int bestIndex = -1;
    float maxSimilarity = FLT_MIN;

    for (int j = 0; j < (int)galleryFeatures.size(); j++)
    {
        float currentSimilarity = static_cast<float>(queryFeatures.dot(galleryFeatures[j]));
        if (currentSimilarity > maxSimilarity)
        {
            maxSimilarity = currentSimilarity;
            bestIndex = j;
        }
    }
    return bestIndex;
}

static vector<Mat> yoloDetector(Mat &frame, Net &net)
{
    int ht = frame.rows;
    int wt = frame.cols;

    int length = max(ht, wt);

    Mat image = Mat::zeros(Size(length, length), frame.type());

    frame.copyTo(image(Rect(0, 0, wt, ht)));

    // Calculate the scale
    double norm_scale = static_cast<double>(length) / yoloWidth;

    Mat blob;
    blobFromImage(image, blob, yoloScale, Size(yoloWidth, yoloHeight), Scalar(), yoloSwapRB, false, CV_32F);
    net.setInput(blob);

    vector<Mat> outputs;
    net.forward(outputs);
    Mat reshapedMatrix = outputs[0].reshape(0, 84);  // Reshape to 2D (84 rows, 8400 columns)

    Mat outputTransposed;
    transpose(reshapedMatrix, outputTransposed);

    int rows = outputTransposed.rows;

    vector<Rect2d> boxes;
    vector<float> scores;
    vector<int> class_ids;

    for (int i = 0; i < rows; i++) {
        double minScore, maxScore;
        Point minClassLoc, maxClassLoc;
        minMaxLoc(outputTransposed.row(i).colRange(4, outputTransposed.cols), &minScore, &maxScore, &minClassLoc, &maxClassLoc);

        if (maxScore >= 0.25 && maxClassLoc.x == 0) {
            double centerX = outputTransposed.at<float>(i, 0);
            double centerY = outputTransposed.at<float>(i, 1);
            double w = outputTransposed.at<float>(i, 2);
            double h = outputTransposed.at<float>(i, 3);

            Rect2d box(
                centerX - 0.5 * w, // x
                centerY - 0.5 * h, // y
                w, // width
                h // height
            );
            boxes.push_back(box);
            scores.push_back(static_cast<float>(maxScore));
            class_ids.push_back(maxClassLoc.x); // x location gives the index
        }
    }

    // Apply Non-Maximum Suppression
    vector<int> indexes;
    NMSBoxes(boxes, scores, 0.25f, 0.45f, indexes, 0.5f, 0);

    vector<Mat> images;
    for (int index : indexes) {
        int x = static_cast<int>(round(boxes[index].x * norm_scale));
        int y = static_cast<int>(round(boxes[index].y * norm_scale));
        int w = static_cast<int>(round(boxes[index].width * norm_scale));
        int h = static_cast<int>(round(boxes[index].height * norm_scale));
        // Make sure the box is within the frame
        x = max(0, x);
        y = max(0, y);
        w = min(w, frame.cols - x);
        h = min(h, frame.rows - y);

        // Crop the image
        Rect roi(x, y, w, h); // Define a region of interest
        Mat crop_img = frame(roi); // Crop the region from the frame
        images.push_back(crop_img);
        imgDict[crop_img] = roi;
    }
    return images;
}

static void drawRectangle(int event, int x, int y, int, void* param) {
    Mat& img = *(Mat*)param;

    switch (event) {
        case EVENT_LBUTTONDOWN:
            drawing = true;
            ix = x;
            iy = y;
            break;

        case EVENT_MOUSEMOVE:
            if (drawing) {
                Mat img_copy = img.clone();
                rectangle(img_copy, Point(ix, iy), Point(x, y), Scalar(0, 255, 0), 2);
                imshow("TRACKING", img_copy);
            }
            break;

        case EVENT_LBUTTONUP:
            drawing = false;
            rect = Rect(Point(ix, iy), Point(x, y));
            rectangle(img, rect, Scalar(0, 255, 0), 2);
            imshow("TRACKING", img);
            break;
    }
}

void extractFrames(Net& reidNet) {
    VideoCapture cap;
    if (!videoPath.empty()){
        videoPath = findFile(videoPath);
        cap.open(videoPath);
    }
    else
        cap.open(0);

    if (!cap.isOpened()) {
        cerr << "Error: Video could not be opened." << endl;
        return;
    }

    Net net = readNetFromONNX(yoloPath);
    vector<Mat> queryImages;

    Mat queryImg;
    if (!queryImagePath.empty()) {
        queryImg = imread(queryImagePath);
        if (queryImg.empty()) {
            cerr << "Error: Query image could not be loaded." << endl;
            return;
        }
        queryImages.push_back(queryImg);
    } else {
        Mat firstFrame;
        cap.read(firstFrame);
        if (firstFrame.empty()) {
            cerr << "Error reading the video" << endl;
            return;
        }

        putText(firstFrame, "Draw Bounding Box on Target", Point(10, 30), FONT_HERSHEY_SIMPLEX, 0.6, Scalar(255, 0, 0), 2);
        imshow("TRACKING", firstFrame);
        setMouseCallback("TRACKING", drawRectangle, &firstFrame);

        for(;;) {
            if (rect.width > 0 && rect.height > 0) {
                queryImg = firstFrame(rect).clone();
                queryImages.push_back(queryImg);
                break;
            }
            if (waitKey(1) == 'q' || waitKey(1) == 27) {
                return;
            }
        }
    }

    Mat frame;
    vector<Mat> queryFeatures;
    extractFeatures(queryImages, reidNet, queryFeatures);

    for(;;) {
        if (!cap.read(frame) || frame.empty()) {
            break;
        }
        vector<Mat> detectedImages = yoloDetector(frame, net);

        vector<Mat> galleryFeatures;
        extractFeatures(detectedImages, reidNet, galleryFeatures);

        int match_idx = findMatching(queryFeatures[0], galleryFeatures);
        if (match_idx != -1 && static_cast<int>(detectedImages.size()) > match_idx) {
            Mat matchImg = detectedImages[match_idx];
            Rect bbox = imgDict[matchImg];
            rectangle(frame, bbox, Scalar(0, 0, 255), 2);
            putText(frame, "Target", Point(bbox.x, bbox.y - 10), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0, 0, 255), 2);
        }

        putText(frame, "Tracking", Point(10, 30), FONT_HERSHEY_SIMPLEX, 0.6, Scalar(255, 0, 0), 2);
        imshow("TRACKING", frame);
        int key = waitKey(1);
        if (key == 'q' || key == 27) {
            break;
        }
    }

    cap.release();
    destroyAllWindows();
}