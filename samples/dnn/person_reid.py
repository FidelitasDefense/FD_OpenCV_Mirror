#!/usr/bin/env python
'''
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

        python person_reid.py --query=/path/to/query_image(optional) --video=/path/to/video/footage --model=path/to/youtu_reid_baseline_medium.onnx --yolo=/path/to/yolov8.onnx

    You can download a baseline ReID model from:
        https://github.com/ReID-Team/ReID_extra_testdata

    Drive Link: https://drive.google.com/drive/folders/1wFGcuolSzX3_PqNKb4BAV3DNac7tYpc2

'''
import argparse
import os.path
import numpy as np
import cv2 as cv
from common import *

def help():
    print(
        '''
        Use this script for Person Re-identification using OpenCV.

        Firstly, download required models i.e. reid and yolov8 using `download_models.py` (if not already done). Set environment variable OPENCV_DOWNLOAD_CACHE_DIR to specify where models should be downloaded. Also, point OPENCV_SAMPLES_DATA_PATH to opencv/samples/data.

        To run:
        Example: python person_reid.py reid

        Re-identification model path can also be specified using --model argument and detection model can be specified using --yolo_model argument.
        '''
    )

def get_args_parser():
    backends = ("default", "openvino", "opencv", "vkcom", "cuda")
    targets = ("cpu", "opencl", "opencl_fp16", "ncs2_vpu", "hddl_vpu", "vulkan", "cuda", "cuda_fp16")

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--zoo', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models.yml'),
                        help='An optional path to file with preprocessing parameters.')
    parser.add_argument('--query', '-q', help='Path to target image. Skip this argument to select target in the video frame.')
    parser.add_argument('--input', '-i', default=0, help='Path to video file.', required=False)
    parser.add_argument('--backend', default="default", type=str, choices=backends,
            help="Choose one of computation backends: "
            "default: automatically (by default), "
            "openvino: Intel's Deep Learning Inference Engine (https://software.intel.com/openvino-toolkit), "
            "opencv: OpenCV implementation, "
            "vkcom: VKCOM, "
            "cuda: CUDA, "
            "webnn: WebNN")
    parser.add_argument('--target', default="cpu", type=str, choices=targets,
            help="Choose one of target computation devices: "
            "cpu: CPU target (by default), "
            "opencl: OpenCL, "
            "opencl_fp16: OpenCL fp16 (half-float precision), "
            "ncs2_vpu: NCS2 VPU, "
            "hddl_vpu: HDDL VPU, "
            "vulkan: Vulkan, "
            "cuda: CUDA, "
            "cuda_fp16: CUDA fp16 (half-float preprocess)")
    args, _ = parser.parse_known_args()
    add_preproc_args(args.zoo, parser, 'person_reid', prefix="")
    add_preproc_args(args.zoo, parser, 'person_reid', prefix="yolo_")
    parser = argparse.ArgumentParser(parents=[parser],
                                        description='Person Re-identification using OpenCV.',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    return parser.parse_args()

drawing = False
ix, iy = -1, -1
rect = None
img_dict = {} # Dictionary to store bounding boxes for corresponding cropped image

def yolo_detector(frame, net):
    global img_dict
    height, width, _ = frame.shape

    length = max((height, width))
    image = np.zeros((length, length, 3), np.uint8)
    image[0:height, 0:width] = frame

    scale = length/args.yolo_width
    # Create blob from the frame with correct scale factor and size for the model

    blob = cv.dnn.blobFromImage(image, scalefactor=args.yolo_scale, size=(args.yolo_width, args.yolo_height), swapRB=args.yolo_rgb)
    net.setInput(blob)
    outputs = net.forward()

    outputs = np.array([cv.transpose(outputs[0])])
    rows = outputs.shape[1]

    boxes = []
    scores = []
    class_ids = []

    for i in range(rows):
        classes_scores = outputs[0][i][4:]
        (_, maxScore, _, (x, maxClassIndex)) = cv.minMaxLoc(classes_scores)
        if maxScore >= 0.25:
            box = [
                outputs[0][i][0] - (0.5 * outputs[0][i][2]),
                outputs[0][i][1] - (0.5 * outputs[0][i][3]),
                outputs[0][i][2],
                outputs[0][i][3],
            ]
            boxes.append(box)
            scores.append(maxScore)
            class_ids.append(maxClassIndex)


    # Apply Non-Maximum Suppression
    indexes = cv.dnn.NMSBoxes(boxes, scores, 0.25, 0.45, 0.5)

    images = []
    for i in indexes:
        x, y, w, h = boxes[i]
        x = round(x*scale)
        y = round(y*scale)
        w = round(w*scale)
        h = round(h*scale)

        x, y = max(0, x), max(0, y)
        w, h = min(w, frame.shape[1] - x), min(h, frame.shape[0] - y)
        crop_img = frame[y:y+h, x:x+w]
        images.append(crop_img)
        img_dict[crop_img.tobytes()] = (x, y, w, h)
    return images

def draw_rectangle(event, x, y, flags, param):
    global ix, iy, drawing, rect

    if event == cv.EVENT_LBUTTONDOWN:
        drawing = True
        ix, iy = x, y

    elif event == cv.EVENT_MOUSEMOVE:
        if drawing:
            img_copy = param[0].copy()
            cv.rectangle(img_copy, (ix, iy), (x, y), (0, 255, 0), 2)
            cv.imshow('TRACKING', img_copy)

    elif event == cv.EVENT_LBUTTONUP:
        drawing = False
        rect = (ix, iy, x, y)
        cv.rectangle(param[0], (ix, iy), (x, y), (0, 255, 0), 2)
        cv.imshow('TRACKING', param[0])

def extract_frames(query_image_path, model_path, yolo_path, batch_size=32):
    cap = cv.VideoCapture(cv.samples.findFile(args.input) if args.input else 0)
    net = cv.dnn.readNet(yolo_path)
    query_images = []

    if query_image_path:
        query_images = [cv.imread(query_image_path)]
    else:
        ret, first_frame = cap.read()
        if not ret:
            print("Error reading the video")
            return
        cv.putText(first_frame, "Draw Bounding Box on Target", (10, 30), cv.FONT_HERSHEY_SIMPLEX, 0.6, (255, 0, 0), 2)
        cv.imshow('TRACKING', first_frame)
        cv.setMouseCallback('TRACKING', draw_rectangle, [first_frame])

        while True:
            if rect:
                break
            if cv.waitKey(1) & 0xFF == ord('q'):
                return

        x1, y1, x2, y2 = rect
        query_image = first_frame[y1:y2, x1:x2]
        query_images = [query_image]

    while cap.isOpened():
        ret, frame = cap.read()
        if not ret:
            break
        images = yolo_detector(frame, net)
        query_feat = extract_feature(query_images, model_path, batch_size)
        gallery_feat = extract_feature(images, model_path, batch_size)

        match_idx = find_matching(query_feat, gallery_feat)

        match_img = images[match_idx]
        x, y, w, h = img_dict[match_img.tobytes()]
        cv.rectangle(frame, (x, y), (x + w, y + h), (0, 0, 255), 2)
        cv.putText(frame, "Target", (x, y - 10), cv.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 255), 2)

        cv.putText(frame, "Tracking", (10, 30), cv.FONT_HERSHEY_SIMPLEX, 0.6, (255, 0, 0), 2)
        cv.imshow("TRACKING", frame)
        if cv.waitKey(1) & 0xFF in [ord('q'), 27]:
            break

    cap.release()
    cv.destroyAllWindows()
    return

def extract_feature(images, model_path, batch_size = 32):
    """
    Extract features from images
    :param images: the input images
    :param model_path: path to ReID model
    :param batch_size: the batch size for each network inference iteration
    """
    feat_list = []

    for i in range(0, len(images), batch_size):
        batch = images[i : min(i + batch_size, len(images))]
        blob = cv.dnn.blobFromImages(batch, scalefactor=args.scale, size=(args.width, args.height), mean=args.mean, swapRB=args.rgb, crop=False, ddepth=cv.CV_32F)

        for i in range(blob.shape[1]):
            blob[:, i, :, :] /= args.std[i]

        feat = run_net(blob, model_path)

        feat_list.append(feat)

    feats = np.concatenate(feat_list, axis = 0)
    return feats

def run_net(inputs, model_path):
    """
    Forword propagation for a batch of images.
    :param inputs: input batch of images
    :param model_path: path to ReID model
    """
    net = cv.dnn.readNet(model_path)
    net.setPreferableBackend(get_backend_id(args.backend))
    net.setPreferableTarget(get_target_id(args.target))
    net.setInput(inputs)
    out = net.forward()
    out = np.reshape(out, (out.shape[0], out.shape[1]))
    return out

def find_matching(query_feat, gallery_feat):
    """
    Return the index of the gallery image most similar to the query image
    :param query_feat: array of feature vectors of query images
    :param gallery_feat: array of feature vectors of gallery images
    """
    cv.normalize(query_feat, query_feat, 1.0, 0.0, cv.NORM_L2)
    cv.normalize(gallery_feat, gallery_feat, 1.0, 0.0, cv.NORM_L2)

    sim = query_feat.dot(gallery_feat.T)
    index = np.argmax(sim, axis=1)[0]
    return index

if __name__ == '__main__':
    args = get_args_parser()
    if args.alias is None or hasattr(args, 'help'):
        help()
        exit(1)

    args.model = findModel(args.model, args.sha1)
    if not os.path.isfile(args.model):
        raise OSError("Model not exist")

    if args.yolo_model is None:
        print("[ERROR] Please pass path to yolov8.onnx model file using --yolo_model.")
        exit(1)
    else:
        args.yolo_model = findModel(args.yolo_model, args.yolo_sha1)
    extract_frames(args.query, args.model, args.yolo_model)