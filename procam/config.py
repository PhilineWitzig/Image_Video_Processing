"""This module contains all configuration settings for reID module"""

import os


MODEL_WEIGHT_RESNET50_s1 = os.path.join(os.getcwd(), "STE_NVAN", "ckpt", "R50_baseline_mean.pth")
MODEL_MOBILE_NET_SSD_PROTOTEXT = os.path.join(
    os.getcwd(), "models", "mobilenet_ssd", "MobileNetSSD_deploy.prototxt")
MODEL_MOBILE_NET_SSD_CAFFE = os.path.join(
    os.getcwd(), "models", "mobilenet_ssd", "MobileNetSSD_deploy.caffemodel")

# classes for mobile net ssd classification
CLASSES = ["background", "aeroplane", "bicycle", "bird", "boat",
           "bottle", "bus", "car", "cat", "chair", "cow", "diningtable",
           "dog", "horse", "motorbike", "person", "pottedplant", "sheep",
           "sofa", "train", "tvmonitor"]

ATTENTION_WIDTH = 8  # determines the attention width of the STE_NVAN neural network
MAX_PERSON_COUNT = 8  # maximum number of persons which can be in a frame
CONFIDENCE = 0.7  # confidence value for cv2 object detector
FRAME_DROP = 10  # process every FRAME_DROP'th frame for real-time deployment
