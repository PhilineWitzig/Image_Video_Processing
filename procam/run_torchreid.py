"""Module for running person detection -> tracking -> re-identification"""

import numpy as np
import cv2
import config
import torch
import STE_NVAN
import imutils
import dlib
import torch.nn as nn
from torchsummary import summary
from STE_NVAN.net import models
from collections import OrderedDict


def main():

    # load pretrained STE_NVAN
    network = nn.DataParallel(models.CNN(2048, model_type='resnet50_NL_stripe_hr',
                                         num_class=625, non_layers=[0, 2, 3, 0], temporal='Done'))
    state = torch.load('STE_NVAN/ckpt/STE_NVAN.pth', map_location='cpu')
    network.load_state_dict(state)
    network = network.cpu().double()  # load to cpzs for now
    network.eval()

    # load DNN object detector
    object_detector = cv2.dnn.readNetFromCaffe(
        config.MODEL_MOBILE_NET_SSD_PROTOTEXT, config.MODEL_MOBILE_NET_SSD_CAFFE)

    # for reading frames from webcam
    cap = cv2.VideoCapture(0)
    frame_counter = 0
    att = 0

    # keep track of trackers, labels and IDs
    trackers = []
    labels = []
    IDs = []

    # buffering frames for a certain number of persons and for a certain number of frames
    buffer = np.zeros([config.MAX_PERSON_COUNT, config.ATTENTION_WIDTH, 3,
                       256, 128])  # 8 persons for a duration of 8 frames

    with torch.no_grad():
        while(True):

            if frame_counter % config.FRAME_DROP == 0:  # process every 10th frame for real time deployment
                ret, frame = cap.read()
                frame = imutils.resize(frame, width=600)
                frame_copy = frame.copy()
                rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)  # dlib requires rgb

                if len(trackers) == 0 or frame_counter % 300 == 0:  # rerun detection every 10 seconds
                    trackers = []
                    labels = []
                    IDs = []
                    # 8 persons for a duration of 8 frames
                    buffer = np.zeros(
                        [config.MAX_PERSON_COUNT, config.ATTENTION_WIDTH, 3, 256, 128])

                    (h, w) = frame.shape[:2]
                    blob = cv2.dnn.blobFromImage(frame, 0.007843, (w, h), 127.5)
                    object_detector.setInput(blob)
                    detections = object_detector.forward()

                    for i in np.arange(0, detections.shape[2]):

                        confidence = detections[0, 0, i, 2]
                        # only consider objects with confidence value larger than certain value
                        if confidence > config.CONFIDENCE:
                            idx = int(detections[0, 0, i, 1])
                            label = config.CLASSES[idx]
                            # if class label is not a person, ignore it
                            if config.CLASSES[idx] != "person":
                                continue

                            box = detections[0, 0, i, 3:7] * np.array([w, h, w, h])
                            (startX, startY, endX, endY) = box.astype("int")
                            t = dlib.correlation_tracker()
                            rect = dlib.rectangle(startX, startY, endX, endY)
                            # start tracking the object
                            t.start_track(rgb, rect)
                            labels.append(label)
                            trackers.append(t)
                            cv2.rectangle(frame_copy, (startX, startY), (endX, endY),
                                          (0, 255, 0), 2)

                else:
                    # update tracker psoitions only
                    for (t, l) in zip(trackers, labels):
                        t.update(rgb)
                        pos = t.get_position()
                        startX = int(pos.left())
                        startY = int(pos.top())
                        endX = int(pos.right())
                        endY = int(pos.bottom())

                        cv2.rectangle(frame_copy, (startX, startY), (endX, endY),
                                      (0, 255, 0), 2)

                for i in range(len(trackers)):  # iterate over all persons
                    pos = trackers[i].get_position()
                    startX = int(pos.left())
                    startY = int(pos.top())
                    endX = int(pos.right())
                    endY = int(pos.bottom())

                    person = frame[startX:endX, startY:endY]
                    person = cv2.resize(person, (128, 256))
                    buffer[i, att] = np.reshape(person, [3, 256, 128])

                att += 1
                if att % config.ATTENTION_WIDTH == 0:
                    input = torch.from_numpy(buffer)
                    output = network.forward(input.double())
                    ids = torch.max(output.data, 1)
                    # for i in range(len(output)):
                    #    if output[i] != 236:  # 236 seems to be the label for 0 input
                    #        IDs.append(output[i])
                    print(ids)
                    att = 0

                # TODO: display IDs in frame
                cv2.imshow('frame', frame_copy)
                if cv2.waitKey(1) & 0xFF == ord('q'):
                    break

            frame_counter += 1
            if frame_counter == 10000:  # start counting from 0 again after the 10k'th frame to avoid overflow
                frame_counter = 0

    cap.release()
    cv2.destroyAllWindows()


if __name__ == "__main__":
    main()
