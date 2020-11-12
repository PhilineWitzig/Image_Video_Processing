import numpy as np
import cv2
import config
import torch
import STE_NVAN
import imutils
import torch.nn as nn
from torchsummary import summary
from STE_NVAN.net import models
from collections import OrderedDict
from imutils.object_detection import non_max_suppression


def detect_pedestrian(img):
    hog = cv2.HOGDescriptor()
    hog.setSVMDetector(cv2.HOGDescriptor_getDefaultPeopleDetector())
    # detect people in the image
    (rects, weights) = hog.detectMultiScale(img, winStride=(4, 4),
                                            padding=(8, 8), scale=1.05)

    rects = np.array([[x, y, x + w, y + h] for (x, y, w, h) in rects])
    pick = non_max_suppression(rects, probs=None, overlapThresh=0.65)

    return pick


def main():

    network = nn.DataParallel(models.CNN(2048, model_type='resnet50_s1', num_class=625))
    # summary(network, (256, 128, 3))

    state_dict = torch.load(config.MODEL_WEIGHT_RESNET50_s1,
                            map_location=torch.device('cpu'))
    network.load_state_dict(state_dict)

    cap = cv2.VideoCapture(0)
    frame_counter = 0

    while(True):
        if frame_counter % 5 == 0:  # process every 5th frame
            ret, frame = cap.read()
            frame = imutils.resize(frame, width=min(400, frame.shape[1]))
            rects = detect_pedestrian(frame)
            # draw the final bounding boxes
            for (xA, yA, xB, yB) in rects:
                cv2.rectangle(frame, (xA, yA), (xB, yB), (0, 255, 0), 2)
                #person = frame[yA:yB, xA:xB]
                #person = cv2.resize(person, (128, 256))

                #input = torch.from_numpy(person)

                # TOOO fix
                #output = network.forward(input)
                #ff = outputs.data.cpu()

            # Display the resulting frame
            cv2.imshow('frame', frame)
            if cv2.waitKey(1) & 0xFF == ord('q'):
                break
        else:
            pass

        frame_counter += 1

        if frame_counter == 10000:  # start counting from 0 again after the 10k'th frame
            frame_counter = 0

    cap.release()
    cv2.destroyAllWindows()


if __name__ == "__main__":
    main()
