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


def main():

    network = nn.DataParallel(models.CNN(2048, model_type='resnet50_NL_stripe_hr',
                                         num_class=625, non_layers=[0, 2, 3, 0], temporal='Done'))
    # load weights
    state = torch.load('STE_NVAN/ckpt/STE_NVAN.pth', map_location='cpu')
    network.load_state_dict(state)
    network = network.cpu().double()
    network.eval()
    # summary(network, (256, 128, 3))

    # state_dict = torch.load(config.MODEL_WEIGHT_RESNET50_s1,
    #                        map_location=torch.device('cpu'))
    # network.load_state_dict(state_dict)

    hog = cv2.HOGDescriptor()
    hog.setSVMDetector(cv2.HOGDescriptor_getDefaultPeopleDetector())

    cap = cv2.VideoCapture(0)
    frame_counter = 0
    buffer = np.zeros([8, 3, 256, 128])
    i = 0

    with torch.no_grad():
        while(True):

            if frame_counter % 3 == 0:  # process every 3rd frame
                ret, frame = cap.read()
                #frame = imutils.resize(frame, width=min(400, frame.shape[1]))
                #frame = cv2.resize(frame, (256, 128))
                # grayscale for faster detection?
                # gray = cv2.cvtColor(frame, cv2.COLOR_RGB2GRAY)
                # detect people in the image
                rects, weights = hog.detectMultiScale(frame, winStride=(4, 4),
                                                      padding=(8, 8), scale=1.05)

                rects = np.array([[x, y, x + w, y + h] for (x, y, w, h) in rects])
                rects = non_max_suppression(rects, probs=None, overlapThresh=0.65)
                # draw the final bounding boxes
                for (xA, yA, xB, yB) in rects:
                    cv2.rectangle(frame, (xA, yA), (xB, yB), (0, 255, 0), 2)
                    person = frame[yA:yB, xA:xB]
                    person = cv2.resize(person, (128, 256))
                    print("Found person.")
                    buffer[i] = np.reshape(person, [3, 256, 128])
                    i += 1
                    if i % 8 == 0:
                        i = 0
                        # add batch dimension
                        input = torch.from_numpy(buffer)
                        B, C, H, W = input.shape
                        input = input.reshape(B//8, 8, C, H, W)  # TODO: don't hardcode 8
                        # TOOO fix
                        output = network.forward(input.double())
                        print(output.size())
                        print(torch.max(output.data, 1))
                        #ff = outputs.data.cpu()

                # Display the resulting frame
                cv2.imshow('frame', frame)
                if cv2.waitKey(1) & 0xFF == ord('q'):
                    break

            frame_counter += 1

            if frame_counter == 10000:  # start counting from 0 again after the 10k'th frame
                frame_counter = 0

    cap.release()
    cv2.destroyAllWindows()


if __name__ == "__main__":
    main()
