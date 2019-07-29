import numpy as np
from dynmeans.dynmeans import DynMeans
import argparse
from data.dataset import get_dataloader
from evaluations import maxmatching
import matplotlib.pyplot as plt
import cv2
import seaborn
import argparse


def parser():
    arg_parse = argparse.ArgumentParser(description="inputs to MULTITUDE")
    arg_parse.add_argument("--dataset", type=str, default='zara01', choices=['zara01', 'zara02', 'uni', 'eth', 'hotel'],
                           help='which dataset you need to tracking?')
    arg_parse.add_argument("--verbose", default='trueLabel', type=str, help='need visualize?')
    args = arg_parse.parse_args()
    return args


def generateData(targets, num_points, clusterStdDev):
    """
    :param targets:         N x 3,  [trackid , pos_x, pos_y], where both pos_x and pos_y are normalized to [-1, 1]
    :return:
    """
    trueLabels = []
    clusterData = []

    for centerData in targets:
        label = int(centerData[0])
        center = np.array(centerData)[np.newaxis, 1:]
        num_new_data = (np.random.randint(num_points)+5)
        ang = 2.0*np.random.rand(num_new_data)* np.pi
        len = np.abs(np.random.randn(num_new_data))*clusterStdDev
        tmp_data = center + np.stack((np.cos(ang) * len, np.sin(ang) * len), axis=1)
        tmp_data = np.clip(tmp_data, -1, 1)
        trueLabels.extend([label]*num_new_data)
        clusterData.append(tmp_data)
        trueLabels.append(label)
        clusterData.append(center)
    clusterData = np.concatenate(clusterData, axis=0).flatten().tolist()
    return clusterData, trueLabels


def computeAccuracy(labels1, labels2, matchings):
    assert len(labels1) == len(labels2), "Error: computeAccuracy requires labels1/labels2 to have the same size"

    acc = 0
    for i in range(len(labels1)):
        if matchings[labels1[i]] == labels2[i]:
            acc += 1.0

    return 100.0*acc/ len(labels1)


def tracking(dataloader, verbose=''):

    clusterStdDev = 0.03        # used to simulate points around center
    nDataPerClusterPerStep = 25     # maximumdatapoints generated for each cluster per timestep

    # dynmeans_cluster super-parameters
    lam = 0.05
    T_Q = 6.8
    K_tau = 1.01
    Q = lam/T_Q
    tau = (T_Q * (K_tau - 1.0) + 1.0)/ (T_Q - 1.0)
    nRestarts = 10
    DMCluster = DynMeans(lam, Q, tau)

    cumulativeAccuracy = 0.0

    for step, data_piece in enumerate(dataloader):
        file_path = data_piece[0][0]
        targets = data_piece[1][0]    # N x 3, the first one are trackid and the last two are x, y

        clusterData, trueLabels = generateData(targets, nDataPerClusterPerStep, clusterStdDev)

        learnedLabels = DMCluster.cluster_wrapper(clusterData, nRestarts)

        # evaluate the performance
        matchings = maxmatching.getMaxMatching(list(learnedLabels), trueLabels)  # dict
        acc = computeAccuracy(learnedLabels, trueLabels, matchings)
        print('Step {} : Accuracy = {}'.format(step, acc))
        cumulativeAccuracy += acc
        # ####################################################3
        #  vis
        # ####################################################
        if verbose != '':
            img = cv2.imread('data/{}'.format(file_path))
            plt.imshow(img[:, :, [2, 1, 0]], alpha=0.5)
            clusterData_2d = np.array(clusterData).reshape(-1, 2)
            colors_list = list(seaborn.xkcd_rgb.keys())
            if verbose=='trueLabel':
                colors = [seaborn.xkcd_rgb[colors_list[int(k)%200*3]] for k in trueLabels]
            else:
                colors = [seaborn.xkcd_rgb[colors_list[int(k) % 200 * 3]] for k in learnedLabels]
            plt.scatter((1 + clusterData_2d[:, 1]) * img.shape[1]/2.0,
                        (1 + clusterData_2d[:, 0]) * img.shape[0]/2.0,
                        marker='x', label='true label', linewidths=2, color=colors)
            plt.pause(3)
            plt.cla()
        # ####################################################
    print('Average Accuracy: {}%'.format(cumulativeAccuracy*1.0/len(dataloader)))
    print('Done!')

if __name__=='__main__':
    args = parser()
    dataloader = get_dataloader(args.dataset)
    tracking(dataloader, verbose=args.verbose)
