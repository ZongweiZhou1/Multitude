import numpy as np
import random
from time import time
from dynmeans_py.dynmeans import DynMeans
from utils.matching import matching


class DMeans():
    def __init__(self, v_lambda, T_Q, K_tau, nRestarts=10, match_func=matching):
        Q = v_lambda / T_Q
        tau = (T_Q * (K_tau - 1.0) + 1.0) / (T_Q - 1.0)
        self.dynmeans = DynMeans(v_lambda, Q, tau)

        if nRestarts <= 0:
            raise ValueError('libdynmeans: ERROR: Cannot have nRestarts <=0')

        self.nRestarts = nRestarts
        self.match_func = match_func

    def cluster(self, newobs, ref_obs, verbose=False):
        """
        :param newobs:      N x 2, double
        :param ref_obs:     M x 2, double
        :param verbose:
        :return:
        """
        tStart = time()

        if len(newobs) == 0:
            raise ValueError('libdynmeans: ERROR: newobservations is empty')

        newobservations = newobs.flatten().tolist()
        ref_observations = ref_obs.flatten().tolist()
        self.dynmeans.set_data(newobservations, ref_observations)
        observations_num = len(newobs)
        randOrderings = list(range(observations_num))
        finalObj = np.inf

        if verbose:
            print("libdynmeans: Clustering {} datapoints with {}"
                  " restarts.".format(observations_num, self.nRestarts))
        for i in range(self.nRestarts):
            self.dynmeans.set_tmpVariables(observations_num)
            obj, prevobj = np.inf, np.inf
            random.shuffle(randOrderings)
            for j in range(100):
                prevobj = obj
                self.dynmeans.assignObservations(randOrderings)
                prms = np.array(self.dynmeans.first_updateParams()).reshape(-1, 2)
                mars, mrds = self.match_func(prms, ref_obs)  # input args are all 2d
                obj = self.dynmeans.setParameters(mars, mrds)
                if obj > prevobj:
                    print("Error: obj > prevobj - monotonicity violated! Check your distance/set parameter functions...")
                    print("obj: {},  prevobj:{}".format(obj, prevobj))
                    break
                if verbose:
                    print("libdymeans: Trail: [{}/{}] objective: {}".format(i+1, self.nRestarts, obj))
                    self.dynmeans.pin_debug(1)

                if obj == prevobj:
                    break

            if obj < finalObj:
                finalObj = obj
                self.dynmeans.set_finalPrms()
        if verbose:
            print('libdynmeans: Done clustering. Min Objective: {}'.format(finalObj))
            self.dynmeans.pin_debug(2)
        finalLabels = self.dynmeans.updateState()
        tTaken = time() - tStart
        return finalLabels, finalObj, tTaken


if __name__=='__main__':
    v_lambda = 0.05
    T_Q = 6.8
    K_tau = 1.01
    dmeans = DMeans(v_lambda, T_Q, K_tau)
    for i in range(100):
        newobs = np.random.rand(50, 2)
        ref_obs = np.random.rand(8, 2)
        finalLabels, finalObj, tTaken = dmeans.cluster(newobs, ref_obs, verbose=True)
        print(finalLabels)
        print(finalObj)
        print(tTaken)