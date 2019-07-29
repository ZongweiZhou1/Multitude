import dynmeans
import numpy as np

v_lambda = 0.05
T_Q = 6.8
K_tau = 1.01
Q = v_lambda / T_Q
tau = (T_Q * (K_tau - 1.0) + 1.0) / (T_Q - 1.0)
nRestarts = 10

object_d = dynmeans.DynMeans(v_lambda, Q, tau)

clusterCenters = np.random.rand(10, 2)
centerData = []
for center in clusterCenters:
    newdata = center + (np.random.rand(2) - 0.5) * 0.2
    centerData.extend(newdata.tolist())

learnedLabels = object_d.cluster_wrapper(centerData, nRestarts)
print(learnedLabels)
