import numpy as np
import maxmatching

labels1 = [1,2,3,4,5,6,5,4,3,1]
labels2 = [2,1,2,3,4,6,5,4,3,2]

matching_map = maxmatching.getMaxMatching(labels1, labels2)

print(isinstance(matching_map, maxmatching.map_int2int))

for k, v in matching_map.items():
    print('{} : {}'.format(k, v))
