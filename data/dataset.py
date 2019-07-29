import torch
import numpy as np
import json
import torch.utils.data as data


class Traj_dataset(data.Dataset):
    def __init__(self, dataset, n_noisy=2):
        self.dataset = json.load(open('data/{}.json'.format(dataset)))
        # add some noise trajectories
        for i in range(n_noisy):
            self.generate_noisy_traj(1000+i)

    def generate_noisy_traj(self, tid):
        traj_len = np.random.randint(100)
        traj_start_point = np.random.rand(1, 2) * 2 - 1.0
        traj_end_point = (np.random.rand(1, 2) * 2 - 1.0) * 0.2 + traj_start_point
        delta_v = (traj_end_point - traj_start_point) / traj_len * np.arange(traj_len).reshape(-1, 1) + traj_start_point
        timestep_point = np.concatenate((np.ones((traj_len, 1))*tid, delta_v), axis=1)
        traj_start_timestep = np.random.randint(len(self.dataset) - 100)

        for timestep in range(traj_start_timestep, traj_start_timestep + traj_len):
            self.dataset[timestep][1] = np.concatenate((self.dataset[timestep][1],
                                                        timestep_point[timestep - traj_start_timestep].reshape(1, -1)),
                                                       axis=0)

    def __len__(self):
        return len(self.dataset)

    def __getitem__(self, item):
        return self.dataset[item][0], self.dataset[item][1]

    @staticmethod
    def collate_fn(batch):
        img_path_list, objects_list = zip(*batch)
        return img_path_list, objects_list



def get_dataloader(dataset_name):
    dataset = Traj_dataset(dataset_name)
    dataloader = data.DataLoader(dataset, batch_size=1, shuffle=False, pin_memory=True, collate_fn=dataset.collate_fn,
                                num_workers=20, drop_last=False)

    return dataloader


if __name__=='__main__':
    dataloader = get_dataloader('zara01')
    for step, v in enumerate(dataloader):
        print(v[0][0])
        print(v[1][0])
        break
