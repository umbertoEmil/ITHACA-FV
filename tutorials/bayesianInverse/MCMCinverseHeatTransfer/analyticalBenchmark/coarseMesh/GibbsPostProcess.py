import h5py
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

filename = "GibbsResults"

with h5py.File(filename, "r") as f:
    # List all groups
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]

    # Get the data
    data = list(f[a_group_key])
    p = plt.figure(1,figsize=(12,8))
    plt.grid()
    for samp in f.get('samples'):
        np.delete(samp, 0)
    	plt.plot(samp)

    d = plt.figure(2,figsize=(12,8))
    sns.distplot(f.get('samples')[(0)]);
    plt.grid()

    d = plt.figure(3,figsize=(12,8))
    sns.distplot(f.get('samples')[(1)]);
    plt.grid()
    plt.show()
