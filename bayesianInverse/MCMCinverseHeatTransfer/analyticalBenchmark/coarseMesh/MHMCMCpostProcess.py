import h5py
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size/2:]

filename = "MHresults"

with h5py.File(filename, "r") as f:
    # List all groups
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]

    # Get the data
    data = list(f[a_group_key])
    p = plt.figure(1,figsize=(12,8))
    plt.grid()
    for samp in f.get('samples'):
    	plt.plot(samp)
    
   # g = plt.figure(13,figsize=(12,8))
   # plt.grid()
   # for LogTarget in f.get("LogTarget"):
   # 	plt.plot(LogTarget)
    
    g = plt.figure(11,figsize=(12,8))
    plt.grid()
    for samp in f.get('samples'):
        mean = samp * 0.0
	for idx in range(len(samp)):
		mean[idx] = samp[:idx].mean() 
    	plt.plot(mean)
    
    g = plt.figure(10,figsize=(12,8))
    plt.grid()
    for samp in f.get('samples'):
    	plt.plot(autocorr(samp))

    d = plt.figure(2,figsize=(12,8))
    sns.distplot(f.get('samples')[(0)]);
    plt.grid()

    d = plt.figure(3,figsize=(12,8))
    sns.distplot(f.get('samples')[(1)]);
    plt.grid()
    plt.show()
