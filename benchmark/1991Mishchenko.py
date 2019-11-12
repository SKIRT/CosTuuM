import numpy as np

fname = "1991Mishchenko_prolate_pdg.txt"

rawdata = np.loadtxt(fname)

data = np.zeros((rawdata.shape[0], rawdata.shape[1]+1))
data[:,1:6] = rawdata[:,:]
data[::4,0] = 0.
data[1::4,0] = np.pi / 6.
data[2::4,0] = np.pi / 3.
data[3::4,0] = 0.5 * np.pi

print(data)
