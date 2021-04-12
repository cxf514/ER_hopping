
# author=cxf
# date=2020-8-8
# file for take a look

# subplot1 RNCU distribution
import numpy as np
import matplotlib.gridspec as mg
import matplotlib.pyplot as mp
# import warnings filter
from warnings import simplefilter

# ignore all future warnings
simplefilter(action='ignore')

import argparse

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('-i', type=str, default=None)
parser.add_argument('-v', type=str, default=None)
args = parser.parse_args()

###########################################################################
data1 = []
name = []
RNCU_value=[]
# Read RNCU index
with open('machine_X_index.txt', 'r') as fx1:
    for line in fx1.readlines():
        each_sample = np.array(line[0:-1].split(','))[1:].astype('int32')
        data1.append(each_sample)
        name.append(np.array(line[0:-1].split(','))[0])
x1 = np.array(data1)
# Read RNCU value
with open('machine_X_values.txt', 'r') as fy1:
    for line in fy1.readlines():
        each_sample = np.array(line[0:-1].split(','))[1:].astype('int32')
        RNCU_value.append(each_sample)
y1 = np.array(RNCU_value)


# draw pictures
for i in range(len(name)):
    # 1
    mp.figure(figsize=(10,5))
    mp.grid(ls=':')
    mp.title(f'{name[i]}')
    mp.plot(x1[i][:40],y1[i][:40], label='RNCU')
    ax = mp.gca()
    ax.xaxis.set_major_locator(mp.MultipleLocator(1))
    mp.legend()
    mp.savefig(f'distribution/{name[i]}.png')

