# author=cxf
# date=2020-8-8
# file for extract features of RNCU distribution by FFT 

import pandas as pd
import numpy as np
import numpy.fft as nf
import scipy.interpolate as si

import argparse

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('-v', type=str, default=None)
parser.add_argument('-i', type=str, default=None)
parser.add_argument('-p', type=int, default=95)
args = parser.parse_args()
# get RNCU distribution ->x



def extract_feature(Value, Index):
    sample_list = []
    total_list = []
    freq1 = []
    freq2 = []
    freq3 = []
    freq4 = []
    freq5 = []
    freq6 = []
    freq7 = []
    freq8 = []
    freq9 = []
    freq10 = []
    pow1 = []
    pow2 = []
    pow3 = []
    pow4 = []
    pow5 = []
    pow6 = []
    pow7 = []
    pow8 = []
    pow9 = []
    pow10 = []
    while True:
        value = Value.readline()
        if value:
            sample = value[:-1].split(',')[0]
            sample_list.append(sample)
            y = value[:-1].split(',')[1:]
            y = [int(i) for i in y]
            total = sum(y)
            y = [i / total for i in y]
            total_list.append(total)
            index = Index.readline()
            x = index[:-1].split(',')[1:]
            x = [int(i) for i in x]
            linear = si.interp1d(x, y, kind='linear')
            x = [i for i in range(1, int(x[-1]) + 1)]
            y = linear(x)
            complex_ary = nf.fft(y)
            y_ = nf.ifft(complex_ary)
            # get freqency of all sinusoidals
            freqs = nf.fftfreq(y_.size, 1)
            # get amplitude of all sinusoidals
            pows = np.abs(complex_ary)
            # get sinusoidal with top ten amplitude and record their frequency and amplitude
            freq_top_ten = freqs[freqs > 0][np.argsort(-pows[freqs > 0])][0:10]
            pow_top_ten = -np.sort(-pows[freqs > 0])[0:10]
            freq_top_ten = list(freq_top_ten)
            pow_top_ten = list(pow_top_ten)
            for i in range(10 - len(pow_top_ten)):
                pow_top_ten.append(0)
                freq_top_ten.append(0)
            freq1.append(freq_top_ten[0])
            freq2.append(freq_top_ten[1])
            freq3.append(freq_top_ten[2])
            freq4.append(freq_top_ten[3])
            freq5.append(freq_top_ten[4])
            freq6.append(freq_top_ten[5])
            freq7.append(freq_top_ten[6])
            freq8.append(freq_top_ten[7])
            freq9.append(freq_top_ten[8])
            freq10.append(freq_top_ten[9])
            pow1.append(pow_top_ten[0])
            pow2.append(pow_top_ten[1])
            pow3.append(pow_top_ten[2])
            pow4.append(pow_top_ten[3])
            pow5.append(pow_top_ten[4])
            pow6.append(pow_top_ten[5])
            pow7.append(pow_top_ten[6])
            pow8.append(pow_top_ten[7])
            pow9.append(pow_top_ten[8])
            pow10.append(pow_top_ten[9])
        else:
            break
    distribution = pd.DataFrame()
    distribution['sample'] = sample_list
    distribution['freq1'] = freq1
    distribution['freq2'] = freq2
    distribution['freq3'] = freq3
    distribution['freq4'] = freq4
    distribution['freq5'] = freq5
    distribution['freq6'] = freq6
    distribution['freq7'] = freq7
    distribution['freq8'] = freq8
    distribution['freq9'] = freq9
    distribution['freq10'] = freq10
    distribution['pow1'] = pow1
    distribution['pow2'] = pow2
    distribution['pow3'] = pow3
    distribution['pow4'] = pow4
    distribution['pow5'] = pow5
    distribution['pow6'] = pow6
    distribution['pow7'] = pow7
    distribution['pow8'] = pow8
    distribution['pow9'] = pow9
    distribution['pow10'] = pow10
    distribution['total'] = total_list
    return distribution


# from RNCU distribution of each samples, we consider there are two main distribution,
# the first whose center at 1 RNCU is caused by hopping and 
# the second is normal distribution of sequencing.
# the RNCU distance between centers of these two distribution we consider is important for cutoff setting
# the center of former we consider is always 1 RNCU
# so here we add the position of the center of latter although may be it has been included in sinusoidals
if __name__ == '__main__':
    Value = open(args.v)
    Index = open(args.i)
    df_input = extract_feature(Value,Index)
    df_input['precise']=args.p
    # write to files
    df_input.to_csv('model_input.csv', header=True, index=0)
    print('Job is done')

