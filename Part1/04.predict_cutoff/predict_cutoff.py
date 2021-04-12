# author=cxf
# date=2020-8-8
# file for model training
# import warnings filter
from warnings import simplefilter

# ignore all future warnings
simplefilter(action='ignore', category=DeprecationWarning)
import argparse
import pandas as pd
import pickle
import numpy as np
import sklearn.ensemble.forest
import sklearn.model_selection as ms
import sklearn.ensemble as se
import matplotlib.pyplot as mp
import sklearn.utils as su

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('-i', type=str, default=None)
args = parser.parse_args()



# get features and expected cutoff
df_input = pd.read_csv(args.i, index_col=0)

with open('RF_pick_cutoff.pkl', 'rb') as f:
    model=pickle.load(f)
    pred_y = model.predict(df_input)
    df_output_0=df_input[['precise']]
    df_output=df_output_0.copy()
    df_output['predict_cutoff']=pred_y
    df_output.to_csv('predict_cutoff.csv')
    print('Job is done')