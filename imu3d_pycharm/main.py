import deadreckoning as dr
import skinematics_test as skin_test
import operator
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import noctisfile

def pull_and_sort_data():
    df = pd.read_csv('Data/Self/Accelerometer.csv', sep='\s*,\s*',
                     header=0, encoding='ascii', engine='python')
    # print(acc_df)
    df = df.sort_values(by='time')

    df = df.drop_duplicates(subset='time', keep='last')

    print(df.shape[0])
    # for i in range(np.math.floor(df.shape[0])):
    #     # print(df['aX'][i])
    #     if (df['aX'].iloc[i] == 0.0 and df['aY'].iloc[i] == 0.0 and df['aZ'].iloc[i] == 0.0) or (
    #             df['wX'].iloc[i] == 0.0 and df['wY'].iloc[i] == 0.0 and df['wZ'].iloc[i] == 0.0) \
    #             or (df['mX'].iloc[i] == 0.0 and df['mY'].iloc[i] == 0.0 and df['mZ'].iloc[i] == 0.0):
    #         # print("dropped", i)
    #         df = df.drop([i])
    # ten_percent = np.math.floor(df.shape[0] * 0.50)
    # for i in range(ten_percent):
    #     df.drop(df.index[i])
    # for i in range(np.math.floor(df.shape[0])-1, np.math.floor(df.shape[0]) - ten_percent, -1):
    #     df.drop(df.index[i])

    return df


def main():
    acc_df = pull_and_sort_data()
    # noctisfile.dead_reckon()
    # plt.plot(acc_df['wX'])
    # plt.show()
    # print(acc_df[0:10])
    dr.dead_reckon(acc_df)
    skin_test.skin_dead_reckon(acc_df)
    # skin_test.rom_elbow()


main()
