import csv

import skinematics as skin
import numpy as np
from skinematics.sensors.manual import MyOwnSensor
import matplotlib.pyplot as plt

RATE = 1000000000
def skin_dead_reckon():
    # load data
    f = open('Data/Rectangle/Accelerometer.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time = [float(i) for i in T[0]]
    aX = [float(i) for i in T[1]]
    aY = [float(i) for i in T[2]]
    aZ = [float(i) for i in T[3]]
    f = open('Data/Rectangle/Gyroscope.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    wX = [float(i) for i in T[1]]
    wY = [float(i) for i in T[2]]
    wZ = [float(i) for i in T[3]]
    f = open('Data/Rectangle/MagneticField.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    mX = [float(i) for i in T[1]]
    mY = [float(i) for i in T[2]]
    mZ = [float(i) for i in T[3]]

    # determine which of the sensors has the longest array
    acc_len = len(aX)
    omega_len = len(wX)
    mag_len = len(mX)

    max_len = max([acc_len, omega_len, mag_len])

    # interpolate all of these data points so they match up
    aX_interp = np.interp(np.linspace(0, 1, max_len), np.linspace(0, 1, len(aX)), aX)
    aY_interp = np.interp(np.linspace(0, 1, max_len), np.linspace(0, 1, len(aY)), aY)
    aZ_interp = np.interp(np.linspace(0, 1, max_len), np.linspace(0, 1, len(aZ)), aZ)
    wX_interp = np.interp(np.linspace(0, 1, max_len), np.linspace(0, 1, len(wX)), wX)
    wY_interp = np.interp(np.linspace(0, 1, max_len), np.linspace(0, 1, len(wY)), wY)
    wZ_interp = np.interp(np.linspace(0, 1, max_len), np.linspace(0, 1, len(wZ)), wZ)
    mX_interp = np.interp(np.linspace(0, 1, max_len), np.linspace(0, 1, len(mX)), mX)
    mY_interp = np.interp(np.linspace(0, 1, max_len), np.linspace(0, 1, len(mY)), mY)
    mZ_interp = np.interp(np.linspace(0, 1, max_len), np.linspace(0, 1, len(mZ)), mZ)

    # create N x 3 np arrays
    acc_array = np.zeros((max_len, 3))
    omega_array = np.zeros((max_len, 3))
    mag_array = np.zeros((max_len, 3))

    for i in range(len(aX)):
        acc_array[i, 0] = aX_interp[i]
        acc_array[i, 1] = aY_interp[i]
        acc_array[i, 2] = aZ_interp[i]
    for i in range(len(aX)):
        omega_array[i, 0] = wX_interp[i]
        omega_array[i, 1] = wY_interp[i]
        omega_array[i, 2] = wZ_interp[i]
    for i in range(len(mX)):
        mag_array[i, 0] = mX_interp[i]
        mag_array[i, 1] = mY_interp[i]
        mag_array[i, 2] = mZ_interp[i]

    # format as a dictionary
    in_data = {
        'rate': RATE,
        'acc': acc_array,
        'omega': omega_array,
        'mag': mag_array
    }

    # create your own sensor (not used here)
    my_sensor = MyOwnSensor(in_data=in_data)

    print(np.shape(np.eye(3)), np.shape(acc_array), np.shape(omega_array), np.shape(np.zeros((3, 1))))

    # dead reckon
    q1, pos1, vel1 = skin.imus.analytical(R_initialOrientation=np.eye(3), \
                                          omega=omega_array, initialPosition=np.array([0,0,0]), \
                                          accMeasured=acc_array, rate=1000000000)
    # plot
    print(np.shape(pos1))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pos1[:, 0], pos1[:, 1], pos1[:, 2], c='r', marker='o')
    plt.show()
