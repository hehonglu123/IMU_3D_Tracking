import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from imu3d import grab_data_from_df

def dead_reckon():
    acc_array, omega_array, mag_array, time_array = grab_data_from_df(False, 'Rectangle')
    time = time_array
    aX = acc_array[:,0]
    aY = acc_array[:,1]
    aZ = acc_array[:,2]

    wX = omega_array[:,0]
    wY = omega_array[:,1]
    wZ = omega_array[:,2]


    X = [0]
    Y = [0]
    Z = [0]
    R = np.eye(3)
    vx = 0
    vy = 0
    vz = 0
    for i in range(len(time) - 2):
        dt = (time[i + 2] - time[i + 1]) * 10 **(+3)
        tx = dt * wX[i + 1]
        ty = dt * wY[i + 1]
        tz = dt * wZ[i + 1]
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(tx), -np.sin(tx)],
                       [0, np.sin(tx), np.cos(tx)]])

        Ry = np.array([[np.cos(ty), 0, np.sin(ty)],
                       [0, 1, 0],
                       [-np.sin(ty), 0, np.cos(ty)]])

        Rz = np.array([[np.cos(tz), -np.sin(tz), 0],
                       [np.sin(tz), np.cos(tz), 0],
                       [0, 0, 1]])
        R = np.matmul(Rz, R)
        R = np.matmul(Ry, R)
        R = np.matmul(Rx, R)
        R_temp = np.linalg.inv(R)
        dtt = (time[i + 1] - time[i]) * 10 ** (+3)

        dv = np.array([[dtt * aX[i]],
                       [dtt * aY[i]],
                       [dtt * (aZ[i] - aZ[0])]])
        dv_world = np.matmul(R_temp, dv)

        vx += dv_world[0][0]
        vy += dv_world[1][0]
        vz += dv_world[2][0]
        X.append(X[-1] + vx * dt)
        Y.append(Y[-1] + vy * dt)
        Z.append(Z[-1] + vz * dt)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, Z, c='r', marker='o')
    plt.xlabel('X')
    plt.ylabel('Y')

    plt.show()
    return 0

