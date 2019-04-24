import csv

import skinematics as skin
import numpy as np
from skinematics.sensors.manual import MyOwnSensor
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.signal import butter, lfilter
from scipy.signal import freqs

RATE = 1000000000
GRAVITY = np.r_[0, 0, 9.8]

def rom_elbow():
    l_elbow = 50
    print("rom elbow")

    point_cloud_elbow = []
    for theta1 in np.linspace(np.deg2rad(-60), np.deg2rad(180), 20):
        for theta2 in np.linspace(np.deg2rad(-40), np.deg2rad(120), 20):
            # theta1 = np.deg2rad(theta1)
            # theta2 = np.deg2rad(theta2)
            # if np.cos(theta2) * np.sin(theta1) > 0:
            point = ([
                    np.cos(theta2) * np.sin(theta1),
                    np.sin(theta2),
                    -np.cos(theta1) * np.cos(theta2)
                    ])
            # point = point * l_elbow

            point_cloud_elbow.append(point)
    point_cloud_elbow = np.array(point_cloud_elbow)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(point_cloud_elbow[:, 0], point_cloud_elbow[:, 1], point_cloud_elbow[:, 2], c='g', marker='*')
    # plt.show()

    A1 = np.zeros((4, 4))
    A2 = np.zeros((4, 4))

    point_cloud_elbow2 = []
    # np.linspace(np.deg2rad(-180), np.deg2rad(180)):
    for t1 in np.linspace(np.deg2rad(-60), np.deg2rad(180), 20):
        for t2 in np.linspace(np.deg2rad(-40), np.deg2rad(120), 20):
            theta1 = t1
            theta2 = t2-t1
            a1 = 0
            # alpha1 = np.pi/2
            alpha1 = np.pi/2
            d1 = 0
            A1[0,0] = np.cos(theta1)
            A1[0,1] = -np.sin(theta1)* np.cos(alpha1)
            A1[0,2] = np.sin(theta1)* np.sin(alpha1)
            A1[0,3] = a1 * np.cos(theta1)
            A1[1,0] = np.sin(theta1)
            A1[1,1] = np.cos(theta1) * np.cos(alpha1)
            A1[1,2] = -np.cos(theta1) * np.sin(alpha1)
            A1[1,3] = a1 * np.sin(theta1)
            A1[2,0] = 0
            A1[2,1] = np.sin(alpha1)
            A1[2,2] = np.cos(alpha1)
            A1[2,3] = d1
            A1[3,0] = 0
            A1[3,1] = 0
            A1[3,2] = 0
            A1[3,3] = 1

            a2 = 0
            alpha2 = -np.pi/2
            d2 = 0
            A2[0, 0] = np.cos(theta2)
            A2[0, 1] = -np.sin(theta2)* np.cos(alpha2)
            A2[0, 2] = np.sin(theta2)* np.sin(alpha2)
            A2[0, 3] = a2 * np.cos(theta2)
            A2[1, 0] = np.sin(theta2)
            A2[1, 1] = np.cos(theta2) * np.cos(alpha2)
            A2[1, 2] = -np.cos(theta2) * np.sin(alpha2)
            A2[1, 3] = a2 * np.sin(theta2)
            A2[2, 0] = 0
            A2[2, 1] = np.sin(alpha2)
            A2[2, 2] = np.cos(alpha2)
            A2[2, 3] = d2
            A2[3, 0] = 0
            A2[3, 1] = 0
            A2[3, 2] = 0
            A2[3, 3] = 1

            M = np.matmul(A1,A2)
            R = M[:3,:3]
            # print
            T = M[3:, :3]
            # print(T.shape, T, M.shape, M)
            # point_cloud_elbow2.append(T)
            base = np.r_[0,-1 * 1,0]

            # point = np.matmul(base,A)


            # if np.matmul(base,R)[2] > 0:
            point_cloud_elbow2.append(np.matmul(base,R))
                # print(M[3:,:3].shape, R.shape)
                # T = np.r_[M[3,0],M[3,1],M[3,2]]
                # point_cloud_elbow2.append(np.matmul(T, R))
    point_cloud_elbow2 = np.array(point_cloud_elbow2)
    # print(point_cloud_elbow.shape, point_cloud_elbow2.shape)
    ax.scatter(point_cloud_elbow2[:, 0], point_cloud_elbow2[:, 1], point_cloud_elbow2[:, 2], c='r', marker='o')
    ax.set_xlabel('$X$', fontsize=10)
    ax.set_ylabel('$Y$', fontsize=10)
    ax.set_zlabel('$Z$', fontsize=10)
    plt.show()

    # rom_wrist(point_cloud_elbow2)

def rom_wrist(point_cloud_elbow):
    l_forearm = 50
    point_cloud_wrist = []

    point_cloud_wrist = np.array(point_cloud_wrist)

    # dir_elbow

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(point_cloud_elbow[:, 0], point_cloud_elbow[:, 1], point_cloud_elbow[:, 2], c='r', marker='*')
    # ax.scatter(point_cloud_wrist[:, 0], point_cloud_wrist[:, 1], point_cloud_wrist[:, 2], c='r', marker='*')
    plt.show()
    # for theta1 in np.linspace(np.deg2rad(-60), np.deg2rad(180), 20):
    #     for theta2 in np.linspace(np.deg2rad(-40), np.deg2rad(120), 20):
    #         for theta3 in np.linspace(np.deg2rad(-60), np.deg2rad(180), 20):
    #             for theta4 in np.linspace(np.deg2rad(-40), np.deg2rad(120), 20):
    #                 a1 = 0
    #                 # alpha1 = np.pi/2
    #                 alpha1 = np.pi / 2
    #                 d1 = 0
    #                 A1[0, 0] = np.cos(theta1)
    #                 A1[0, 1] = -np.sin(theta1) * np.cos(alpha1)
    #                 A1[0, 2] = np.sin(theta1) * np.sin(alpha1)
    #                 A1[0, 3] = a1 * np.cos(theta1)
    #                 A1[1, 0] = np.sin(theta1)
    #                 A1[1, 1] = np.cos(theta1) * np.cos(alpha1)
    #                 A1[1, 2] = -np.cos(theta1) * np.sin(alpha1)
    #                 A1[1, 3] = a1 * np.sin(theta1)
    #                 A1[2, 0] = 0
    #                 A1[2, 1] = np.sin(alpha1)
    #                 A1[2, 2] = np.cos(alpha1)
    #                 A1[2, 3] = d1
    #                 A1[3, 0] = 0
    #                 A1[3, 1] = 0
    #                 A1[3, 2] = 0
    #                 A1[3, 3] = 1
    #
    #                 a2 = 0
    #                 alpha2 = -np.pi / 2
    #                 d2 = 0
    #                 A2[0, 0] = np.cos(theta2)
    #                 A2[0, 1] = -np.sin(theta2) * np.cos(alpha2)
    #                 A2[0, 2] = np.sin(theta2) * np.sin(alpha2)
    #                 A2[0, 3] = a2 * np.cos(theta2)
    #                 A2[1, 0] = np.sin(theta2)
    #                 A2[1, 1] = np.cos(theta2) * np.cos(alpha2)
    #                 A2[1, 2] = -np.cos(theta2) * np.sin(alpha2)
    #                 A2[1, 3] = a2 * np.sin(theta2)
    #                 A2[2, 0] = 0
    #                 A2[2, 1] = np.sin(alpha2)
    #                 A2[2, 2] = np.cos(alpha2)
    #                 A2[2, 3] = d2
    #                 A2[3, 0] = 0
    #                 A2[3, 1] = 0
    #                 A2[3, 2] = 0
    #                 A2[3, 3] = 1
    #
    #                 a3 = 0
    #                 alpha3 = np.pi / 2
    #                 d3 = 0
    #                 A3[0, 0] = np.cos(theta3)
    #                 A3[0, 1] = -np.sin(theta3) * np.cos(alpha3)
    #                 A3[0, 2] = np.sin(theta3) * np.sin(alpha3)
    #                 A3[0, 3] = a3 * np.cos(theta3)
    #                 A3[1, 0] = np.sin(theta3)
    #                 A3[1, 1] = np.cos(theta3) * np.cos(alpha3)
    #                 A3[1, 2] = -np.cos(theta3) * np.sin(alpha3)
    #                 A3[1, 3] = a3 * np.sin(theta3)
    #                 A3[2, 0] = 0
    #                 A3[2, 1] = np.sin(alpha3)
    #                 A3[2, 2] = np.cos(alpha3)
    #                 A3[2, 3] = d3
    #                 A3[3, 0] = 0
    #                 A3[3, 1] = 0
    #                 A3[3, 2] = 0
    #                 A3[3, 3] = 1
    #
    #                 a4 = 0
    #                 alpha4 = 0
    #                 d4 = l_elbow
    #                 A4[0, 0] = np.cos(theta4)
    #                 A4[0, 1] = -np.sin(theta4) * np.cos(alpha4)
    #                 A4[0, 2] = np.sin(theta4) * np.sin(alpha4)
    #                 A4[0, 3] = a4 * np.cos(theta4)
    #                 A4[1, 0] = np.sin(theta4)
    #                 A4[1, 1] = np.cos(theta4) * np.cos(alpha4)
    #                 A4[1, 2] = -np.cos(theta4) * np.sin(alpha4)
    #                 A4[1, 3] = a4 * np.sin(theta4)
    #                 A4[2, 0] = 0
    #                 A4[2, 1] = np.sin(alpha4)
    #                 A4[2, 2] = np.cos(alpha4)
    #                 A4[2, 3] = d4
    #                 A4[3, 0] = 0
    #                 A4[3, 1] = 0
    #                 A4[3, 2] = 0
    #                 A4[3, 3] = 1
    #
    #                 R = np.matmul(A1, A2)
    #                 R = np.matmul(A2, A3)
    #                 R = np.matmul(R, A4)
    #                 R = R[:3, :3]
    #
    #                 base = np.r_[0, -1, 0]
    #                 point_cloud_wrist.append(np.matmul(base, R))
    # for theta1 in np.linspace(np.deg2rad(-60), np.deg2rad(180), 20):
    #     for theta2 in np.linspace(np.deg2rad(-40), np.deg2rad(120), 20):

def butter_lowpass(cutOff, fs, order=5):
    nyq = 0.5 * fs
    normalCutoff = cutOff / nyq
    b, a = butter(order, normalCutoff, btype='low', analog = True)
    return b, a

def butter_lowpass_filter(data, cutOff, fs, order=4):
    b, a = butter_lowpass(cutOff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_highpass(cutOff, fs, order=5):
    nyq = 0.5 * fs
    normalCutoff = cutOff / nyq
    b, a = butter(order, normalCutoff, btype='high', analog = True)
    return b, a

def butter_highpass_filter(data, cutOff, fs, order=4):
    b, a = butter_lowpass(cutOff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def lp_filter(data):
    from scipy.fftpack import fft
    # Number of sample points
    N = len(data)
    # sample spacing
    T = 1.0 / RATE
    x = np.linspace(0.0, N * T, N)
    y = data[:, 0]
    yf = fft(y)
    from scipy.signal.windows import blackman
    w = blackman(N)
    ywf = fft(y * w)
    xf = np.linspace(0.0, 1.0 / (2.0 * T), N / 2)

    # plt.semilogy(xf[1:N // 2], 2.0 / N * np.abs(yf[1:N // 2]), '-g')
    # plt.semilogy(xf[1:N // 2], 2.0 / N * np.abs(ywf[1:N // 2]), '-o')
    # plt.legend(['FFT', 'FFT w. window'])
    # plt.grid()
    # plt.plot()
    # plt.show()


    cutOff = 4000 / 6283185300 #cutoff frequency in rad/s
    # 10^9 in rad/s
    fs = 6283185300 #sampling frequency in rad/s
    order = 10 #order of filter

    #print sticker_data.ps1_dxdt2
    lped = butter_lowpass_filter(data, cutOff, fs, order)
    # fig = plt.figure()
    # plt.plot(data)
    # fig2 = plt.figure()
    # plt.plot(data-y)
    # plt.plot(y)

    return lped

def hp_filter(data):
    cutOff = 10000000 / 6283185300 #cutoff frequency in rad/s
    # 10^9 in rad/s
    fs = 6283185300 #sampling frequency in rad/s
    order = 10 #order of filter

    #print sticker_data.ps1_dxdt2
    hped = butter_highpass_filter(data, cutOff, fs, order)

    return hped

def clean_acc(acc_array):
    fig = plt.figure()
    plt.plot(acc_array)

    new_acc_array = np.array()
    # for i in range(0, len(acc_array), 3):
        # np.average()


def skin_dead_reckon():
    # load data
    f = open('Data/Triangle/Accelerometer.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time = [float(i) for i in T[0]]
    aX = [float(i) for i in T[1]]
    aY = [float(i) for i in T[2]]
    aZ = [float(i) for i in T[3]]
    f = open('Data/Triangle/Gyroscope.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    wX = [float(i) for i in T[1]]
    wY = [float(i) for i in T[2]]
    wZ = [float(i) for i in T[3]]
    f = open('Data/Triangle/MagneticField.csv', 'r')
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

    # Baseline
    q1, pos1, vel1 = skin.imus.analytical(R_initialOrientation=skin.quat.convert( \
        skin.imus.kalman(rate=RATE, acc=acc_array, omega=omega_array, mag=mag_array)[0], to='rotmat'), \
        omega=omega_array, initialPosition=np.array([0, 0, 0]), \
        accMeasured=acc_array, rate=1000000000)

    # low pass the accelerometer value
    acc_array = lp_filter(acc_array)
    # high pass the gyro value
    omega_array = hp_filter(omega_array)

    # calculate the location of gravity at time 0
    initial_orientation = skin.quat.convert(
        skin.imus.kalman(rate=RATE, acc=acc_array, omega=omega_array, mag=mag_array)[0], to='rotmat')
    g0 = np.linalg.inv(initial_orientation).dot(GRAVITY)

    # calculate the quaternion that rotates the first accelerometer value into the
    # frame of gravity
    q0 = skin.vector.q_shortest_rotation(acc_array[0], g0)

    # obtain the initial quaternion (aka the quaternion of the initial orientation)
    q_initial = skin.rotmat.convert(initial_orientation, to='quat')

    # rotate the initial orientation to the orientation to the direction the first
    # accelerometer value is facing
    q_ref = skin.quat.q_mult(q_initial, q0)

    # calculate quaternions step by step using the omega this is essentially the
    # "integration" of the gyroscope values this is because when you multiply
    # two quaternions together it is essentially the same as "rotating" from one
    # orientation to a second orientation we do this in the frame of the smartwatch
    # as the values read from the IMU is always in the reference frame of the smartwatch
    # q = skin.quat.calc_quat(omega_array, q_ref, RATE, 'bf')
    q = skin.imus.kalman(RATE, acc_array, omega_array, mag_array) #skin.quat.calc_quat(omega_array, q_ref, RATE, 'bf')
    q = skin.quat.q_mult(q_ref, q)

    # At every timestep define where the gravity vector is
    bf_grav = skin.vector.rotate_vector(GRAVITY, skin.quat.q_inv(q))
    # Remove gravity from the accelerometer readings
    acc_no_gravity = acc_array - bf_grav
    # Now re-rotate the cleaned accelerometer readings so they line up with the
    # location of the sensors
    acc_bf = skin.vector.rotate_vector(acc_no_gravity, q)

    # Move quaternions to the first orientation
    q = skin.quat.q_mult(q, skin.quat.q_inv(q[0]))

    # clean quaternions

    # Position and Velocity through integration, assuming 0-velocity at t=0
    vel = np.nan*np.ones_like(acc_bf)
    pos = np.nan*np.ones_like(acc_bf)

    for ii in range(acc_bf.shape[1]):
        vel[:,ii] = cumtrapz(acc_bf[:,ii], dx=1./RATE, initial=0)
        pos[:,ii] = cumtrapz(vel[:,ii],        dx=1./RATE, initial=np.array([0, 0, 0])[ii])



    # plot
    # print(q1)
    # print(pos1)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pos1[:, 0], pos1[:, 1], pos1[:, 2], c='r', marker='o')
    ax.scatter(pos[:, 0], pos[:, 1], pos1[:, 2], c='g', marker='o')
    plt.show()
