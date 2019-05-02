import csv

import skinematics as skin
import numpy as np
from skinematics.sensors.manual import MyOwnSensor
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.signal import butter, lfilter
from scipy.signal import medfilt
from scipy.signal import freqs

RATE = 10**9
GRAVITY = np.r_[0, 0, 9.8]


def analytical(R_initialOrientation=np.eye(3),
               omega=np.zeros((5, 3)),
               initialPosition=np.zeros(3),
               accMeasured=np.column_stack((np.zeros((5, 2)), 9.81 * np.ones(5))),
               rate=225, time_array=np.zeros(6000)):
    ''' Reconstruct position and orientation with an analytical solution,
    from angular velocity and linear acceleration.
    Assumes a start in a stationary position. No compensation for drift.

    Parameters
    ----------
    R_initialOrientation: ndarray(3,3)
        Rotation matrix describing the initial orientation of the sensor,
        except a mis-orienation with respect to gravity
    omega : ndarray(N,3)
        Angular velocity, in [rad/s]
    initialPosition : ndarray(3,)
        initial Position, in [m]
    accMeasured : ndarray(N,3)
        Linear acceleration, in [m/s^2]
    rate : float
        sampling rate, in [Hz]

    Returns
    -------
    q : ndarray(N,3)
        Orientation, expressed as a quaternion vector
    pos : ndarray(N,3)
        Position in space [m]
    vel : ndarray(N,3)
        Velocity in space [m/s]

    Example
    -------

    >>> q1, pos1 = analytical(R_initialOrientation, omega, initialPosition, acc, rate)

    '''

    # Transform recordings to angVel/acceleration in space --------------

    # Orientation of \vec{g} with the sensor in the "R_initialOrientation"
    g = constants.g
    g0 = np.linalg.inv(R_initialOrientation).dot(np.r_[0, 0, g])

    # for the remaining deviation, assume the shortest rotation to there
    q0 = vector.q_shortest_rotation(accMeasured[0], g0)

    q_initial = rotmat.convert(R_initialOrientation, to='quat')

    # combine the two, to form a reference orientation. Note that the sequence
    # is very important!
    q_ref = quat.q_mult(q_initial, q0)

    # Calculate orientation q by "integrating" omega -----------------
    q = quat.calc_quat(omega, q_ref, rate, 'bf', time_array)

    # Acceleration, velocity, and position ----------------------------
    # From q and the measured acceleration, get the \frac{d^2x}{dt^2}
    g_v = np.r_[0, 0, g]
    accReSensor = accMeasured - vector.rotate_vector(g_v, quat.q_inv(q))
    accReSpace = vector.rotate_vector(accReSensor, q)

    # Make the first position the reference position
    q = quat.q_mult(q, quat.q_inv(q[0]))

    # compensate for drift
    # drift = np.mean(accReSpace, 0)
    # accReSpace -= drift*0.7

    # Position and Velocity through integration, assuming 0-velocity at t=0
    vel = np.nan * np.ones_like(accReSpace)
    pos = np.nan * np.ones_like(accReSpace)

    for i in range(len(accReSpace)):
        if (i < len(accReSpace) - 1):
            dt = (time_array[i + 1] - time_array[i]) * 10 ** (-9)
        for ii in range(accReSpace.shape[1]):
            if (i == 0):
                vel[i, :] = accReSpace[i, :] * dt
                pos[i, :] = vel[i, :] * dt
            else:
                vel[i, :] = vel[i - 1, :] + accReSpace[i, :] * dt
                pos[i, :] = pos[i - 1, :] + vel[i, :] * dt
            print(vel[i, 0])

    # for ii in range(accReSpace.shape[1]):
    #     vel[:,ii] = cumtrapz(accReSpace[:,ii], dx=1./rate, initial=0)
    #     pos[:,ii] = cumtrapz(vel[:,ii],        dx=1./rate, initial=initialPosition[ii])

    return (q, pos, vel)


def calc_quat(omega, q0, rate, CStype, time_array):
    '''
    Take an angular velocity (in rad/s), and convert it into the
    corresponding orientation quaternion.

    Parameters
    ----------
    omega : array, shape (3,) or (N,3)
        angular velocity [rad/s].
    q0 : array (3,)
        vector-part of quaternion (!!)
    rate : float
        sampling rate (in [Hz])
    CStype:  string
        coordinate_system, space-fixed ("sf") or body_fixed ("bf")

    Returns
    -------
    quats : array, shape (N,4)
        unit quaternion vectors.

    Notes
    -----
    1) The output has the same length as the input. As a consequence, the last velocity vector is ignored.
    2) For angular velocity with respect to space ("sf"), the orientation is given by

      .. math::
          q(t) = \\Delta q(t_n) \\circ \\Delta q(t_{n-1}) \\circ ... \\circ \\Delta q(t_2) \\circ \\Delta q(t_1) \\circ q(t_0)

      .. math::
        \\Delta \\vec{q_i} = \\vec{n(t)}\\sin (\\frac{\\Delta \\phi (t_i)}{2}) = \\frac{\\vec \\omega (t_i)}{\\left| {\\vec \\omega (t_i)} \\right|}\\sin \\left( \\frac{\\left| {\\vec \\omega ({t_i})} \\right|\\Delta t}{2} \\right)

    3) For angular velocity with respect to the body ("bf"), the sequence of quaternions is inverted.

    4) Take care that you choose a high enough sampling rate!

    Examples
    --------
    >>> v0 = np.r_[0., 0., 100.] * np.pi/180.
    >>> omega = np.tile(v0, (1000,1))
    >>> rate = 100
    >>> out = quat.calc_quat(omega, [0., 0., 0.], rate, 'sf')
    array([[ 1.        ,  0.        ,  0.        ,  0.        ],
       [ 0.99996192,  0.        ,  0.        ,  0.00872654],
       [ 0.9998477 ,  0.        ,  0.        ,  0.01745241],
       ...,
       [-0.74895572,  0.        ,  0.        ,  0.66262005],
       [-0.75470958,  0.        ,  0.        ,  0.65605903],
       [-0.76040597,  0.        ,  0.        ,  0.64944805]])
    '''

    omega_05 = np.atleast_2d(omega).copy()

    # The following is (approximately) the quaternion-equivalent of the trapezoidal integration (cumtrapz)
    if omega_05.shape[1] > 1:
        omega_05[:-1] = 0.5 * (omega_05[:-1] + omega_05[1:])

    omega_t = np.sqrt(np.sum(omega_05 ** 2, 1))
    omega_nonZero = omega_t > 0

    # initialize the quaternion
    q_delta = np.zeros(omega_05.shape)
    q_pos = np.zeros((len(omega_05), 4))
    q_pos[0, :] = unit_q(q0)

    # magnitude of position steps

    for ii in range(len(omega_05) - 1):
        if (ii > 0):
            rate = float(10 ** 9) / (time_array[ii] - time_array[ii - 1])
        else:
            rate = 225
        dq_total = np.sin(omega_t[omega_nonZero] / (2. * rate))

        q_delta[omega_nonZero, :] = omega_05[omega_nonZero, :] * np.tile(dq_total / omega_t[omega_nonZero], (3, 1)).T

        q1 = unit_q(q_delta[ii, :])
        q2 = q_pos[ii, :]
        if CStype == 'sf':
            qm = q_mult(q1, q2)
        elif CStype == 'bf':
            qm = q_mult(q2, q1)
        else:
            print('I don''t know this type of coordinate system!')
        q_pos[ii + 1, :] = qm

    return q_pos

***************************************************************************************************************
def calc_rate(df):
    length=len(df['time'])
    end=df['time'][length-1]
    begin=df['time'][0]
    rate=1/((end-begin)/(10**9*length))
    return rate
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
                    l_elbow * np.cos(theta2) * np.sin(theta1),
                    l_elbow * np.sin(theta2),
                    l_elbow * -np.cos(theta1) * np.cos(theta2)
                    ])
            # point = point * l_elbow

            point_cloud_elbow.append(point)
    point_cloud_elbow = np.array(point_cloud_elbow)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(point_cloud_elbow[:, 0], point_cloud_elbow[:, 1], point_cloud_elbow[:, 2], c='g', marker='*')
    # plt.show()

    A1 = np.zeros((4, 4))
    A2 = np.zeros((4, 4))

    point_cloud_elbow2 = []
    # np.linspace(np.deg2rad(-180), np.deg2rad(180)):
    for t1 in np.linspace(np.deg2rad(-60), np.deg2rad(180), 20):
        for t2 in np.linspace(np.deg2rad(-40), np.deg2rad(120), 20):
            theta1 = t2
            theta2 = t1
            a1 = 0
            alpha1 = -np.pi/2
            # alpha1 = 0
            d1 = l_elbow
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
            alpha2 = 0
            d2 = l_elbow
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
            base = np.r_[0,-1,0]

            # point = np.matmul(base,A)


            # if np.matmul(base,R)[2] > 0:
            point_cloud_elbow2.append(np.matmul(R, l_elbow * base))
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


    cutOff = 428.654824 #cutoff frequency in rad/s
    # 10^9 in rad/s
    fs = 628.31853 #sampling frequency in rad/s
    order = 1 #order of filter

    #print sticker_data.ps1_dxdt2
    lped = butter_lowpass_filter(data, cutOff, fs, order)
    # fig = plt.figure()
    # plt.plot(data)
    # fig2 = plt.figure()
    # plt.plot(data-y)
    # plt.plot(y)

    return lped

def hp_filter(data):
    cutOff = 100.654824 #cutoff frequency in rad/s
    # 10^9 in rad/s
    fs = 628.31853 #sampling frequency in rad/s
    order = 1 #order of filter

    #print sticker_data.ps1_dxdt2
    hped = butter_highpass_filter(data, cutOff, fs, order)

    return hped

def clean_acc(acc_array):
    fig = plt.figure()
    plt.plot(acc_array)

    new_acc_array = np.array()
    # for i in range(0, len(acc_array), 3):
        # np.average()

def dead_reckon(aX, aY, aZ, wX, wY, wZ, time):
    X = [0]
    Y = [0]
    Z = [0]
    R = np.eye(3)
    vx = 0
    vy = 0
    vz = 0
    for i in range(len(time) - 2):
        dt = (time[i + 2] - time[i + 1]) * (1/RATE)
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
        dtt = (time[i + 1] - time[i]) * (1/RATE)

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
    return X, Y, Z


def running_mean(l, N):
    sum = 0
    result = list(0 for x in l)

    for i in range(0, N):
        sum = sum + l[i]
        result[i] = sum / (i + 1)

    for i in range(N, len(l)):
        sum = sum - l[i - N] + l[i]
        result[i] = sum / N

    return result

def old_data_triangle():
    f = open('Data/Shape3/Accelerometer.csv', 'r')

    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time_acc = [float(i) for i in T[0]]
    aX = [float(i) for i in T[1]]
    aY = [float(i) for i in T[2]]
    aZ = [float(i) for i in T[3]]
    f = open('Data/Shape3/Gyroscope.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time_omega = [float(i) for i in T[0]]
    wX = [float(i) for i in T[1]]
    wY = [float(i) for i in T[2]]
    wZ = [float(i) for i in T[3]]
    f = open('Data/Shape3/MagneticField.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time_mag = [float(i) for i in T[0]]
    mX = [float(i) for i in T[1]]
    mY = [float(i) for i in T[2]]
    mZ = [float(i) for i in T[3]]

    return time_acc, aX, aY, aZ, wX, wY, wZ, mX, mY, mZ

def old_data_rectangle():
    f = open('Data/Rectangle/Accelerometer.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time_acc = [float(i) for i in T[0]]
    aX = [float(i) for i in T[1]]
    aY = [float(i) for i in T[2]]
    aZ = [float(i) for i in T[3]]
    f = open('Data/Rectangle/Gyroscope.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time_omega = [float(i) for i in T[0]]
    wX = [float(i) for i in T[1]]
    wY = [float(i) for i in T[2]]
    wZ = [float(i) for i in T[3]]
    f = open('Data/Rectangle/MagneticField.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time_mag = [float(i) for i in T[0]]
    mX = [float(i) for i in T[1]]
    mY = [float(i) for i in T[2]]
    mZ = [float(i) for i in T[3]]

    return time_acc, aX, aY, aZ, wX, wY, wZ, mX, mY, mZ

def skin_dead_reckon(df, old_data, type):

    time_acc = list(df['time'])

    aX = list(df['aX'])
    aY = list(df['aY'])
    aZ = list(df['aZ'])

    wX = list(df['wX'])
    wY = list(df['wY'])
    wZ = list(df['wZ'])

    mX = list(df['mX'])
    mY = list(df['mY'])
    mZ = list(df['mZ'])

    if old_data and type == 'triangle':
        time_acc, aX, aY, aZ, wY, wX, wZ, mX, mY, mZ = old_data_triangle()
    if old_data and type == 'rectangle':
        time_acc, aX, aY, aZ, wY, wX, wZ, mX, mY, mZ = old_data_rectangle()




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


    # low pass the accelerometer value
    # acc_array = lp_filter(acc_array)
    # high pass the gyro value
    # omega_array = hp_filter(omega_array)
    q1, pos1, vel1 = skin.imus.analytical(R_initialOrientation=np.eye(3), \
        omega=omega_array, initialPosition=np.array([0, 0, 0]), \
        accMeasured=acc_array, rate=225, time_array=time_acc)
    aX=running_mean(aX,5)
    aY=running_mean(aY,5)
    aZ=running_mean(aZ,5)
    wX=running_mean(wX,5)
    wY=running_mean(wY,5)
    wZ=running_mean(wZ,5)
    # calculate the location of gravity at time 0
    initial_orientation = np.eye(3) #skin.quat.convert(
        # skin.imus.kalman(rate=10**3, acc=acc_array, omega=omega_array, mag=mag_array)[0], to='rotmat')
    print("Initial orientation: ", initial_orientation)
    g0 = np.linalg.inv(initial_orientation).dot(GRAVITY)
    print("G0", g0)

    # calculate the quaternion that rotates the first accelerometer value into the
    # frame of gravity
    q0 = skin.vector.q_shortest_rotation(acc_array[0], g0)
    print(q0)

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
    q = []
    print("_____________",acc_array[0])
    for i in range(len(time_acc)-1):
        time_diff = time_acc[i+1] - time_acc[i]
        q_temp = skin.imus.kalman(time_diff/(10**9), acc_array[i:i+1], omega_array[i:i+1], mag_array[i:i+1])
        q.append(skin.quat.q_mult(q_temp, q_ref))
        # q_ref = q_temp
        # q.append(q_temp)
    # print(skin.imus.kalman(RATE, acc_array, omega_array, mag_array).shape) #skin.quat.calc_quat(omega_array, q_ref, RATE, 'bf')
    q = np.asarray(q)
    q = np.squeeze(q)
    # print(q)
    # print(q.shape)
    # q = skin.quat.q_mult(q_ref, q)

    # At every timestep define where the gravity vector is
    bf_grav = skin.vector.rotate_vector(GRAVITY, skin.quat.q_inv(q))
    # Remove gravity from the accelerometer readings
    # acc_no_gravity = acc_array - bf_grav
    acc_no_gravity = acc_array[:-1] - bf_grav
    # Now re-rotate the cleaned accelerometer readings so they line up with the
    # location of the sensors
    acc_bf = skin.vector.rotate_vector(acc_no_gravity, q)
    aX_new = acc_bf[:,0]
    aY_new = acc_bf[:,1]
    aZ_new = acc_bf[:,2]

    # Move quaternions to the first orientation
    q = skin.quat.q_mult(q, skin.quat.q_inv(q[0]))

    # Position and Velocity through integration, assuming 0-velocity at t=0
    vx = 0
    vy = 0
    vz = 0
    X = [0]
    Y = [0]
    Z = [0]
    gravity_vector=np.array([[np.average(aX[0:2])],[np.average(aY[0:2])],[np.average(aZ[0:2])]])


    for i in range(len(time_acc) - 2):
        dt = (time_acc[i + 2] - time_acc[i + 1]) * 10**(-9)
        dtt = (time_acc[i + 1] - time_acc[i]) * 10**(-9)
        cur_gravity = np.matmul(np.linalg.inv(skin.quat.convert(q[i], to='rotmat')), gravity_vector)
        aX_cur=aX[i] - cur_gravity[0][0]
        aY_cur=aY[i] - cur_gravity[1][0]
        aZ_cur=aZ[i] - cur_gravity[2][0]
        if (np.sqrt(aX_cur**2+aY_cur**2+aZ_cur**2)>10):
            aX_cur*=0.9
            aY_cur*=0.9
            aZ_cur*=0.9
        dv = np.array([[dtt * aX_cur],
                       [dtt * aY_cur],
                       [dtt * aZ_cur]])
        dv_world = np.matmul(skin.quat.convert(q[i], to='rotmat'), dv)

        vx += dv_world[0][0]
        vy += dv_world[1][0]
        vz += dv_world[2][0]
        # print(vy,aY_cur)
        # if (np.sqrt(vx**2+vy**2+vz**2)>0.5):
        #     vx*=0.9
        #     vy*=0.9
        #     vz*=0.9

        X.append(X[-1] + vx * dt)
        Y.append(Y[-1] + vy * dt)
        Z.append(Z[-1] + vz * dt)

    # for i in range(len(time_acc) - 2):
    #     dt = (time_acc[i + 2] - time_acc[i + 1]) * 10**(-9)
    #     dtt = (time_acc[i + 1] - time_acc[i]) * 10**(-9)
    #
    #     dv = np.array([[dtt * aX_new[i]],
    #                    [dtt * aY_new[i]],
    #                    [dtt * (aZ_new[i]-aZ_new[0])]])
    #     # dv_world = np.matmul(skin.quat.convert(q[i], to='rotmat'), dv)
    #
    #     vx += dv[0][0]
    #     vy += dv[1][0]
    #     vz += dv[2][0]
    #     X.append(X[-1] + vx * dt * 10)
    #     Y.append(Y[-1] + vy * dt * 10)
    #     Z.append(Z[-1] + vz * dt * 10)

    # Position and Velocity through integration, assuming 0-velocity at t=0
    vel = np.nan * np.ones_like(acc_bf)
    pos = np.nan * np.ones_like(acc_bf)
    initialPosition = np.zeros(3)

    # for i in range(len(time_acc)):
    #     time_acc[i] = time_acc[i] / 10**9
    for ii in range(acc_array.shape[1]):
        vel[:, ii] = cumtrapz(acc_bf[:, ii], time_acc[:-1], initial=0)
        pos[:, ii] = cumtrapz(vel[:, ii], time_acc[:-1], initial=initialPosition[ii])

    for j in range(3):
        for i in range(pos.shape[0]):
            pos[i][j] = pos[i][j]/10**17


    X_base, Y_base, Z_base = dead_reckon(aX, aY, aZ, wX, wY, wZ, time_acc)
    # plot
    # print(q1)
    # print(pos1)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # t = np.linspace(0, 2 * np.pi, len(pos1[:,0]))
    # ax.scatter(X_base, Y_base, Z_base, c='r', marker='o')
    # ax.scatter(X, Y, Z, c='b', marker='o')
    # ax.scatter(pos[:,0], pos[:,1], pos[:,2], c='g', marker='*')
    ax.scatter(pos1[:, 0], pos1[:,1], pos1[:, 2], c='r', marker='*')
    # plt.plot(np.linspace(0, 1, len(aX_new)), aX_new)
    # plt.plot(np.linspace(0, 1, len(aX_new)), aY_new)
    # plt.plot(np.linspace(0, 1, len(aX_new)), aZ_new)
    # plt.plot(np.linspace(0, 1, len(wX_interp)), wX_interp)
    # plt.plot(np.linspace(0, 1, len(wX_interp)), wY_interp)
    # plt.plot(np.linspace(0, 1, len(wX_interp)), wZ_interp)
    plt.show()

    # # Baseline
    # q1, pos1, vel1 = skin.imus.analytical(R_initialOrientation=np.eye(3), \
    #     omega=omega_array, initialPosition=np.array([0, 0, 0]), \
    #     accMeasured=acc_array, rate=1)
    #
    # # low pass the accelerometer value
    # acc_array = lp_filter(acc_array)
    # # high pass the gyro value
    # omega_array = hp_filter(omega_array)
    #
    # # calculate the location of gravity at time 0
    # initial_orientation = skin.quat.convert(
    #     skin.imus.kalman(rate=RATE, acc=acc_array, omega=omega_array, mag=mag_array)[0], to='rotmat')
    # print("Initial orientation: ", initial_orientation)
    # g0 = np.linalg.inv(initial_orientation).dot(GRAVITY)
    # print("G0", g0)
    #
    # # calculate the quaternion that rotates the first accelerometer value into the
    # # frame of gravity
    # q0 = skin.vector.q_shortest_rotation(acc_array[0], g0)
    # print(q0)
    #
    # # obtain the initial quaternion (aka the quaternion of the initial orientation)
    # q_initial = skin.rotmat.convert(initial_orientation, to='quat')
    #
    # # rotate the initial orientation to the orientation to the direction the first
    # # accelerometer value is facing
    # q_ref = skin.quat.q_mult(q_initial, q0)
    #
    # # calculate quaternions step by step using the omega this is essentially the
    # # "integration" of the gyroscope values this is because when you multiply
    # # two quaternions together it is essentially the same as "rotating" from one
    # # orientation to a second orientation we do this in the frame of the smartwatch
    # # as the values read from the IMU is always in the reference frame of the smartwatch
    # # q = skin.quat.calc_quat(omega_array, q_ref, RATE, 'bf')
    # q = []
    # print("_____________",acc_array[0])
    # for i in range(len(time_acc)-1):
    #     time_diff = time_acc[i+1] - time_acc[i]
    #     q_temp = skin.imus.kalman(time_diff, acc_array[i:i+1], omega_array[i:i+1], mag_array[i:i+1])
    #     q.append(q_temp)
    # # print(skin.imus.kalman(RATE, acc_array, omega_array, mag_array).shape) #skin.quat.calc_quat(omega_array, q_ref, RATE, 'bf')
    # q = np.asarray(q)
    # q = np.squeeze(q)
    # # print(q)
    # # print(q.shape)
    # q = skin.quat.q_mult(q_ref, q)
    #
    # # At every timestep define where the gravity vector is
    # bf_grav = skin.vector.rotate_vector(GRAVITY, skin.quat.q_inv(q))
    # # Remove gravity from the accelerometer readings
    # # acc_no_gravity = acc_array - bf_grav
    # acc_no_gravity = acc_array[:-1] - bf_grav
    # # Now re-rotate the cleaned accelerometer readings so they line up with the
    # # location of the sensors
    # acc_bf = skin.vector.rotate_vector(acc_no_gravity, q)
    # aX_new = acc_bf[:,0]
    # aY_new = acc_bf[:,1]
    # aZ_new = acc_bf[:,2]
    #
    # # Move quaternions to the first orientation
    # q = skin.quat.q_mult(q, skin.quat.q_inv(q[0]))
    #
    # # Position and Velocity through integration, assuming 0-velocity at t=0
    # # vx = 0
    # # vy = 0
    # # vz = 0
    # # X = [0]
    # # Y = [0]
    # # Z = [0]
    # # for i in range(len(time_acc) - 2):
    # #     dt = (time_acc[i + 2] - time_acc[i + 1]) * (1/100)
    # #     dtt = (time_acc[i + 1] - time_acc[i]) * (1/100)
    # #
    # #     dv = np.array([[dtt * aX_new[i]],
    # #                    [dtt * aY_new[i]],
    # #                    [dtt * (aZ_new[i] - aZ_new[0])]])
    # #     dv_world = np.matmul(skin.quat.convert(q[i], to='rotmat'), dv)
    # #
    # #     vx += dv_world[0][0]
    # #     vy += dv_world[1][0]
    # #     vz += dv_world[2][0]
    # #     X.append(X[-1] + vx * dt * 10)
    # #     Y.append(Y[-1] + vy * dt * 10)
    # #     Z.append(Z[-1] + vz * dt * 10)
    # initialPosition = np.zeros(3)
    # vel = np.nan*np.ones_like(acc_bf)
    # pos = np.nan*np.ones_like(acc_bf)
    #
    # print(acc_bf.shape, len(time_acc))
    #
    # for ii in range(acc_bf.shape[1]):
    #     vel[:,ii] = cumtrapz(acc_bf[:,ii], time_acc[:-1], initial=0)
    #     pos[:,ii] = cumtrapz(vel[:,ii], time_acc[:-1], initial=initialPosition[ii])
    #
    #
    # X_base, Y_base, Z_base = dead_reckon(aX, aY, aZ, wX, wY, wZ, time_acc)
    #
    # # plot
    # # print(q1)
    # # print(pos1)
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # # ax.scatter(X_base, Y_base, Z_base, c='r', marker='o')
    # # ax.scatter(X, Y, Z, c='b', marker='o')
    # ax.scatter(pos1[:,0], pos1[:,1], pos1[:,2], c='g', marker='o')
    # # ax.scatter(pos[:,0], pos[:,1], pos[:,2], c='r', marker='o')
    #
    # plt.show()
