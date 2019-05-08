import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
import csv

from scipy.signal import lfilter, butter

import madgwickahrs as mad

from skinematics.sensors.manual import MyOwnSensor


def butter_lowpass(cutOff, fs, order=5):
    nyq = 0.5 * fs
    normalCutoff = cutOff / nyq
    b, a = butter(order, normalCutoff, btype='low', analog=True)
    return b, a


def butter_lowpass_filter(data, cutOff, fs, order=4):
    b, a = butter_lowpass(cutOff, fs, order=order)
    y = lfilter(b, a, data)
    return y


def butter_highpass(cutOff, fs, order=5):
    nyq = 0.5 * fs
    normalCutoff = cutOff / nyq
    b, a = butter(order, normalCutoff, btype='high', analog=True)
    return b, a


def butter_highpass_filter(data, cutOff, fs, order=4):
    b, a = butter_lowpass(cutOff, fs, order=order)
    y = lfilter(b, a, data)
    return y


def lp_filter(data, cutOff=400, fs=628):
    # cutOff = 428.654824 #cutoff frequency in rad/s
    # 10^9 in rad/s
    # fs = 628.31853 #sampling frequency in rad/s
    order = 1  # order of filter

    # print sticker_data.ps1_dxdt2
    lped = butter_lowpass_filter(data, cutOff, fs, order)

    return lped


def hp_filter(data, cutOff=400, fs=628):
    # cutOff = 100.654824 #cutoff frequency in rad/s
    # 10^9 in rad/s
    # fs = 628.31853 #sampling frequency in rad/s
    order = 1  # order of filter

    # print sticker_data.ps1_dxdt2
    hped = butter_highpass_filter(data, cutOff, fs, order)

    return hped


def old_data_rectangle():
    f = open('Data/Triangle/Accelerometer.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time_acc = [float(i) for i in T[0]]
    aX = [float(i) for i in T[1]]
    aY = [float(i) for i in T[2]]
    aZ = [float(i) for i in T[3]]
    f = open('Data/Triangle/Gyroscope.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time_omega = [float(i) for i in T[0]]
    wX = [float(i) for i in T[1]]
    wY = [float(i) for i in T[2]]
    wZ = [float(i) for i in T[3]]
    f = open('Data/Triangle/MagneticField.csv', 'r')
    readfile = csv.reader(f)
    T = list(map(list, zip(*readfile)))
    time_mag = [float(i) for i in T[0]]
    mX = [float(i) for i in T[1]]
    mY = [float(i) for i in T[2]]
    mZ = [float(i) for i in T[3]]

    return time_acc, aX, aY, aZ, wX, wY, wZ, mX, mY, mZ


def grab_data_from_df():
    time_acc, aX, aY, aZ, wY, wX, wZ, mX, mY, mZ = old_data_rectangle()

    # determine which of the sensors has the longest array
    acc_len = len(aX)
    omega_len = len(wX)
    mag_len = len(mX)
    time_len = len(time_acc)

    max_len = max([acc_len, omega_len, mag_len, time_len])

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
    time_acc_interp = np.interp(np.linspace(0, 1, max_len), np.linspace(0, 1, len(time_acc)), time_acc)

    # create N x 3 np arrays
    acc_array = np.zeros((max_len, 3))
    omega_array = np.zeros((max_len, 3))
    mag_array = np.zeros((max_len, 3))
    time_array = np.zeros((max_len))

    #averaging filter
    N=30
    aX_temp=np.convolve(aX_interp, np.ones((N,)) / N, mode='valid')
    aY_temp=np.convolve(aY_interp, np.ones((N,)) / N, mode='valid')
    aZ_temp=np.convolve(aZ_interp, np.ones((N,)) / N, mode='valid')
    wX_temp = np.convolve(wX_interp, np.ones((N,)) / N, mode='valid')
    wY_temp = np.convolve(wY_interp, np.ones((N,)) / N, mode='valid')
    wZ_temp = np.convolve(wZ_interp, np.ones((N,)) / N, mode='valid')


    aX_interp[int(N/2)-1:int(-N/2)]=aX_temp
    aY_interp[int(N/2)-1:int(-N/2)]=aY_temp
    aZ_interp[int(N/2)-1:int(-N/2)]=aZ_temp
    wX_interp[int(N / 2) - 1:int(-N / 2)] = wX_temp
    wY_interp[int(N / 2) - 1:int(-N / 2)] = wY_temp
    wZ_interp[int(N / 2) - 1:int(-N / 2)] = wZ_temp


    for i in range(max_len):
        acc_array[i, 0] = aX_interp[i]
        acc_array[i, 1] = aY_interp[i]
        acc_array[i, 2] = aZ_interp[i]
    for i in range(max_len):
        omega_array[i, 0] = wX_interp[i]
        omega_array[i, 1] = wY_interp[i]
        omega_array[i, 2] = wZ_interp[i]
    for i in range(max_len):
        mag_array[i, 0] = mX_interp[i]
        mag_array[i, 1] = mY_interp[i]
        mag_array[i, 2] = mZ_interp[i]
    for i in range(max_len):
        time_array[i] = time_acc_interp[i]

    #gyro drift
    mean=np.mean(omega_array,0)
    omega_array-=0.7*mean

    return acc_array, omega_array, mag_array, time_array


def normalize(n):
    # normalize
    from numpy.linalg import norm

    if np.array(n).ndim == 1:
        vectorFlag = True
    else:
        vectorFlag = False

    n = np.double(np.atleast_2d(n))
    length = norm(n, axis=1)
    n[length != 0] = (n[length != 0].T / length[length != 0]).T
    if vectorFlag:
        n = n.ravel()
    return n


def plot_trajectory3d(pos_analytical, pos_kalman):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pos_analytical[:, 0], pos_analytical[:, 1], pos_analytical[:, 2], c='g', marker='*')

    plt.show()


def plot_trajectory2d(pos_analytical, pos_kalman):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pos_analytical[:, 0], pos_analytical[:, 1], c='r', marker='*')

    plt.show()


def quat_2_vec(v1, v2):
    # calculate the perp vec to v1 v2
    # quaternion lives here
    n = normalize(np.cross(v1, v2))

    # calculate angle between the two vectors v1 and v2
    v1 = np.array(v1)
    v2 = np.array(v2)

    n1 = normalize(v1)
    n2 = normalize(v2)
    if v2.ndim == 1:
        angle = np.arccos(n1.dot(n2))

    # calculate the quaternion
    q = (n.T * np.sin(angle / 2.)).T

    return q


def convert_rotmat_to_quat(rotmat):
    # assume we input a 3x3 matrix
    # http: // www.euclideanspace.com / maths / geometry / rotations / conversions / matrixToQuaternion /
    q = np.zeros(4)

    q[0] = np.sqrt(1.0 + rotmat[0, 0] + rotmat[1, 1] + rotmat[2, 2]) / 2.0
    q[1] = (rotmat[2, 1] - rotmat[1, 2]) / (q[0] * 4.0)
    q[2] = (rotmat[0, 2] - rotmat[2, 0]) / (q[0] * 4.0)
    q[3] = (rotmat[1, 0] - rotmat[0, 1]) / (q[0] * 4.0)

    return q
def convert_quat_to_rotmat(quat):
    q = unit_q(quat).T

    R = np.zeros((9, q.shape[1]))
    R[0] = q[0] ** 2 + q[1] ** 2 - q[2] ** 2 - q[3] ** 2
    R[1] = 2 * (q[1] * q[2] - q[0] * q[3])
    R[2] = 2 * (q[1] * q[3] + q[0] * q[2])
    R[3] = 2 * (q[1] * q[2] + q[0] * q[3])
    R[4] = q[0] ** 2 - q[1] ** 2 + q[2] ** 2 - q[3] ** 2
    R[5] = 2 * (q[2] * q[3] - q[0] * q[1])
    R[6] = 2 * (q[1] * q[3] - q[0] * q[2])
    R[7] = 2 * (q[2] * q[3] + q[0] * q[1])
    R[8] = q[0] ** 2 - q[1] ** 2 - q[2] ** 2 + q[3] ** 2

    if R.shape[1] == 1:
        return np.reshape(R, (3, 3))
    else:
        return R.T

def multiply_quaternions(p, q):


    p = unit_q(p).T
    q = unit_q(q).T

    if np.prod(np.shape(p)) > np.prod(np.shape(q)):
        r = np.zeros(np.shape(p))
    else:
        r = np.zeros(np.shape(q))

    r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3]
    r[1] = p[1] * q[0] + p[0] * q[1] + p[2] * q[3] - p[3] * q[2]
    r[2] = p[2] * q[0] + p[0] * q[2] + p[3] * q[1] - p[1] * q[3]
    r[3] = p[3] * q[0] + p[0] * q[3] + p[1] * q[2] - p[2] * q[1]

    r = r.T
    return r
def unit_q(inData):

    inData = np.atleast_2d(inData)
    (m, n) = inData.shape

    if n == 3:
        qLength = 1 - np.sum(inData ** 2, 1)
        numLimit = 1e-12
        # Check for numerical problems
        if np.min(qLength) < -numLimit:
            raise ValueError('Quaternion is too long!')
        else:
            # Correct for numerical problems
            qLength[qLength < 0] = 0
        outData = np.hstack((np.c_[np.sqrt(qLength)], inData))

    else:
        outData = inData

    return outData


def calc_quat(omega, q0, time_array):

    omega_05 = np.atleast_2d(omega).copy()
    if omega_05.shape[1] > 1:
        omega_05[:-1] = 0.5 * (omega_05[:-1] + omega_05[1:])

    omega_t = np.sqrt(np.sum(omega_05 ** 2, 1))
    omega_nonZero = omega_t > 0

    # initialize the quaternion
    q_delta = np.zeros(omega_05.shape)
    q_pos = np.zeros((len(omega_05), 4))
    q_pos[0, :] = unit_q(q0)


    for ii in range(len(omega_05) - 1):
        if (ii > 0):
            rate = float(10 ** 9) / (time_array[ii] - time_array[ii - 1])
        else:
            rate = 225
        dq_total = np.sin(omega_t[omega_nonZero] / (2. * rate))

        q_delta[omega_nonZero, :] = omega_05[omega_nonZero, :] * np.tile(dq_total / omega_t[omega_nonZero], (3, 1)).T

        q1 = unit_q(q_delta[ii, :])
        q2 = q_pos[ii, :]
        qm = multiply_quaternions(q2, q1)

        q_pos[ii + 1, :] = qm

    return q_pos


def angle(v1, v2):
    # make sure lists are handled correctly
    v1 = np.array(v1)
    v2 = np.array(v2)

    if v1.ndim < v2.ndim:
        v1, v2 = v2, v1
    n1 = normalize(v1)
    n2 = normalize(v2)
    if v2.ndim == 1:
        angle = np.arccos(n1.dot(n2))
    else:
        angle = np.arccos(list(map(np.dot, n1, n2)))
    return angle

def q_shortest_rotation(v1, v2):

    # calculate the direction
    n = normalize(np.cross(v1, v2))

    # make sure vectors are handled correctly
    n = np.atleast_2d(n)

    # handle 0-quaternions
    nanindex = np.isnan(n[:, 0])
    n[nanindex, :] = 0

    # find the angle, and calculate the quaternion
    angle12 = angle(v1, v2)
    q = (n.T * np.sin(angle12 / 2.)).T

    # if you are working with vectors, only return a vector
    if q.shape[0] == 1:
        q = q.flatten()

    return q
def q_inv(q):
    q = np.atleast_2d(q)
    if q.shape[1]==3:
        return -q
    else:
        qLength = np.sum(q**2, 1)
        qConj = q * np.r_[1, -1,-1,-1]
    return (qConj.T / qLength).T


def rotate_vector(vector, q):
    vector = np.atleast_2d(vector)
    qvector = np.hstack((np.zeros((vector.shape[0], 1)), vector))
    vRotated = multiply_quaternions(q, multiply_quaternions(qvector, q_inv(q)))
    vRotated = vRotated[:, 1:]

    if min(vRotated.shape) == 1:
        vRotated = vRotated.ravel()

    return vRotated
def self():
    # setup
    R_initialOrientation = np.eye(3)
    initialPosition = np.zeros(3)
    accMeasured, omega, mag_array, time_array = grab_data_from_df()

    # low pass the accelerometer value
    # accMeasured = lp_filter(accMeasured)
    # high pass the gyro value
    # omega = hp_filter(omega)

    in_data = {
        'rate': 256,
        'acc': accMeasured,
        'omega': omega,
        'mag': mag_array
    }

    # create your own sensor (not used here)
    my_sensor = MyOwnSensor(in_data=in_data)

    # Orientation of \vec{g} with the sensor in the "R_initialOrientation"
    g = constants.g
    g0 = np.linalg.inv(R_initialOrientation).dot(np.r_[0, 0, g])

    # for the remaining deviation, assume the shortest rotation to there
    q0 = q_shortest_rotation(accMeasured[0], g0)

    q_initial = convert_rotmat_to_quat(R_initialOrientation)

    # combine the two, to form a reference orientation. Note that the sequence
    # is very important!
    q_ref = multiply_quaternions(q_initial, q0)

    # Calculate orientation q by "integrating" omega -----------------
    q = calc_quat(omega, q_ref, time_array)

    md = mad.MadgwickAHRS()
    tq = []
    for i in range(len(accMeasured)):
        if i > 0:
            rate = (time_array[i] - time_array[i - 1]) / 10 ** -9
        else:
            rate = 225
        rate = 1 / rate
        md.set_rate(sampleRate=rate)
        md.update_imu(gyroscope=omega[i, :], accelerometer=accMeasured[i, :])
        # md.update(omega[i, :], accMeasured[i, :], mag_array[i, :])
        tq.append(md.quaternion.q)
    tq = np.asarray(tq)

    q = normalize(q * 0.6 + tq * 0.4)
    # q = q * 0.5 + tq * 0.5
    # Acceleration, velocity, and positicaon ----------------------------
    # From q and the measured acceleration, get the \frac{d^2x}{dt^2}
    g_v = np.r_[0, 0, g]

    # Make the first position the reference position
    q = multiply_quaternions(q, q_inv(q[0]))

    gravity_vector = np.array(
        [np.average(accMeasured[0:2, 0]), np.average(accMeasured[0:2, 1]), np.average(accMeasured[0:2, 2])])
    accReSensor = accMeasured - rotate_vector(g_v, q_inv(q))
    for i in range(len(q)):
        cur_gravity = np.matmul(np.linalg.inv(convert_quat_to_rotmat(q[i])), gravity_vector)
        accReSensor[i] = (accMeasured[i] - cur_gravity)
    accReSpace = rotate_vector(accReSensor, q)

    # compensate for drift
    drift = np.mean(accReSpace, 0)
    accReSpace -= drift * 0.7

    # Position and Velocity through integration, assuming 0-velocity at t=0
    vel = np.nan * np.ones_like(accReSpace)
    pos = np.nan * np.ones_like(accReSpace)

    time_diff_arr = []
    for i in range(len(accReSpace)):
        if (i < len(accReSpace) - 1):
            dt = (time_array[i + 1] - time_array[i]) * 10 ** (-9)
            time_diff_arr.append(dt)
        for ii in range(accReSpace.shape[1]):
            if (i == 0):
                vel[i, :] = accReSpace[i, :] * dt
                pos[i, :] = vel[i, :] * dt
            else:
                vel[i, :] = vel[i - 1, :] + accReSpace[i, :] * dt
                pos[i, :] = pos[i - 1, :] + vel[i, :] * dt

    return (q, pos, vel)


q, pos, vel = self()
plot_trajectory3d(pos, [])
# plot_trajectory2d(pos, [])
