%% Housekeeping
 
addpath('ximu_matlab_library');	% include x-IMU MATLAB library
addpath('quaternion_library');	% include quatenrion library
close all;                     	% close all figures
clear;                         	% clear all variables
clc;                          	% clear the command terminal
 
% Import data

xIMUdata = xIMUdataClass('LoggedData/LoggedData');

samplePeriod = 1/256;

acc = csvread('C:\Users\pk\Code\SKIN_TEST\Data\3DShape1\Accelerometer.csv');
time = acc(:,1);
acc = acc(:,2:4);
gyr = csvread('C:\Users\pk\Code\SKIN_TEST\Data\3DShape1\Gyroscope.csv');
gyr = gyr(:,2:4);
% 
time_diff = []
for i=1:length(time)-1
    time_diff = [time_diff ; time(i+1) - time(i)];
end
samplePeriod = mean(time_diff) * 10^-9

% Process data through AHRS algorithm (calcualte orientation)
% See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/

R = zeros(3,3,length(gyr));     % rotation matrix describing sensor relative to Earth

ahrs = MahonyAHRS('SamplePeriod', samplePeriod, 'Kp', 1);

for i = 1:length(gyr)
    if (i ==1)
        sp = samplePeriod;
    else
        sp = time(i) - time(i-1);
        sp = sp * 10^-9;
    end
    ahrs.UpdateIMU(gyr(i,:), acc(i,:), sp);	% gyroscope units must be radians
    R(:,:,i) = quatern2rotMat(ahrs.Quaternion)';    % transpose because ahrs provides Earth relative to sensor
end

% Calculate 'tilt-compensated' accelerometer

tcAcc = zeros(size(acc));  % accelerometer in Earth frame
gravity_vector = [mean(acc(1:3, 1)) mean(acc(1:3, 2)) mean(acc(1:3, 3))];
for i = 1:length(acc)
    cur_gravity = R(:,:,i) * gravity_vector';
    tcAcc(i,:) = acc(i,:) - cur_gravity';
%     tcAcc(i,:) = R(:,:,i) * acc(i,:)';
end

% Calculate linear acceleration in Earth frame (subtracting gravity)

linAcc = tcAcc - [zeros(length(tcAcc), 1), zeros(length(tcAcc), 1), ones(length(tcAcc), 1)];
% linAcc = linAcc * 9.81;     % convert from 'g' to m/s/s

% Calculate linear velocity (integrate acceleartion)

linVel = zeros(size(linAcc));

for i = 2:length(linAcc)
    if i==0
        sp = 1/256;
    else
        sp = time(i) - time(i-1);
        sp = sp * 10^-9;
    end
    linVel(i,:) = linVel(i-1,:) + linAcc(i,:) * sp; %samplePeriod;
end

% High-pass filter linear velocity to remove drift

order = 1;
filtCutOff = 0.1;
[b, a] = butter(order, (2*filtCutOff)/(1/samplePeriod), 'high');
linVelHP = filtfilt(b, a, linVel);

% Calculate linear position (integrate velocity)

linPos = zeros(size(linVelHP));

for i = 2:length(linVelHP)
    if i==0
        sp = 1/256;
    else
        sp = time(i) - time(i-1);
        sp = sp * 10^-9;
    end
    linPos(i,:) = linPos(i-1,:) + linVelHP(i,:) * sp; % samplePeriod;
end

% High-pass filter linear position to remove drift

order = 1;
filtCutOff = 0.1;
[b, a] = butter(order, (2*filtCutOff)/(1/samplePeriod), 'high');
linPosHP = filtfilt(b, a, linPos);

scatter3(linPosHP(:,1), linPosHP(:,2), linPosHP(:,3))
%% End of script