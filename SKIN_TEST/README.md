#IMU3D Dead Reckoning
This project does IMU MEMS sensor 3D dead reckoning.
This is a PyCharm project but feel free to install all dependencies and just run from CLI.

# Dependencies
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import constants
    import csv
    import pandas as pd
    from scipy.signal import lfilter, butter
    import madgwickahrs as mad
    from mpl_toolkits.mplot3d import axes3d

Please install matplotlib, numpy, scipy, pandas.
Madgwick is provided in the project.

#Folder Structure
* Root: Contains everything including all the Python scripts
* PICTURES: Contains our plots of all the various provided data
    * Note the format of the name is SHAPE_CALCQUAT%_MADGWICK%.png
* Data: Contains all the data
    * Self: Contains self collected data, the format of this data is different from provided data
    * All other folders contain provided data which is Accelerometer.csv, MagneticField.csv, and Gyroscope.csv per shape
* data_collection_4_final: Contains an Android Studio project used to collect data from a phone.

#Usage
Please run imu3d.py it will prompt you to enter values.

    Would you like to run the baseline model? (BL) Or the final model? (F): F
    Do you wish do use data collected? (Y) Or provided data? (N): N
    Enter the name of the shape (Self: B, Provided: Rectangle): Rectangle
    Enter the % calc_quat you would like to use: 0.5
    Enger the % tq_scale you would like to use: 0.5
    
This is an example usage case of the code.