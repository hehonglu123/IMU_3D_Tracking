package com.pk.data_collection;

import android.content.Context;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.os.Bundle;
import android.os.Environment;
import android.support.design.widget.FloatingActionButton;
import android.support.design.widget.Snackbar;
import android.support.v7.app.AppCompatActivity;
import android.support.v7.widget.Toolbar;
import android.util.Log;
import android.view.View;
import android.view.Menu;
import android.view.MenuItem;
import android.widget.Button;
import android.widget.TextView;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

public class MainActivity extends AppCompatActivity
    implements SensorEventListener {

    ArrayList<String> accTime = new ArrayList<>();
    ArrayList<String> accArrayX = new ArrayList<>();
    ArrayList<String> accArrayY = new ArrayList<>();
    ArrayList<String> accArrayZ = new ArrayList<>();
    ArrayList<String> magTime = new ArrayList<>();
    ArrayList<String> magArrayX = new ArrayList<>();
    ArrayList<String> magArrayY = new ArrayList<>();;
    ArrayList<String> magArrayZ = new ArrayList<>();
    ArrayList<String> gyroTime = new ArrayList<>();
    ArrayList<String> gyroArrayX = new ArrayList<>();
    ArrayList<String> gyroArrayY = new ArrayList<>();
    ArrayList<String> gyroArrayZ = new ArrayList<>();
    long write_timestamp;
    float[] mag = new float[3];
    float[] gyro = new float[3];
    float[] acc = new float[3];

    SensorManager mSensorManager;
    private FileOutputStream outFD_acc;
    private FileOutputStream outFD_mag;
    private FileOutputStream outFD_gyro;

    int file_counter;
    int count;

    Button writeButton;
    TextView text;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
        Toolbar toolbar = findViewById(R.id.toolbar);
        setSupportActionBar(toolbar);

        mSensorManager = (SensorManager) this.getSystemService(Context.SENSOR_SERVICE);
        mSensorManager.registerListener(this, mSensorManager.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD), SensorManager.SENSOR_DELAY_GAME);
        mSensorManager.registerListener(this, mSensorManager.getDefaultSensor(Sensor.TYPE_GYROSCOPE), SensorManager.SENSOR_DELAY_GAME);
        mSensorManager.registerListener(this, mSensorManager.getDefaultSensor(Sensor.TYPE_ACCELEROMETER), SensorManager.SENSOR_DELAY_GAME);

        file_counter = 0;
        writeButton = findViewById(R.id.button);
        text = findViewById(R.id.textView);
        writeButton.setOnClickListener(new clickListener());

        count = 0;
    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        // Inflate the menu; this adds items to the action bar if it is present.
        getMenuInflater().inflate(R.menu.menu_main, menu);
        return true;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        // Handle action bar item clicks here. The action bar will
        // automatically handle clicks on the Home/Up button, so long
        // as you specify a parent activity in AndroidManifest.xml.
        int id = item.getItemId();

        //noinspection SimplifiableIfStatement
        if (id == R.id.action_settings) {
            return true;
        }

        return super.onOptionsItemSelected(item);
    }

    @Override
    public void onSensorChanged(SensorEvent sensorEvent) {
        count = count + 1;
//        write_timestamp = (new Date()).getTime()
//                + (sensorEvent.timestamp - System.nanoTime()) / 1000000L;
        write_timestamp = sensorEvent.timestamp;

        if(count>500) {
            text.setText("Ready");
            if (sensorEvent.sensor.getType() == Sensor.TYPE_ACCELEROMETER) {
                acc[0] = sensorEvent.values[0];
                acc[1] = sensorEvent.values[1];
                acc[2] = sensorEvent.values[2];

                accTime.add(String.valueOf(write_timestamp));
                accArrayX.add(String.valueOf(sensorEvent.values[0]));
                accArrayY.add(String.valueOf(sensorEvent.values[1]));
                accArrayZ.add(String.valueOf(sensorEvent.values[2]));

            }
            else if (sensorEvent.sensor.getType() == Sensor.TYPE_GYROSCOPE) {
                gyro[0] = sensorEvent.values[0];
                gyro[1] = sensorEvent.values[1];
                gyro[2] = sensorEvent.values[2];

                gyroTime.add(String.valueOf(write_timestamp));
                gyroArrayX.add(String.valueOf(sensorEvent.values[0]));
                gyroArrayY.add(String.valueOf(sensorEvent.values[1]));
                gyroArrayZ.add(String.valueOf(sensorEvent.values[2]));

            }
            else if (sensorEvent.sensor.getType() == Sensor.TYPE_MAGNETIC_FIELD) {
                mag[0] = sensorEvent.values[0];
                mag[1] = sensorEvent.values[1];
                mag[2] = sensorEvent.values[2];

                magTime.add(String.valueOf(write_timestamp));
                magArrayX.add(String.valueOf(sensorEvent.values[0]));
                magArrayY.add(String.valueOf(sensorEvent.values[1]));
                magArrayZ.add(String.valueOf(sensorEvent.values[2]));

            }
        }

    }

    @Override
    public void onAccuracyChanged(Sensor sensor, int i) {

    }

    private class clickListener implements View.OnClickListener{
        @Override
        public void onClick(View v) {
            switch (v.getId()){
                case R.id.button:
                    while(!getBaseContext().getFileStreamPath("ECE498_sensor_file_acc"+file_counter+".csv").exists()) {
                        Log.d("Opening file", "asdf");
                        try {
                            outFD_acc = openFileOutput("ECE498_sensor_file_acc" + file_counter + ".csv", MODE_APPEND);
                            outFD_acc.write("time,aX,aY,aZ\n".getBytes());
                            outFD_gyro = openFileOutput("ECE498_sensor_file_gyro" + file_counter + ".csv", MODE_APPEND);
                            outFD_gyro.write("time,gX,gY,gZ\n".getBytes());
                            outFD_mag = openFileOutput("ECE498_sensor_file_gyro" + file_counter + ".csv", MODE_APPEND);
                            outFD_mag.write("time,gX,gY,gZ\n".getBytes());
                        } catch (FileNotFoundException e) {
                            Log.d("file", "failed to open file");
                            e.printStackTrace();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }


                    // write all data to file
                    for(int i = 0 ; i < accArrayX.size() ; i++) {
                        String accT = accTime.get(i);
                        String accX = accArrayX.get(i);
                        String accY = accArrayY.get(i);
                        String accZ = accArrayZ.get(i);
                        String magT = magTime.get(i);
                        String magX = magArrayX.get(i);
                        String magY = magArrayY.get(i);
                        String magZ = magArrayZ.get(i);
                        String gyroT = gyroTime.get(i);
                        String gyroX = gyroArrayX.get(i);
                        String gyroY = gyroArrayY.get(i);
                        String gyroZ = gyroArrayZ.get(i);

                        if(accT != magT || accT != gyroT){
                            Log.e("FAILED SAME TIME CHECK", "Times were different in arrays");
                        }

                        Log.d("acc", accT+" "+accX+" "+accY+" "+accZ);
                        try {
                            outFD_acc.write(accT.getBytes());
                            outFD_gyro.write(gyroT.getBytes());
                            outFD_mag.write(magT.getBytes());

                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        try {
                            outFD_acc.write(",".getBytes());
                            outFD_gyro.write(",".getBytes());
                            outFD_mag.write(",".getBytes());

                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        try {
                            outFD_acc.write(accX.getBytes());
                            outFD_gyro.write(gyroX.getBytes());
                            outFD_mag.write(magX.getBytes());

                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        try {
                            outFD_acc.write(",".getBytes());
                            outFD_gyro.write(",".getBytes());
                            outFD_mag.write(",".getBytes());
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        try {
                            outFD_acc.write(accY.getBytes());
                            outFD_gyro.write(gyroY.getBytes());
                            outFD_mag.write(magY.getBytes());

                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        try {
                            outFD_acc.write(",".getBytes());
                            outFD_gyro.write(",".getBytes());
                            outFD_mag.write(",".getBytes());
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        try {
                            outFD_acc.write(accZ.getBytes());
                            outFD_gyro.write(gyroZ.getBytes());
                            outFD_mag.write(magZ.getBytes());

                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        try {
                            outFD_acc.write(",".getBytes());
                            outFD_gyro.write(",".getBytes());
                            outFD_mag.write(",".getBytes());
                        } catch (IOException e) {
                            e.printStackTrace();
                        }

                        try {
                            outFD_acc.write("\n".getBytes());
                            outFD_gyro.write("\n".getBytes());
                            outFD_mag.write("\n".getBytes());

                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }

                    try {
                        outFD_acc.close();
                        outFD_gyro.close();
                        outFD_mag.close();

                    } catch (IOException e) {
                        e.printStackTrace();
                    }
//                    try {
//                        outFD_mag.close();
//                    } catch (IOException e) {
//                        e.printStackTrace();
//                    }


                    accArrayX.clear();
                    accArrayY.clear();
                    accArrayZ.clear();

                    magArrayX.clear();
                    magArrayY.clear();
                    magArrayZ.clear();

                    gyroArrayX.clear();
                    gyroArrayY.clear();
                    gyroArrayZ.clear();
                    file_counter += 1;
                    count = 0;
                    text.setText("Not Ready");
                    break;
            }

        }
    }
}
