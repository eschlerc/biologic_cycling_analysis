import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import os
import navani.echem as ec


os.makedirs('plots', exist_ok=True)
filename = glob.glob('*.mpr')[0]
print('Filename: ', filename)
df = ec.echem_file_loader(filename)
print(df.columns)

voltage_limits = [1.4, 4.5]
start_indices = []
end_indices = [] 


# Filter initial OCV part
# data.drop(data[data['Ns']== 0].index, inplace=True)


# Plot discharge curve
plt.figure(figsize=(12,4))
plt.plot(-1*data['(Q-Qo)/mA.h'], data['Ewe/V'])
plt.ylim(voltage_limits)
plt.xlabel('Discharge capacity (mAh)')
plt.ylabel('Voltage (V)')
plt.title('Discharge curve')
plt.savefig('plots/discharge.png')
plt.close()

# Plot voltage time series
plt.figure(figsize=(12,4))
plt.plot(data['time/s'], data['Ewe/V'], '-*')
# plt.ylim([2.5, 3.5])
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.title('Voltage time series')
plt.savefig('plots/voltage_time_series.png')
plt.close()


# Split data into sections based on current to figure out when relaxation is happening
i = 0
while i < len(data):

    while np.array(data['control/V/mA'])[i] < 0.0:
        i=i+1
    start_indices.append(i)

    while np.array(data['control/V/mA'])[i] == 0.0 and i <= len(data):
        i=i+1
        if i == len(data):
            break
    end_indices.append(i)


# Plot relaxation curves for each GITT relaxation
os.makedirs('relaxations', exist_ok=True)
for i in range(len(start_indices)):
    x_data = data['time/s'][start_indices[i]:end_indices[i]]/3600 - data['time/s'][start_indices[i]]/3600
    y_data = data['Ewe/V'][start_indices[i]:end_indices[i]]
    np.savetxt('relaxations/'+str(i)+'.csv', np.transpose([x_data, y_data]), delimiter=',')
    # popt, pcov = curve_fit(func, x_data, y_data, p0=popt, method='lm', maxfev=10000)
    # print(popt)
    # plt.plot(x_data, func(x_data, *popt), 'r-', label='fit')
    plt.plot(x_data, y_data, label='Data', color='black')
    plt.ylim(voltage_limits)
    plt.xlabel('Time (hr)')
    plt.ylabel('Voltage (v)')
    plt.title(str(i+1)+'th relaxation')
    plt.savefig('relaxations/'+str(i)+'.png')
    plt.close()

# Plot pseudo-OCV curve
plt.plot(-1*data['(Q-Qo)/mA.h'][np.array(end_indices)], data['Ewe/V'][np.array(end_indices)])
plt.ylim(voltage_limits)
plt.xlabel('Discharge capacity (mAh)')
plt.ylabel('Voltage (V)')
plt.title('Pseudo-OCV')
plt.savefig('plots/pseudo_OCV.png')
plt.close()
