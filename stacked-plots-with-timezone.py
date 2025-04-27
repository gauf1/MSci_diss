import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import os
import glob
import seaborn as sns
from scipy.optimize import curve_fit
import pytz
utc=pytz.timezone('UTC')
local_tz=pytz.timezone('Canada/Yukon')

# Set figure style
plt.style.use('seaborn-v0_8-whitegrid')

# Ffunction to add sunrise/sunset shading to a plot
def add_sunrise_sunset(axis, sunrises, sunsets, min_time, max_time):
    
  
    sunrises = sorted(sunrises)
    sunsets = sorted(sunsets)
    

    if not sunrises or sunrises[0] > min_time:
        # finds the last sunset before min_time
        last_sunset_before = None
        for ss in reversed(sunsets):
            if ss < min_time:
                last_sunset_before = ss
                break
        
      
        if last_sunset_before is not None:
            # finds the first sunrise after this sunset
            for sr in sunrises:
                if sr > last_sunset_before:
                    # add night shading from plot start to this sunrise
                    axis.axvspan(min_time, sr, alpha=0.1, color='#808080')
                    break
        else:
          
            if sunrises and sunrises[0] > min_time:
                # shade as night from min_time to first sunrise
                axis.axvspan(min_time, sunrises[0], alpha=0.1, color='#808080')
    
    #  process all the complete day/night cycles
    for i in range(len(sunrises)):
        sr = sunrises[i]
        if sr > max_time:
            break
            
        # finds the next sunset
        next_sunset = None
        for ss in sunsets:
            if ss > sr:
                next_sunset = ss
                break
                
        if next_sunset is None:
            axis.axvspan(sr, max_time, alpha=0.1, color='#e2e2e2')
            break
            
        # adds day shading
        end_point = min(next_sunset, max_time)
        axis.axvspan(sr, end_point, alpha=0.1, color='#e2e2e2')
        
        if next_sunset >= max_time:
            break
            
     
        next_sr = None
        for sunrise in sunrises:
            if sunrise > next_sunset:
                next_sr = sunrise
                break
                
        if next_sr is None:
            axis.axvspan(next_sunset, max_time, alpha=0.1, color='#808080')
            break
            
        end_point = min(next_sr, max_time)
        axis.axvspan(next_sunset, end_point, alpha=0.1, color='#808080')
    
    if sunsets and sunsets[-1] < max_time:

        last_sr_after = None
        for sr in sunrises:
            if sr > sunsets[-1] and sr < max_time:
                last_sr_after = sr
                break
                
        if last_sr_after is not None:
            axis.axvspan(sunsets[-1], last_sr_after, alpha=0.1, color='#808080')
            axis.axvspan(last_sr_after, max_time, alpha=0.1, color='#e2e2e2')
        else:
            axis.axvspan(sunsets[-1], max_time, alpha=0.1, color='#808080')

# function for EC calibration fitting
def fit_func(x, a, b, c):
    return a / (x + b) + c
    
    

def get_fixed_y_limits():
    

    limits = {
        # air temperature
        'air_temp': [-18, 25],
        
        # air pressure
        'air_pressure': [785, 825],
        
        # datalogger voltage
        'voltage': [2650, 4000],
        
        # water pressure
        'water_pressure': [30, 160],
        
        # electrical conductivity
        'conductivity': [0, 230],
        
        # wurst temperature
        'cryo_temp': [-0.15, 0.05]
    }
    
    return limits


# function to get date range for a specific week
def get_week_date_range(start_date, week_number):

 
    week_start = start_date + datetime.timedelta(days=7 * week_number)
    week_start = week_start.replace(hour=0, minute=0, second=0)
    
    week_end = week_start + datetime.timedelta(days=7) - datetime.timedelta(seconds=1)
    week_end = week_end.replace(hour=23, minute=59, second=59)
    
    return week_start, week_end

# function to load and fit EC calibration curves
def load_ec_calibration_curves(working_directory):
    calibration_params = {}
    wursts = [2, 4, 7, 8]
    
    for wurst in wursts:
        calibration_file = f"{working_directory}/data/calibration/wurst{wurst}_ec_calibration.csv"
        
        if not os.path.exists(calibration_file):
            print(f"Warning: Calibration file for Wurst {wurst} not found: {calibration_file}")
            calibration_params[wurst] = None
            continue
            
        try:
            # load calibration data
            calibration_data = pd.read_csv(calibration_file)
            

            x = calibration_data.iloc[:, 0].values  # first column (mA)
            y = calibration_data.iloc[:, 1].values  # second column (µS)
            
            # no negative or zero values in x
            if np.any(x <= 0):
                x += np.abs(np.min(x)) + 1.01  
            

            a_guess = (max(y) - min(y)) * np.mean(x)
            b_guess = np.mean(x) * 0.1
            c_guess = np.mean(y)
            initial_guess = [a_guess, b_guess, c_guess]
            
            #bounds
            lower_bounds = [0.01, -max(x), min(y) - 10]
            upper_bounds = [np.inf, max(x), max(y) + 10]
            
            #fit model
            popt, _ = curve_fit(fit_func, x, y, p0=initial_guess, 
                               bounds=(lower_bounds, upper_bounds), maxfev=50000)
            

            calibration_params[wurst] = popt
            print(f"Fitted calibration for Wurst {wurst}: y={popt[0]:.2f}/(x+{popt[1]:.2f})+{popt[2]:.2f}")
            
        except Exception as e:
            print(f"Error fitting calibration curve for Wurst {wurst}: {e}")
            calibration_params[wurst] = None
    
    return calibration_params

#apply EC calibration to a data frame
def apply_ec_calibration(df, uid_to_wurst, calibration_params):
    if 'ec' not in df.columns:
        return df
    

    df_cal = df.copy()
    
    # add a new column for calibrated EC
    df_cal['ec_calibrated'] = df_cal['ec']
    
    # apply calibration for each wurst
    for uid, wurst_num in uid_to_wurst.items():
        if wurst_num in calibration_params and calibration_params[wurst_num] is not None:
            # get the calibration parameters
            a, b, c = calibration_params[wurst_num]
            

            mask = df_cal['UID'] == uid
            

            df_cal.loc[mask, 'ec_calibrated'] = fit_func(df_cal.loc[mask, 'ec'], a, b, c)
    
    return df_cal

# function to create the stacked plot 
def create_stacked_plot(weather_data, all_data, wurst2, wurst4, wurst7, wurst8, min_time, max_time, y_limits, sunrises=None, sunsets=None):

    fig, axs = plt.subplots(6, 1, figsize=(14, 18), sharex=True)
    

    for ax in axs:
        ax.xaxis_date(local_tz)
    

    min_time_local = min_time.astimezone(local_tz) if min_time.tzinfo else utc.localize(min_time).astimezone(local_tz)
    max_time_local = max_time.astimezone(local_tz) if max_time.tzinfo else utc.localize(max_time).astimezone(local_tz)
    
    title = f'Environmental Data Visualization (Local Time - Yukon)\n{min_time_local.strftime("%Y-%m-%d")} to {max_time_local.strftime("%Y-%m-%d")}'
    fig.suptitle(title, fontsize=16, y=0.98)  # Moved title up to avoid overlap
    
    # apply colors
    colours = sns.color_palette("husl", 8)
    wurst2_colour = colours[1]
    wurst4_colour = colours[3]
    wurst7_colour = colours[6]
    wurst8_colour = colours[7]
    temp_colour = colours[4]
    humidity_colour = colours[5]
    
    # hour locator for x-axis
    hour_locator = 24  # 1 day - will show a tick for each day
 
    ms = 3
    

    weather_time_utc = []
    
    if not weather_data.empty:
        for timestamp in weather_data['datetime']:

            if isinstance(timestamp, datetime.datetime):
                dt_utc = utc.localize(timestamp) if timestamp.tzinfo is None else timestamp
            else:
                dt_utc = utc.localize(datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S'))
            weather_time_utc.append(dt_utc)
    
    wurst2_time_utc = []
    wurst4_time_utc = []
    wurst7_time_utc = []
    wurst8_time_utc = []
    
    if not wurst2.empty:
        for timestamp in wurst2['time']:
            if isinstance(timestamp, datetime.datetime):
                dt_utc = utc.localize(timestamp) if timestamp.tzinfo is None else timestamp
            else:
                dt_utc = utc.localize(datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S'))
            wurst2_time_utc.append(dt_utc)
    
    if not wurst4.empty:
        for timestamp in wurst4['time']:
            if isinstance(timestamp, datetime.datetime):
                dt_utc = utc.localize(timestamp) if timestamp.tzinfo is None else timestamp
            else:
                dt_utc = utc.localize(datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S'))
            wurst4_time_utc.append(dt_utc)
    
    if not wurst7.empty:
        for timestamp in wurst7['time']:
            if isinstance(timestamp, datetime.datetime):
                dt_utc = utc.localize(timestamp) if timestamp.tzinfo is None else timestamp
            else:
                dt_utc = utc.localize(datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S'))
            wurst7_time_utc.append(dt_utc)
    
    if not wurst8.empty:
        for timestamp in wurst8['time']:
            if isinstance(timestamp, datetime.datetime):
                dt_utc = utc.localize(timestamp) if timestamp.tzinfo is None else timestamp
            else:
                dt_utc = utc.localize(datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S'))
            wurst8_time_utc.append(dt_utc)
    

    if sunrises is not None and sunsets is not None and len(sunrises) > 0 and len(sunsets) > 0:
        padded_min_time = min_time_local - datetime.timedelta(days=7)  
        padded_max_time = max_time_local + datetime.timedelta(days=7)  
        
        week_sunrises = [sr for sr in sunrises if padded_min_time <= sr <= padded_max_time]
        week_sunsets = [ss for ss in sunsets if padded_min_time <= ss <= padded_max_time]
        
        week_sunrises.sort()
        week_sunsets.sort()
        
        # only proceed if we have data points
        if len(week_sunrises) > 0 or len(week_sunsets) > 0:
            # Create custom patches for the legend
            from matplotlib.patches import Patch
            day_patch = Patch(color='#e2e2e2', alpha=0.1, label='Day')
            night_patch = Patch(color='#808080', alpha=0.1, label='Night')
            
            # add shading to each subplot
            for ax in axs:
                add_sunrise_sunset(ax, week_sunrises, week_sunsets, min_time_local, max_time_local)
    
    # 1. air temperature plot
    if 'Temperature_20339014_°C' in weather_data.columns and not weather_data.empty:
        axs[0].plot(weather_time_utc, weather_data['Temperature_20339014_°C'], 
                   color=temp_colour, linewidth=1.5, marker='.', markersize=ms, label='Air Temperature')
        axs[0].set_ylabel('Air Temp (°C)')
        axs[0].axhline(y=0, color='black', linestyle='--', alpha=0.5)  
        # use fixed y-axis limits for air temperature
        axs[0].set_ylim(y_limits['air_temp'])
        
        # display original min/max if available in the current time window
        if not weather_data.empty:
            min_temp = weather_data['Temperature_20339014_°C'].min()
            max_temp = weather_data['Temperature_20339014_°C'].max()
            axs[0].set_title(f'Air Temperature (min: {min_temp:.2f}°C, max: {max_temp:.2f}°C)', pad=3)  # Reduced padding
        else:
            axs[0].set_title('Air Temperature', pad=3)
            
        # add day/night to legend
        if sunrises is not None and sunsets is not None and len(sunrises) > 0 and len(sunsets) > 0:
            handles, labels = axs[0].get_legend_handles_labels()
            day_patch = Patch(color='#e2e2e2', alpha=0.1, label='Day')
            night_patch = Patch(color='#808080', alpha=0.1, label='Night')
            handles.extend([day_patch, night_patch])
            axs[0].legend(handles=handles, loc='best', fontsize=8)
        else:
            axs[0].legend(loc='best', fontsize=8)
    
    # 2. air Pressure plot
    if 'Pressure_20290338_mbar' in weather_data.columns and not weather_data.empty:
        axs[1].plot(weather_time_utc, weather_data['Pressure_20290338_mbar'], 
                   color='darkgreen', linewidth=1.5, marker='.', markersize=ms, label='Air Pressure')
        axs[1].set_ylabel('Air Pressure (mBar)')
        
        # use fixed y-axis limits for air pressure
        axs[1].set_ylim(y_limits['air_pressure'])
        
        # display original min/max if available 
        if not weather_data.empty:
            min_press = weather_data['Pressure_20290338_mbar'].min()
            max_press = weather_data['Pressure_20290338_mbar'].max()
            axs[1].set_title(f'Air Pressure (min: {min_press:.2f} mBar, max: {max_press:.2f} mBar)', pad=3)
        else:
            axs[1].set_title('Air Pressure', pad=3)
            
        # add legend for air pressure
        axs[1].legend(loc='best', fontsize=8)
    
    # 3. datalogger voltage plot
    if 'logger_voltage' in all_data.columns:
        # plot different wursts with their colours
        if not wurst2.empty:
            axs[2].plot(wurst2_time_utc, wurst2['logger_voltage'], 
                       color=wurst2_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 2 - 138m')
        
        if not wurst4.empty:
            axs[2].plot(wurst4_time_utc, wurst4['logger_voltage'], 
                       color=wurst4_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 4 - 157.8m')
        
        if not wurst7.empty:
            axs[2].plot(wurst7_time_utc, wurst7['logger_voltage'], 
                       color=wurst7_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 7 - 87.19m')
        
        if not wurst8.empty:
            axs[2].plot(wurst8_time_utc, wurst8['logger_voltage'], 
                       color=wurst8_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 8 - 36.89m')
        
        axs[2].set_ylabel('Voltage')
        
        # use fixed y-axis limits for voltage
        axs[2].set_ylim(y_limits['voltage'])
        
        # display original min/max
        if not all_data.empty:
            min_v = all_data['logger_voltage'].min()
            max_v = all_data['logger_voltage'].max()
            axs[2].set_title(f'Datalogger Voltage (min: {min_v:.0f}, max: {max_v:.0f})', pad=3)
        else:
            axs[2].set_title('Datalogger Voltage', pad=3)
            
        # add day/night to legend 
        if sunrises is not None and sunsets is not None and len(sunrises) > 0 and len(sunsets) > 0:
            handles, labels = axs[2].get_legend_handles_labels()
            day_patch = Patch(color='#e2e2e2', alpha=0.1, label='Day')
            night_patch = Patch(color='#808080', alpha=0.1, label='Night')
            handles.extend([day_patch, night_patch])
            axs[2].legend(handles=handles, loc='best', fontsize=8)
        else:
            axs[2].legend(loc='best', fontsize=8)
    
    # 4. mH2O (derived from pressure) plot
    if 'pressure' in all_data.columns:
   
        
        # Plot separate wursts with different colours
        if not wurst2.empty:
            wurst2_mh2o = wurst2['pressure'] * 10.2
            wurst2_filtered = wurst2[wurst2_mh2o > -100]
            if not wurst2_filtered.empty:
                wurst2_filtered_times = []
                for timestamp in wurst2_filtered['time']:
                    if isinstance(timestamp, datetime.datetime):
                        dt_utc = utc.localize(timestamp) if timestamp.tzinfo is None else timestamp
                    else:
                        dt_utc = utc.localize(datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S'))
                    wurst2_filtered_times.append(dt_utc)
                
                axs[3].plot(wurst2_filtered_times, wurst2_filtered['pressure'] * 10.2, 
                           color=wurst2_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 2 - 138m')
        
        if not wurst4.empty:
            wurst4_mh2o = wurst4['pressure'] * 10.2
            wurst4_filtered = wurst4[wurst4_mh2o > -100]
            if not wurst4_filtered.empty:
                wurst4_filtered_times = []
                for timestamp in wurst4_filtered['time']:
                    if isinstance(timestamp, datetime.datetime):
                        dt_utc = utc.localize(timestamp) if timestamp.tzinfo is None else timestamp
                    else:
                        dt_utc = utc.localize(datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S'))
                    wurst4_filtered_times.append(dt_utc)
                
                axs[3].plot(wurst4_filtered_times, wurst4_filtered['pressure'] * 10.2, 
                           color=wurst4_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 4 - 157.8m')
        
        if not wurst7.empty:
            wurst7_mh2o = wurst7['pressure'] * 10.2
            wurst7_filtered = wurst7[wurst7_mh2o > -100]
            if not wurst7_filtered.empty:
                wurst7_filtered_times = []
                for timestamp in wurst7_filtered['time']:
                    if isinstance(timestamp, datetime.datetime):
                        dt_utc = utc.localize(timestamp) if timestamp.tzinfo is None else timestamp
                    else:
                        dt_utc = utc.localize(datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S'))
                    wurst7_filtered_times.append(dt_utc)
                
                axs[3].plot(wurst7_filtered_times, wurst7_filtered['pressure'] * 10.2, 
                           color=wurst7_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 7 - 87.19m')
        
        if not wurst8.empty:
            wurst8_mh2o = wurst8['pressure'] * 10.2
            wurst8_filtered = wurst8[wurst8_mh2o > -100]
            if not wurst8_filtered.empty:
                wurst8_filtered_times = []
                for timestamp in wurst8_filtered['time']:
                    if isinstance(timestamp, datetime.datetime):
                        dt_utc = utc.localize(timestamp) if timestamp.tzinfo is None else timestamp
                    else:
                        dt_utc = utc.localize(datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S'))
                    wurst8_filtered_times.append(dt_utc)
                
                axs[3].plot(wurst8_filtered_times, wurst8_filtered['pressure'] * 10.2, 
                           color=wurst8_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 8 - 36.89m')
                
        axs[3].set_ylabel('Water Pressure (mH2O)')
        
        # use fixed y-axis limits for water pressure
        axs[3].set_ylim(y_limits['water_pressure'])
        
    
        mh2o = all_data['pressure'] * 10.2
        mh2o_filtered = all_data.copy()
        mh2o_filtered['mH2O'] = mh2o
        mh2o_filtered = mh2o_filtered[mh2o_filtered['mH2O'] > -100]  
        
        if not mh2o_filtered.empty:
            min_mh2o = mh2o_filtered['mH2O'].min()
            max_mh2o = mh2o_filtered['mH2O'].max()
            axs[3].set_title(f'Water Pressure (min: {min_mh2o:.2f} mH2O, max: {max_mh2o:.2f} mH2O)', pad=3)
        else:
            axs[3].set_title('Water Pressure', pad=3)
            
        axs[3].legend(loc='best', fontsize=8)
    
    # 5. electrical conductivity plot - using calibrated values if available
    ec_column = 'ec_calibrated' if 'ec_calibrated' in all_data.columns else 'ec'
    
    if ec_column in all_data.columns:
        # Plot separate wursts with different colors
        if not wurst2.empty and ec_column in wurst2.columns:
            axs[4].plot(wurst2_time_utc, wurst2[ec_column], 
                       color=wurst2_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 2 - 138m')
                       
        if not wurst4.empty and ec_column in wurst4.columns:
            axs[4].plot(wurst4_time_utc, wurst4[ec_column], 
                       color=wurst4_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 4 - 157.8m')
                       
        if not wurst7.empty and ec_column in wurst7.columns:
            axs[4].plot(wurst7_time_utc, wurst7[ec_column], 
                       color=wurst7_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 7 - 87.19m')
                       
        if not wurst8.empty and ec_column in wurst8.columns:
            axs[4].plot(wurst8_time_utc, wurst8[ec_column], 
                       color=wurst8_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 8 - 36.89m')
        
        axs[4].set_ylabel('Conductivity (µS)')
        
        # fixed y-axis limits for conductivity
        axs[4].set_ylim(y_limits['conductivity'])
        

        if not all_data.empty and ec_column in all_data.columns:
            min_ec = all_data[ec_column].min()
            max_ec = all_data[ec_column].max()
            
            title_prefix = "Calibrated " if ec_column == 'ec_calibrated' else ""
            axs[4].set_title(f'{title_prefix}Electrical Conductivity (min: {min_ec:.0f} µS, max: {max_ec:.0f} µS)', pad=3)
        else:
            title_prefix = "Calibrated " if ec_column == 'ec_calibrated' else ""
            axs[4].set_title(f'{title_prefix}Electrical Conductivity', pad=3)
            
        axs[4].legend(loc='best', fontsize=8)
    
    # 6. cryowurst temperature plot
    if 'tmp_temp' in all_data.columns:
        # plot separate wursts with different colors
        if not wurst2.empty:
            axs[5].plot(wurst2_time_utc, wurst2['tmp_temp'], 
                       color=wurst2_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 2 - 138m')
                       
        if not wurst4.empty:
            axs[5].plot(wurst4_time_utc, wurst4['tmp_temp'], 
                       color=wurst4_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 4 - 157.8m')
                       
        if not wurst7.empty:
            axs[5].plot(wurst7_time_utc, wurst7['tmp_temp'], 
                       color=wurst7_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 7 - 87.19m')
                       
        if not wurst8.empty:
            axs[5].plot(wurst8_time_utc, wurst8['tmp_temp'], 
                       color=wurst8_colour, linewidth=1.5, marker='.', markersize=ms, label='Wurst 8 - 36.89m')
        
        axs[5].set_ylabel('Cryowurst Temp (°C)')
        axs[5].axhline(y=0, color='black', linestyle='--', alpha=0.5)  # Freezing point reference
        
        # use fixed y-axis limits for cryowurst temperature
        axs[5].set_ylim(y_limits['cryo_temp'])
        
        # display original min/max 
        if not all_data.empty:
            min_cryo = all_data['tmp_temp'].min()
            max_cryo = all_data['tmp_temp'].max()
            axs[5].set_title(f'Cryowurst Temperature (min: {min_cryo:.2f}°C, max: {max_cryo:.2f}°C)', pad=3)
        else:
            axs[5].set_title('Cryowurst Temperature', pad=3)
            
        axs[5].legend(loc='best', fontsize=8)
    
    # set the x-axis limits with offset
    axis_offset = datetime.timedelta(hours=24)
    min_time_local = min_time.astimezone(local_tz) if min_time.tzinfo else utc.localize(min_time).astimezone(local_tz)
    max_time_local = max_time.astimezone(local_tz) if max_time.tzinfo else utc.localize(max_time).astimezone(local_tz)
    
    for ax in axs:
        ax.set_xlim(min_time_local - axis_offset, max_time_local + axis_offset)
    
    # format x-axis with local timezone - include hours when appropriate
    for ax in axs[:-1]:  
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b', tz=local_tz))
        ax.xaxis.set_major_locator(mdates.DayLocator())  # Show daily ticks
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 6, 12, 18]))  # Show 6-hour minor ticks
        ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H', tz=local_tz))  # Hour format for minor ticks
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)  # Smaller font size
        ax.grid(True, linestyle='--', alpha=0.7)
    
    # for bottom subplot, show more details including date and hour
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%d %b', tz=local_tz))
    axs[-1].xaxis.set_major_locator(mdates.DayLocator())
    axs[-1].xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 6, 12, 18]))
    axs[-1].xaxis.set_minor_formatter(mdates.DateFormatter('%H', tz=local_tz))
    plt.setp(axs[-1].get_xticklabels(), rotation=45, ha='right', fontsize=8)  # Smaller font size
    axs[-1].set_xlabel("Date and Time (Local Time - Yukon, UTC-7)", fontsize=10)
    
    # note about data sources and timezone - moved to top right
    fig.text(0.75, 0.01, 
             'Data sources:\n- Air Temperature and Pressure: DataGarrison Weather Station\n' +
             '- Voltage, Water Pressure, Conductivity, Cryowurst Temperature: Satellite & SD data\n' +
             '- Sunrise/Sunset data: NOAA Solar Calculator\n' +
             'Note: All timestamps in local time (Yukon/Mountain Time, UTC-7)',
             fontsize=8, ha='right', va='bottom')
    
    # adjust layout with more bottom space
    plt.tight_layout()
    plt.subplots_adjust(top=0.94, bottom=0.12, hspace=0.3)  # Increased bottom margin
    
    return fig

# function which analyses the entire data set to find its min and max values
def find_global_limits(all_data, weather_data):
    """
    global min and max values for the variables in dataset
    Parameters:
    all_data 
    weather_data 
    
    """
    limits = {}
    
    # air temperature
    if 'Temperature_20339014_°C' in weather_data.columns:
        min_temp = weather_data['Temperature_20339014_°C'].min()
        max_temp = weather_data['Temperature_20339014_°C'].max()

        padding = (max_temp - min_temp) * 0.1
        limits['air_temp'] = [min_temp - padding, max_temp + padding]
    else:
        limits['air_temp'] = [-5, 25]  
    
    # air pressure
    if 'Pressure_20290338_mbar' in weather_data.columns:
        min_press = weather_data['Pressure_20290338_mbar'].min()
        max_press = weather_data['Pressure_20290338_mbar'].max()
        padding = (max_press - min_press) * 0.05
        limits['air_pressure'] = [min_press - padding, max_press + padding]
    else:
        limits['air_pressure'] = [800, 820]  
    # datalogger voltage
    if 'logger_voltage' in all_data.columns:
        min_v = all_data['logger_voltage'].min()
        max_v = all_data['logger_voltage'].max()
        padding = (max_v - min_v) * 0.05
        limits['voltage'] = [min_v - padding, max_v + padding]
    else:
        limits['voltage'] = [3000, 3700]  
    # water pressure (mH2O)
    if 'pressure' in all_data.columns:
        # Convert pressure to mH2O and filter out extreme values
        mh2o = all_data['pressure'] * 10.2
        mh2o_filtered = mh2o[mh2o > -100]  
        
        if not mh2o_filtered.empty:
            min_mh2o = mh2o_filtered.min()
            max_mh2o = mh2o_filtered.max()
            padding = (max_mh2o - min_mh2o) * 0.05
            limits['water_pressure'] = [min_mh2o - padding, max_mh2o + padding]
        else:
            limits['water_pressure'] = [0, 150]  
    else:
        limits['water_pressure'] = [0, 150]  
    
    # EC
    ec_column = 'ec_calibrated' if 'ec_calibrated' in all_data.columns else 'ec'
    
    if ec_column in all_data.columns:
        min_ec = all_data[ec_column].min()
        max_ec = all_data[ec_column].max()
        padding = (max_ec - min_ec) * 0.05
        limits['conductivity'] = [min_ec - padding, max_ec + padding]
    else:
        limits['conductivity'] = [0, 250]  
    
    # wurst temperature
    if 'tmp_temp' in all_data.columns:
        min_cryo = all_data['tmp_temp'].min()
        max_cryo = all_data['tmp_temp'].max()
        padding = (max_cryo - min_cryo) * 0.1
        limits['cryo_temp'] = [min_cryo - padding, max_cryo + padding]
    else:
        limits['cryo_temp'] = [-0.15, 0.05]  
    return limits

# main function
def main():
    # working directory for code to run
    working_directory = 'C:/python/donjek_data_process-main'
    output_path = working_directory+'/plots/'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    # define the timezones we're working with
    utc = pytz.timezone('UTC')
    local_tz = pytz.timezone('Canada/Yukon')
    
    # set up colours
    colours = sns.color_palette("husl", 8)
    
    # import and concatinate both data sets
    satellite_file = working_directory+'/data/processed/satellite_data_processed.csv'
    sd_file = working_directory+'/data/processed/sd_data_july24_processed.csv'
    all_data = pd.read_csv(satellite_file)
    sd_data = pd.read_csv(sd_file)
    all_data = pd.concat([all_data, sd_data], ignore_index=True)
    
    # format datetime
    good_time=[]
    for i in range(0, len(all_data)):
        good_time.append(datetime.datetime.strptime(all_data['time'][i], '%Y-%m-%d %H:%M:%S'))  
    all_data['time']=good_time
    
    # import weather station data 
    weather_list = glob.glob(working_directory+'/data/raw/300234068884730*.txt')
    weather_list.sort()
    weather_data = pd.read_csv(weather_list[0], delimiter='\t', header=2, skipfooter=2, engine='python')
    for index in range(1, len(weather_list)):
        filename = weather_list[index]
        data = pd.read_csv(filename, delimiter='\t', header=2, skipfooter=2, engine='python')
        weather_data = pd.concat([weather_data, data], ignore_index=True)
    
    # format weather data datetime
    good_time=[]
    for index in range(len(weather_data['Date_Time'])):
        good_time.append(datetime.datetime.strptime(weather_data['Date_Time'][index], '%m/%d/%y %H:%M:%S'))
    weather_data['datetime']=good_time
    
    # get sunrise and sunset times at deployment location from NOAA
    # https://gml.noaa.gov/grad/solcalc/calcdetails.html
    print("Loading solar data...")
    solar_data = pd.read_csv(working_directory+'/data/danzhur_solar_2024.csv')    
    
    # convert to datetime and localize to Yukon timezone
    sunrises = [local_tz.localize(datetime.datetime.strptime(solar_data['Date'][i]+' '+solar_data['Sunrise Time (LST)'][i],
                                                          "%Y-%m-%d %H:%M:%S")) for i in range(len(solar_data))]
    sunsets = [local_tz.localize(datetime.datetime.strptime(solar_data['Date'][i]+' '+solar_data['Sunset Time (LST)'][i],
                                                         "%Y-%m-%d %H:%M:%S")) for i in range(len(solar_data))]
    
    # data processing
    # sort by time
    all_data = all_data.sort_values(by=['time'])
    weather_data = weather_data.sort_values(by=['datetime'])
    
    # add wurst depths to dataframe
    all_data['depth'] = np.zeros(len(all_data))
    all_data.loc[all_data['UID']=='cf240002', 'depth'] = 138
    all_data.loc[all_data['UID']=='cf240004', 'depth'] = 157.8
    all_data.loc[all_data['UID']=='cf240007', 'depth'] = 87.19 
    all_data.loc[all_data['UID']=='cf240008', 'depth'] = 36.89 
    
    # correct pressure
    all_data['pressure']=all_data['pressure']-(all_data['logger_pressure']/1e9)
    
    # load EC calibration curves
    print("Loading EC calibration curves...")
    calibration_params = load_ec_calibration_curves(working_directory)
    
    # map UID to wurst number
    uid_to_wurst = {
        'cf240002': 2,
        'cf240004': 4,
        'cf240007': 7,
        'cf240008': 8
    }
    
    # apply EC calibration to the full dataset
    print("Applying EC calibration to the data...")
    all_data = apply_ec_calibration(all_data, uid_to_wurst, calibration_params)
    
    # get fixed y-axis limits that are appropriate for the data
    print("Setting appropriate fixed y-axis limits...")
    fixed_y_limits = get_fixed_y_limits()
    
    # define the overall date range
    start_date = datetime.datetime(2024, 7, 22, 12, 0, 0)  # Original start date
    end_date = datetime.datetime(2024, 11, 30, 23, 59, 59)  # End of November
    
    # make dates timezone-aware
    start_date = utc.localize(start_date)
    end_date = utc.localize(end_date)
    
    # calculate how many weeks we have from start_date to end_date
    days_difference = (end_date - start_date).days
    num_weeks = days_difference // 7 + 1  # +1 to include any partial weeks
    
    print(f"Processing {num_weeks} weeks of data from {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}...")
    
    # process each week
    for week in range(num_weeks):
        # Get date range for this week
        week_start, week_end = get_week_date_range(start_date, week)
        
        # Make week dates timezone-aware
        week_start = utc.localize(week_start) if week_start.tzinfo is None else week_start
        week_end = utc.localize(week_end) if week_end.tzinfo is None else week_end
        
        # skip if the week is beyond our end date
        if week_start > end_date:
            break
            
        print(f"Processing Week {week+1}: {week_start.strftime('%Y-%m-%d')} to {week_end.strftime('%Y-%m-%d')}")
        
        # filter data for this week
        week_all_data = all_data[(all_data['time'] >= week_start.replace(tzinfo=None)) & 
                                (all_data['time'] <= week_end.replace(tzinfo=None))]
        week_weather_data = weather_data[(weather_data['datetime'] >= week_start.replace(tzinfo=None)) & 
                                       (weather_data['datetime'] <= week_end.replace(tzinfo=None))]
        
        # separate out different wursts for this week
        week_wurst2 = week_all_data[week_all_data['UID']=='cf240002']
        week_wurst4 = week_all_data[week_all_data['UID']=='cf240004']
        week_wurst7 = week_all_data[week_all_data['UID']=='cf240007']
        week_wurst8 = week_all_data[week_all_data['UID']=='cf240008']
        
        # output file path for this week
        week_str = week_start.strftime('%Y%m%d')
        output_file = f"{output_path}stacked_plots_week{week+1}_{week_str}.png"
        
        # print data counts
        print(f"  Weather station records: {len(week_weather_data)}")
        print(f"  Combined data records: {len(week_all_data)}")
        print(f"  Wurst 2 (138m) records: {len(week_wurst2)}")
        print(f"  Wurst 4 (157.8m) records: {len(week_wurst4)}")
        print(f"  Wurst 7 (87.19m) records: {len(week_wurst7)}")
        print(f"  Wurst 8 (36.89m) records: {len(week_wurst8)}")
        
        # Filter sunrise and sunset data for this week
        print(f"  Filtering sunrise/sunset data for Week {week+1}...")
        week_start_local = week_start.astimezone(local_tz)
        week_end_local = week_end.astimezone(local_tz)
        
        # create the plot for this week using the fixed y-axis limits
        print(f"  Creating stacked plot for Week {week+1}...")
        fig = create_stacked_plot(
            week_weather_data, 
            week_all_data, 
            week_wurst2, 
            week_wurst4, 
            week_wurst7, 
            week_wurst8, 
            week_start, 
            week_end,
            fixed_y_limits,
            sunrises,
            sunsets
        )
        
        
        # save the figure with error handling
        print(f"  Saving plot to {output_file}...")
        try:
            # make sure output directory exists
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)  # Close the figure to free memory
            print(f"  Successfully saved plot to {output_file}")
        except Exception as e:
            print(f"  Error saving plot: {e}")
            # try with an alternative filename
            alt_output_file = f"{output_path}week{week+1}.png"
            print(f"  Trying alternative filename: {alt_output_file}")
            try:
                fig.savefig(alt_output_file, dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"  Successfully saved plot to {alt_output_file}")
            except Exception as e2:
                print(f"  Error saving with alternative filename: {e2}")
                plt.close(fig)
    
    print("All weeks processed successfully!")

if __name__ == "__main__":
    main()