import matplotlib.pyplot as plt
from fmpy import *
from fmpy.fmi2 import FMU2Slave
from matplotlib.widgets import Slider
from matplotlib.animation import FuncAnimation
import shutil
import math
import random as rd
import numpy as np
import csv

FMU_FILENAME = "ControledHeaterAuto.fmu"

OUPUTS = {'real': ["ThermalMass.T", "prescribedHeatFlow.Q_flow", "thermalConductor.Q_flow", "SetHeater","SetTemp"], 
          'int': []}

INPUTS = {'real': ["SetHeater","SetTemp"], 'int': []}
INITVAL = {'real': [0.0,270.0], 'int': []}
setVal = INITVAL['real']

start_time = 0.0
stop_time = 1e3
STEP_SIZE = 1.0  

N = 5 #number of training data sets
training_data = np.full((N, 7), None, dtype=object)


for i in range(N):
    train_temp = []
    #definition of outside temperature 
    def outside_temperature(t, params=[(0,100,270,280)]): #intial params is defined as (0,100,270,280), function adds new values every 100 s 
        global train_temp
        segment = int(t // 100)
        start_segment = segment * 100
        end_segment = start_segment + 100
        
        # Generate new set of temperature values
        if segment > len(params):
            params.append((start_segment,end_segment,params[-1][3],rd.randint(270, 300))) 
    
        start_segment, end_segment, start_temp_f, end_temp_f = params[-1]

        
        train_temp.append(start_temp_f + (end_temp_f - start_temp_f) * (t - start_segment) / (end_segment - start_segment))
        
        # Linear interpolation, calculate the temperature for the current time
        return start_temp_f + (end_temp_f - start_temp_f) * (t - start_segment) / (end_segment - start_segment)


    # Initialize the outside temperature
    outside_temp = outside_temperature(start_time)

    #set the desired temperature and hysteresis
    desiredTemp = 293.0
    hyst = 3.0



    def get_value_ref(var_name):
        for variable in model_desc.modelVariables:
            #print(variable.name)
            if variable.name == var_name:
                return variable.valueReference
        raise ValueError(f"Variable '{var_name}' not in the FMU.")


    def controlUnitOutput():

        #set the values for the heater
        on_value = 700.
        off_value = 0.

        #control of the heater (hysteresis)       
        if fmu.getReal(outvarref['real'])[0] < desiredTemp - hyst:
            setVal[0] = on_value
        elif fmu.getReal(outvarref['real'])[0] > desiredTemp + hyst:
            setVal[0] = off_value
            
        return setVal[0]

    def setup_simulation():
        """Initialize and start the FMU simulation."""
        # Load and instantiate FMU
        unzip_dir = extract(FMU_FILENAME)
        fmu = FMU2Slave(guid=model_desc.guid,
                        modelIdentifier=model_desc.coSimulation.modelIdentifier,
                        unzipDirectory=unzip_dir)

        # Initialize FMU
        fmu.instantiate()
        fmu.setupExperiment(startTime=start_time, stopTime=stop_time)
        fmu.enterInitializationMode()
        fmu.exitInitializationMode()
        
        # Variables references
        invarref = {'real': [get_value_ref(var) for var in INPUTS['real']],
                    'int': [get_value_ref(var) for var in INPUTS['int']],
                }
        outvarref = {'real': [get_value_ref(var) for var in OUPUTS['real']],
                    'int': [get_value_ref(var) for var in OUPUTS['int']],
                    }
        
        # Set initial input value
        if len(INPUTS['real'])>0:
            fmu.setReal(vr=invarref['real'], value=INITVAL['real'])
        if len(INPUTS['int'])>0:
            fmu.setInteger(vr=invarref['int'], value=INITVAL['int'])
        
        return fmu, unzip_dir, invarref, outvarref

    def step_simulation(fmu, current_time, dt, invarref, outvarref):
        """Advance the simulation by one step and return new output values."""
        # Update input value
        setVal[0] = controlUnitOutput()
        setVal[1] = outside_temperature(current_time)
        if len(INPUTS['real'])>0:
            fmu.setReal(vr=invarref['real'], value=setVal)
        if len(INPUTS['int'])>0:
            fmu.setInteger(vr=invarref['int'], value=int(setVal))

        # Perform simulation step
        
        fmu.doStep(current_time, dt)
        current_time += dt

        # Read output values
        out = [current_time]
        if len(OUPUTS['real'])>0:
            out.append(fmu.getReal(outvarref['real']))
        else:
            out.append([])

        if len(OUPUTS['int'])>0:
            out.append(fmu.getInteger(outvarref['int']))
        else:
            out.append([])

        # print(out)
        return out


    def runSimulation(fmu, current_time):
        # Initialize data storage
        time_data = []
        output_data = {'real': [[] for _ in OUPUTS['real']],
                    'int': [[] for _ in OUPUTS['int']],
                    }
        
        while current_time < stop_time:
            dt = min(STEP_SIZE, stop_time-current_time)
            current_time, routput_values, ioutput_values= step_simulation(fmu, current_time, dt, invarref, outvarref)
            
            # print(current_time)
            # print(routput_values)
            # print(ioutput_values)
            # Update outside temperature
            outside_temperature(current_time)
            
            # Append new values to data arrays
            time_data.append(current_time)
            for i, var in enumerate(OUPUTS['real']):
                output_data['real'][i].append(routput_values[i])
            

    
    
        return time_data, output_data, 

    if __name__ == "__main__":
        flag_fmuready = False
        try:
            model_desc = read_model_description(FMU_FILENAME)
            
            # Start FMU simulation
            fmu, unzip_dir, invarref, outvarref = setup_simulation()
            flag_fmuready = True

        
            
            # run whole simulation
            current_time = 0.0
            time_data, output_data = runSimulation(fmu, current_time)

            """""
            # Set up figure
            fig, axs = plt.subplots(4, 1)

            # Clear and draw the plots
            iax = 0
            axs[iax].clear()
            axs[iax].plot(time_data, output_data['real'][3], label=OUPUTS['real'][3], linewidth=2)
            axs[iax].set_ylabel(f"{OUPUTS['real'][3]} [W]")
            axs[iax].legend()
            axs[iax].grid(True)
            
            iax += 1
            axs[iax].clear()
            axs[iax].plot(time_data, output_data['real'][0], label=OUPUTS['real'][0], linewidth=2)
            axs[iax].set_ylabel(f"{OUPUTS['real'][0]} [K]")
            axs[iax].legend()
            axs[iax].grid(True)
            axyy = axs[iax].twinx()
            axyy.plot(time_data, [Ti - 273.15 for Ti in output_data['real'][0]], label=None, linewidth=2)
            axyy.set_ylabel(f"{OUPUTS['real'][0]} [C]")
            
            iax += 1
            axs[iax].clear()
            axs[iax].plot(time_data, output_data['real'][1], label=OUPUTS['real'][1], linewidth=2)
            axs[iax].plot(time_data, output_data['real'][2], label=OUPUTS['real'][2], linewidth=2)
            axs[iax].set_xlabel("Time [s]")
            axs[iax].set_ylabel("$Q_{flow} [W]$")
            axs[iax].legend()
            axs[iax].grid(True)
            
            
            iax += 1
            axs[iax].clear()
            axs[iax].plot(time_data, output_data['real'][4],label=OUPUTS['real'][4], linewidth=2)
            axs[iax].set_xlabel("Time [s]")
            axs[iax].set_ylabel("Outside Temperature [K]")
            axs[iax].legend()
            axs[iax].grid(True)
            
            plt.show()
            """
            # Save the training data
            # Create a dictionary to store the training data
            # Append the time data and output data to the training data
            "ThermalMass.T", "prescribedHeatFlow.Q_flow", "thermalConductor.Q_flow", "SetHeater","SetTemp"
           
            training_data[i][0] = time_data
            training_data[i][1] = output_data['real'][0]
            training_data[i][2] = output_data['real'][1]
            training_data[i][3] = output_data['real'][2]
            training_data[i][4] = output_data['real'][3]
            training_data[i][5] = output_data['real'][4]
            training_data[i][6] = train_temp
            # Append the training data to the list
  
        finally:
            # clean up
            if flag_fmuready:
                fmu.terminate()
                fmu.freeInstance()
            shutil.rmtree(unzip_dir, ignore_errors=True)


outside_temp_train = training_data[0][5]
inside_temp_train = training_data[0][1]
# Save the training data to a file


# Save the training data to a .npy
np.save("training_data_out.npy", outside_temp_train)
np.save("training_data_in.npy", inside_temp_train)



        