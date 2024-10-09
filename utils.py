import numpy as np
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import math
import random
import pywt

def slicing(df, starttime, stoptime, date_column, date_index=True):
    if date_index:
        outdf = df[(df.index>=starttime)&(df.index<=stoptime)]
    else: 
        outdf = df[(df[date_column]>=starttime)&(df[date_column]<=stoptime)]
    return outdf


def keep_pos(x):
    if x < 0:
        x=0
    return x

def conversion(level):
    if level == 0:
        Qp = 0
    else:
        D = 1.25 # pipe diameter [m]
        ks = 0.5/1000 # roughness coefficient [m]
        phi = 2*np.pi - 2*np.arccos((level/1000 - D/2)/(D/2))
        Ad = ((phi - np.sin(phi))/8)*D**2
        Pd = phi*D/2
        Rd = (1-np.sin(phi))*D/4
        Ap = (phi - np.sin(phi))/(2*np.pi)
        Rp = 1-np.sin(phi)/phi
        Vp = (1+np.log(Rp)/np.log(3.7*D/ks))*Rp**0.5
        Qp = (1+np.log(Rp)/np.log(3.7*D/ks))*Ap*Rp**0.5

    return Qp

def dataresample(df, index, resolution):
    df.set_index(index, inplace=True)
    df_resampled = df.resample(resolution).asfreq()
    df_interpolated = df_resampled.interpolate(method='linear')
    df_interpolated.reset_index(inplace=True)

    return df_interpolated

def rmse(observed, predicted):
    return math.sqrt(sum((obs - pred) ** 2 for obs, pred in zip(observed, predicted)) / len(observed))

def nash_sutcliffe(observed, simulated):
    """
    To calculate Nash-Sutcliffe coefficient

    parameter:
    observed: stormwater value
    simulated: Qout from simulation

    return:
    Nash-Sutcliffe coefficient
    """
    mean_observed = sum(observed) / len(observed)

    # numerator of Nash-Sutcliffe
    numerator = sum((obs - sim) ** 2 for obs, sim in zip(observed, simulated))

    # denominator of Nash-Sutcliffe
    denominator = sum((obs - mean_observed) ** 2 for obs in observed)

    # calculate Nash-Sutcliffe
    if denominator != 0:
        ns_coefficient = 1 - (numerator / denominator)
    else:
        ns_coefficient = float('nan')  # Avoid error when 0 in denominator

    return ns_coefficient

def count_nonzero_segments(data):
    count = 0
    in_segment = False

    for value in data:
        if value != 0:
            if not in_segment:
                in_segment = True
                count += 1
        else:
            in_segment = False

    return count

def rho(s, n, Z_root, E_max, E_w, beta, s_h, s_w, s_s, s_fc, Ks):

    Z_soil = n * Z_root   # grass root depth[cm]
    eta_w = E_w / Z_soil;        # normalized maximum evapotranspiration [-/day]
    eta = E_max / Z_soil
    m = Ks / (Z_soil * (np.exp(beta * (1 - s_fc)) - 1))

    if s <= s_h:
        return 0
    elif s_h < s <= s_w:
        return eta_w * (s - s_h) / (s_w - s_h)
    elif s_w < s <= s_s:
        return eta_w + (eta - eta_w) * (s - s_w) / (s_s - s_w)
    elif s_s < s <= s_fc:
        return eta
    else:
        return eta + m * (np.exp(beta * (s - s_fc)) - 1)
        # return eta + heavy/Z_soil
    
def model_st(tspan=[0], y0=[0,0,0], precp_values=[0], k=1, frac_rt2s=0.8, frac_tk2rf=0, n=0.45, Z_root=60, \
            E_max=0.5, E_w=0.0625, beta=14.8, s_h=0.19, s_w=0.24, s_s=0.57, s_fc=1, Ks=20,\
            omegas=0, omega0=0, omegat=0, Vtmax=0, heavy=55, lag=0, resolution=10):
    '''
    Simulate the discharge, soil moisture and V_tank given inputs like precipitation, hydraulic conductivity, etc.

    Args:
        tspan (ndarray): Data points to input.
        y0 (ndarray): Initialization of [simulated_discharge, soil_moisture, simulated_Vtank].
        precp_values: 
        k (float): Time constant 
        fraction (float): Part of rooftop surface that serves as tank
        n (float): Soil porosity
        Z_root: Active root depth
        E_max: maximum evapotransporation for a certain plant
        E_w: evapotransporation at wilting point
        beta: soil parameter
        s_h (float): hygroscopic point of soil moisture
        s_w (float): wilting point
        s_s (float): stomatal closure point 
        s_fc (float): soil moisture at field capacity
        Ks: soil hydraulic conductivity
        omegas (float): area of soil surface
        omega0 (float): area of impermeable surface
        Vtmax (float): maximum tank volume

    Returns:
        int or float: The sum of the two numbers.

    '''

     # Initialization
    # y0 = [0, s_s, 0]

    # 0: runoff, 1: soil moisture, 2: 
    y_result = np.zeros((len(tspan), 3))
    y_result[0] = y0
    dydt = np.zeros(3)
    V_runoff  = 0
    V_tank = 0
    # soil reservoir overflow
    Qsoil = np.zeros(len(tspan))
    Qout = np.zeros(len(tspan))
    precp_accu = np.zeros(len(tspan))
    discharge_accu = np.zeros(len(tspan))
    Qsoil_accu = np.zeros(len(tspan))

    discharge_total = 0
    precp_total = 0
    Qp = 0
    precp_lag = precp_values if lag == 0 else np.concatenate((np.zeros(lag), precp_values[:-lag]))

    for i in range(len(tspan)-1):
        # Rainfall intensity at time point i （from mm/10min to m/s)
        p = np.interp(i, tspan, precp_lag)/600000
    

        # Flow produced by soil
        Qs = p * omegas
        # Qs activated when the rain is too heavy or when the soil is saturated
        Qs_function = Qs if (p > heavy) or (y_result[i,1] >= s_fc) else 0
        # Qs_function = Qs if y_result[i,1] >= s_fc else 0
        Qsoil[i] = Qs_function

        # Discharge
        V_runoff = y_result[i,0]*omega0 # Water volume of V0 reservoir at time point i
        Qout[i] = k * V_runoff / 86400

        # Flow produced by tank reservoir
        Vt = frac_tk2rf*omegat*y_result[i, 2]
        Qt = p*frac_tk2rf*omegat if (Vt>Vtmax) else 0

        # Equation 1 for V0
        # dydt[0] = (p*(omega0+(1-frac_rt2s)*omegat) + Qs_function + Qt - Qout[i])/((omega0+omegat)*(1-frac_tk2im))
        dydt[0] = (p*(omega0+(1-frac_tk2rf)*(1-frac_rt2s)*omegat) + Qs_function + Qt - Qout[i])/omega0
        #### Boundary conditions to add ####
        Z_soil = n*Z_root
    
        dydt[1] = p*(1+(1-frac_tk2rf)*frac_rt2s*omegat/omegas)/ (Z_soil / 100) \
                 - rho(y_result[i, 1],n, Z_root, E_max, E_w, beta, s_h, s_w, s_s, s_fc, Ks) / 86400
        
        if frac_tk2rf != 0:
            dydt[2] = p - (Qt+Qp)/(frac_tk2rf*omegat)
        else:
            dydt[2] = 0

        # Next Runoff Head
        y_result[i+1, 0] = y_result[i, 0] + dydt[0]*60*resolution # tspan resolution is 10 min = 600s
        if y_result[i+1, 0] < 0:
            y_result[i+1, 0] = 0
        
        
        # Next Soil Moisture
        y_result[i+1, 1] = y_result[i, 1] + dydt[1]*60*resolution
        if y_result[i+1, 1] >= s_fc:
            y_result[i+1, 1] = s_fc

        # Next Tank Head
        y_result[i+1, 2] = y_result[i, 2] + dydt[2]*60*resolution

        # Accumulated water discharges [m3]
        if i == 0:
            precp_accu[i]= 0 + p*(omega0+omegas+omegat)*60*resolution
            discharge_accu[i] = 0 + Qout[i]*60*resolution
            Qsoil_accu[i] = 0 + Qsoil[i]*60*resolution
        else:
            precp_accu[i] = precp_accu[i-1] + p*(omega0+(1-frac_rt2s)*omegat)*60*resolution
            discharge_accu[i] = discharge_accu[i-1] + Qout[i]*60*resolution
            Qsoil_accu[i] = Qsoil_accu[i-1] + Qsoil[i]*60*resolution

        discharge_total += Qout[i]*60*resolution
        precp_total += p*60*resolution*(omega0+(1-frac_rt2s)*omegat)
    
    Q_peak = np.max(Qout)
    # print('Time constant k = '+str(k)+'/day')
    return y_result, Q_peak, discharge_total, precp_total, discharge_accu, precp_accu, Qsoil_accu, Qsoil, Qout
