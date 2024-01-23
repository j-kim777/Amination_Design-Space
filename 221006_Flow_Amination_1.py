'''Import library'''
import numpy as np
import time
import csv
import pandas as pd
import numba
from numba import jit
from scipy import optimize
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
import joblib
import concurrent.futures

t1 = time.time()

# Initial condition
#ID
v_flowrate = 0.5 #cm3 min-1
Temp_input = 313.0 #[313.0 + 20.0 * float(a) for a in range(6)] #K###########################################################################################################################
C_SM_input = [1000.0]#, 500.0, 750.0, 1000.0] # mol m-3
C_MP_input = [1.0 + 0.5 * float(a) for a in range(9)] #eq
C_DBU_input = [0.0 + 0.2 * float(a) for a in range(11)] #eq
#Inner_diameter = [1.0 + 0.2 * float(a) for a in range(8)]
#temp_input = [313.0 + 20.0 * float(a) for a in range(6)] #K
design_total = [0] * 99
number = 0
for a in C_SM_input:
    for b in C_MP_input:
        for c in C_DBU_input:
            design_total[number] = [a, a * b, a * c]
            number += 1

'''Main Calculation'''
#@jit
def multi_processing(n_multiprocessing):
    def main_calculation():
        for nnn in range(99):
            dpipe_ini = 1.0e-3 * (1.0 + 0.1 * float(n_multiprocessing))
                            
            '''reactor'''
            dpipe = dpipe_ini
            rpipe = dpipe / 2.0
            d_outer = dpipe * 1.3
            r_outer = d_outer / 2.0

            #Residence time: 10 min
            #Required volume
            #Vpipe = v_flowrate * 10.0 #[cm3]
            t_residence = 30.0 #[min]
            surface = np.pi * ((dpipe / 2.0) * (dpipe / 2.0))
            Vpipe = v_flowrate * t_residence / 1.0e+6 #[m3]
            lpipe = Vpipe / surface
            
            '''calculation condition'''
            t = 0.0
            Nl = 100
            dl = lpipe / float(Nl)
            dt = 1.0 ###########################################################################################################################
            #t_final = 20.0 * 60.0
            countnumber = int(1.0 / dt)
    
            '''initial condition'''        
            C_SM_ini = design_total[nnn][0] # [mol m-3]
            C_MP_ini = design_total[nnn][1] # [mol m-3]
            C_DBU_ini = design_total[nnn][2] # [mol m-3]
            temp_ini = Temp_input #design_total[nnn][3] # [K]  #temp_sol_ini = 30.0 + 273.15 #[K]
            temp_cooling_ini = temp_ini
            temp_in = temp_ini
            temp_cooling_in = temp_cooling_ini
            C_DBUH_ini = 0.0 # [mol m-3]
            C_DP_ini = 0.0 # [mol m-3]
            
            '''velocity'''      
            u_cc_min = v_flowrate
            u = u_cc_min * 1.0e-6 / 60.0 / surface #[m/s]
            t_final = lpipe / u * 2.0
            
            '''matrix'''
            innerdiameter = np.full(Nl + 2, dpipe_ini * 1.0e+3) #mm
            time_step = np.zeros(Nl + 2) #min
            for i in range(int(Nl + 2)):
                time_step[i] = float(i) * (dl / u / 60.0)
            temp = np.full(Nl + 2, temp_ini)
            C_SM = np.zeros(Nl + 2)
            C_MP = np.zeros(Nl + 2)
            C_DBU = np.zeros(Nl + 2)
            C_DBUH = np.zeros(Nl + 2)
            C_DP = np.zeros(Nl + 2)
            
            C_SM[0] = C_SM_ini
            C_MP[0] = C_MP_ini
            C_DBU[0] = C_DBU_ini
            
            '''physical property'''
            #https://www.dow.com/content/dam/dcc/documents/en-us/productdatasheet/327/327-00025-01-isopropyl-acetate-tds.pdf?iframe=true
            kb = 1.38 * 1.0e-23
            η_sol = 0.5 * 1.0e-3 # viscosity [Pa s] (iPrOAc, 20℃)
            ρ_sol = 870.0 # density [kg m-3] (iPrOAc)
            Cp_sol = 1.904 * 1000.0 # [J kg-1 K-1] #酢酸エチル #https://www.matweb.com/search/datasheet_print.aspx?matguid=c634566b56e04467bbfc09ffd3434ebb
            k_thermal = 0.137 # [J s-1 m-1 K-1] #酢酸エチル https://wikitech.info/1230
            k_thermal_PTFE = 0.25 #https://www.packing.co.jp/PTFE/ptfe_ippantokusei1.htm#an3
            Nu = 3.66
            h_thermal = Nu * k_thermal / dpipe
            #print(h_thermal)
            R_constant = 8.31
            rb = 800.0 * 1.0e-12
            α = k_thermal / ρ_sol / Cp_sol
            dH = 129.0 * 1000.0 #[J mol-1]
            Re = ρ_sol * u * dpipe / η_sol
            
            '''kinetic parameter'''
            k_coef1 = 0.001308358471035654
            k_coef2 = 0.001280735406004463
            k_coef3 = 8.69785313725282e-09
            k_coef4 = -1.4910056433447654e-08
            k_coef5 = 0.0007642694352623905
            k_coef6 = 0.0007372303821965298
            A_1 = Re * k_coef1 + k_coef2
            A_2 = Re * k_coef3 + k_coef4
            A_3 = Re * k_coef5 + k_coef6
            E_1 = 4.79 * 1.0e+4
            E_2 = 1.98 * 1.0e+4
            E_3 = 4.83 * 1.0e+4
            
            '''time count'''
            t_con = 0
            t_num = 0.0
            t_mod = int(1.0 / dt)
    
            t = 0.0
            t_print = 0
            j = 0
    
            while t <= t_final:        
                Amat_temp = np.zeros((Nl + 2, Nl + 2))
                Amat_SM = np.zeros((Nl + 2, Nl + 2))
                Amat_MP = np.zeros((Nl + 2, Nl + 2))
                Amat_DBU = np.zeros((Nl + 2, Nl + 2))
                Amat_DBUH = np.zeros((Nl + 2, Nl + 2))
                Amat_DP = np.zeros((Nl + 2, Nl + 2))
                
                bvec_temp = np.zeros((Nl + 2))
                bvec_SM = np.zeros((Nl + 2))
                bvec_MP = np.zeros((Nl + 2))
                bvec_DBU = np.zeros((Nl + 2))
                bvec_DBUH = np.zeros((Nl + 2))
                bvec_DP = np.zeros((Nl + 2))
        
                Amat_temp[0, 0] = 1.0
                Amat_temp[Nl + 1, Nl] = - 1.0
                Amat_temp[Nl + 1, Nl + 1] = 1.0
        
                Amat_SM[0, 0] = 1.0
                Amat_SM[Nl + 1, Nl] = - 1.0
                Amat_SM[Nl + 1, Nl + 1] = 1.0
        
                Amat_MP[0, 0] = 1.0
                Amat_MP[Nl + 1, Nl] = - 1.0
                Amat_MP[Nl + 1, Nl + 1] = 1.0
                
                Amat_DBU[0, 0] = 1.0
                Amat_DBU[Nl + 1, Nl] = - 1.0
                Amat_DBU[Nl + 1, Nl + 1] = 1.0
                
                Amat_DBUH[0, 0] = 1.0
                Amat_DBUH[Nl + 1, Nl] = - 1.0
                Amat_DBUH[Nl + 1, Nl + 1] = 1.0
                
                Amat_DP[0, 0] = 1.0
                Amat_DP[Nl + 1, Nl] = - 1.0
                Amat_DP[Nl + 1, Nl + 1] = 1.0
        
                bvec_temp[0] = temp_ini
                bvec_temp[Nl + 1] = 0.0
        
                bvec_SM[0] = C_SM_ini
                bvec_SM[Nl + 1] = 0.0
        
                bvec_MP[0] = C_MP_ini
                bvec_MP[Nl + 1] = 0.0
                
                bvec_DBU[0] = C_DBU_ini
                bvec_DBU[Nl + 1] = 0.0
                
                bvec_DBUH[0] = C_DBUH_ini
                bvec_DBUH[Nl + 1] = 0.0
        
                bvec_DP[0] = C_DP_ini
                bvec_DP[Nl + 1] = 0.0
        
                for i in range(1, Nl + 1):
                    k_1 = A_1 * np.exp(- E_1 / R_constant / temp[i])
                    k_2 = A_2 * np.exp(- E_2 / R_constant / temp[i])
                    k_3 = A_3 * np.exp(- E_3 / R_constant / temp[i])
                    diffusivity = kb * temp[i] / 6.0 / np.pi / η_sol / rb
                    sss = diffusivity * dt / dl / dl
                    xxx = u * dt / dl
                    ttt = α * dt / dl / dl
                    yyy = u * dt / dl
                    # Temp
                    Amat_temp[i, i - 1] = - ttt - yyy
                    Amat_temp[i, i] = 2.0 * ttt + 2.0 + yyy
                    Amat_temp[i, i + 1] = - ttt
                    
                    # SM
                    Amat_SM[i, i - 1] = - sss - xxx
                    Amat_SM[i, i] = 2.0 * sss + 2.0 + xxx
                    Amat_SM[i, i + 1] = - sss
        
                    # MP
                    Amat_MP[i, i - 1] = - sss - xxx
                    Amat_MP[i, i] = 2.0 * sss + 2.0 + xxx
                    Amat_MP[i, i + 1] = - sss
                    
                    # DBU
                    Amat_DBU[i, i - 1] = - sss - xxx
                    Amat_DBU[i, i] = 2.0 * sss + 2.0 + xxx
                    Amat_DBU[i, i + 1] = - sss
                    
                    # DBUH
                    Amat_DBUH[i, i - 1] = - sss - xxx
                    Amat_DBUH[i, i] = 2.0 * sss + 2.0 + xxx
                    Amat_DBUH[i, i + 1] = - sss
        
                    # Product
                    Amat_DP[i, i - 1] = - sss - xxx
                    Amat_DP[i, i] = 2.0 * sss + 2.0 + xxx
                    Amat_DP[i, i + 1] = - sss
                    
                    if C_DBU[i] >= 0.0:
                        bvec_temp[i] = temp[i - 1] * (ttt + yyy) + temp[i] * (2.0 - 2.0 * ttt - yyy) + temp[i + 1] * ttt + 2.0 * (k_1 * C_SM[i] * C_MP[i] + k_2 * C_SM[i] * C_MP[i] * C_DBU[i] + k_3 * C_SM[i] * C_MP[i] * C_DBUH[i]) * dH * dt / ρ_sol / Cp_sol - 2.0 * dt / ρ_sol / Cp_sol * 2.0 * (temp[i] - temp_in) / rpipe / rpipe / (1.0 / h_thermal / rpipe + 1.0 / k_thermal_PTFE * np.log(r_outer / rpipe))
                        bvec_SM[i] = C_SM[i - 1] * (sss + xxx) + C_SM[i] * (2.0 - 2.0 * sss - xxx) + C_SM[i + 1] * sss - 2.0 * k_1 * C_SM[i] * C_MP[i] * dt - 2.0 * k_2 * C_SM[i] * C_MP[i] * C_DBU[i] * dt - 2.0 * k_3 * C_SM[i] * C_MP[i] * C_DBUH[i] * dt
                        bvec_MP[i] = C_MP[i - 1] * (sss + xxx) + C_MP[i] * (2.0 - 2.0 * sss - xxx) + C_MP[i + 1] * sss - 2.0 * k_1 * C_SM[i] * C_MP[i] * dt - 2.0 * k_2 * C_SM[i] * C_MP[i] * C_DBU[i] * dt - 2.0 * k_3 * C_SM[i] * C_MP[i] * C_DBUH[i] * dt
                        bvec_DBU[i] = C_DBU[i - 1] * (sss + xxx) + C_DBU[i] * (2.0 - 2.0 * sss - xxx) + C_DBU[i + 1] * sss - 2.0 * k_1 * C_SM[i] * C_MP[i] * dt - 2.0 * k_2 * C_SM[i] * C_MP[i] * C_DBU[i] * dt - 2.0 * k_3 * C_SM[i] * C_MP[i] * C_DBUH[i] * dt
                        bvec_DBUH[i] = C_DBUH[i - 1] * (sss + xxx) + C_DBUH[i] * (2.0 - 2.0 * sss - xxx) + C_DBUH[i + 1] * sss + 2.0 * k_1 * C_SM[i] * C_MP[i] * dt + 2.0 * k_2 * C_SM[i] * C_MP[i] * C_DBU[i] * dt + 2.0 * k_3 * C_SM[i] * C_MP[i] * C_DBUH[i] * dt
                        bvec_DP[i] = C_DP[i - 1] * (sss + xxx) + C_DP[i] * (2.0 - 2.0 * sss - xxx) + C_DP[i + 1] * sss + 2.0 * k_1 * C_SM[i] * C_MP[i] * dt + 2.0 * k_2 * C_SM[i] * C_MP[i] * C_DBU[i] * dt + 2.0 * k_3 * C_SM[i] * C_MP[i] * C_DBUH[i] * dt
                    
                    else:
                        bvec_temp[i] = temp[i - 1] * (ttt + yyy) + temp[i] * (2.0 - 2.0 * ttt - yyy) + temp[i + 1] * ttt + 2.0 * (k_1 * C_SM[i] * C_MP[i] + k_2 * C_SM[i] * C_MP[i] * C_DBU[i] + k_3 * C_SM[i] * C_MP[i] * C_DBUH[i]) * dH * dt / ρ_sol / Cp_sol - 2.0 * dt / ρ_sol / Cp_sol * 2.0 * (temp[i] - temp_in) / rpipe / rpipe / (1.0 / h_thermal / rpipe + 1.0 / k_thermal_PTFE * np.log(r_outer / rpipe))
                        bvec_SM[i] = C_SM[i - 1] * (sss + xxx) + C_SM[i] * (2.0 - 2.0 * sss - xxx) + C_SM[i + 1] * sss - 2.0 * k_1 * C_SM[i] * C_MP[i] * dt - 2.0 * k_3 * C_SM[i] * C_MP[i] * C_DBUH[i] * dt
                        bvec_MP[i] = C_MP[i - 1] * (sss + xxx) + C_MP[i] * (2.0 - 2.0 * sss - xxx) + C_MP[i + 1] * sss - 2.0 * 2.0 * k_1 * C_SM[i] * C_MP[i] * dt - 2.0 * 2.0 * k_3 * C_SM[i] * C_MP[i] * C_DBUH[i] * dt
                        bvec_DBU[i] = C_DBU[i - 1] * (sss + xxx) + C_DBU[i] * (2.0 - 2.0 * sss - xxx) + C_DBU[i + 1] * sss
                        bvec_DBUH[i] = C_DBUH[i - 1] * (sss + xxx) + C_DBUH[i] * (2.0 - 2.0 * sss - xxx) + C_DBUH[i + 1] * sss
                        bvec_DP[i] = C_DP[i - 1] * (sss + xxx) + C_DP[i] * (2.0 - 2.0 * sss - xxx) + C_DP[i + 1] * sss + 2.0 * k_1 * C_SM[i] * C_MP[i] * dt + 2.0 * k_3 * C_SM[i] * C_MP[i] * C_DBUH[i] * dt
        
                temp = np.linalg.solve(Amat_temp, bvec_temp)
                C_SM = np.linalg.solve(Amat_SM, bvec_SM)
                C_MP = np.linalg.solve(Amat_MP, bvec_MP)
                C_DBU = np.linalg.solve(Amat_DBU, bvec_DBU)
                C_DBUH = np.linalg.solve(Amat_DBUH, bvec_DBUH)
                C_DP = np.linalg.solve(Amat_DP, bvec_DP)
        
                t_old = t; t = t_old + dt
                t_con += 1
                t_print += 1
            
            innerdiameter_df = pd.DataFrame(innerdiameter, columns = ['ID_mm'])
            temp_df = pd.DataFrame(temp, columns = ['temp_K'])
            C_SM_df = pd.DataFrame(C_SM, columns = ['C_SM'])
            C_MP_df = pd.DataFrame(C_MP, columns = ['C_MP'])
            C_DBU_df = pd.DataFrame(C_DBU, columns = ['C_DBU'])
            C_DBUH_df = pd.DataFrame(C_DBUH, columns = ['C_DBUH'])
            C_DP_df = pd.DataFrame(C_DP, columns = ['C_DP'])
            time_step_df = pd.DataFrame(time_step, columns = ['Time_min'])
            
            data_all = pd.concat([innerdiameter_df, time_step_df, temp_df, C_SM_df, C_MP_df, C_DBU_df, C_DBUH_df, C_DP_df], axis = 1)
            
            #print(data_all)
            data_all.to_csv('f_' + str(int(u_cc_min)) + '_ID_' + str(int(n_multiprocessing)) + '_temp_' + str(int(temp_ini)) + '_C_MP_' + str(int(C_MP_ini)) + '_C_DBU_' + str(int(C_DBU_ini)) + '.csv')

    '''Main'''
    main_calculation()

'''Multiprocessing'''
if __name__ == '__main__':
    n_multiprocessing = [i for i in range(15)]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(multi_processing, n_multiprocessing)

    elapsed_time = time.time() - t1
    print("elapsed_time:{0}".format(elapsed_time)+"[sec]")