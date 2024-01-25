'''
Created by Junu Kim (j-kim@pse.t.u-tokyo.ac.jp) at the University of Tokyo
Modified on Jan. 25th, 2024

Publication title: Kinetic study and model-based design space determination for a drug substance flow synthesis using amination reaction via nucleophilic aromatic substitution
'''


'''Import library'''
import numpy as np
import time
import csv
import pandas as pd
import numba
from scipy import optimize
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import joblib
import concurrent.futures

t1 = time.time()

# Initial condition
v_flowrate = 0.5
Temp_input = 313.0
C_SM_input = [1000.0]
C_MP_input = [1.0 + 0.5 * float(a) for a in range(9)]
C_DBU_input = [0.0 + 0.2 * float(a) for a in range(11)]
design_total = [0] * 99
number = 0
for a in C_SM_input:
    for b in C_MP_input:
        for c in C_DBU_input:
            design_total[number] = [a, a * b, a * c]
            number += 1

'''Main Calculation'''
def multi_processing(n_multiprocessing):
    def main_calculation():
        for nnn in range(99):
            # Please insert your calculation conditions and parameters here

            # Please add a code to solve differential equations. To solve the flow model in the publication, the Crank-Nicolson method was used.

            # The following is an example of the Crank-Nicolson method.
            while t <= t_final:
                # Generate matrix for temperature and concentration   
                Amat_temp = np.zeros((Nl + 2, Nl + 2))
                Amat_SM = np.zeros((Nl + 2, Nl + 2))
                
                bvec_temp = np.zeros((Nl + 2))
                bvec_SM = np.zeros((Nl + 2))
                
                # Replace the values of the matrix depending on calculation conditions
        
                for i in range(1, Nl + 1):
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
                    
                    if C_DBU[i] >= 0.0:
                        bvec_temp[i] = # add equations
                        bvec_SM[i] = # add equations
                    
                    else:
                        bvec_temp[i] = # add equations
                        bvec_SM[i] = # add equations
        
                temp = np.linalg.solve(Amat_temp, bvec_temp)
                C_SM = np.linalg.solve(Amat_SM, bvec_SM)
        
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