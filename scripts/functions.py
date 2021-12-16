import numpy as np
import copy
import itertools
from scipy import linalg as sLA
from numpy import linalg as LA

from Gates import *


# возвращает все комбинации элементов листа без учёта перестановок
def get_combinations(val):
    comb_list = []
    for j in range(1,len(val)+1):
        com_set = itertools.combinations(val, j)
        for i in com_set: 
            comb_list.append(i)
    return comb_list


# применяет унитарное преобразование U над кубитами под индексами axis в состоянии с коэффициентами coefs.
# U - матрица унитарного преобразования
# coefs - вектор состояния (коээфициенты)
# axis - индексы выбранных кубитов

def apply_U(U, coefs, axis):
    axis = list(reversed(axis))
    if ((len(U) != 2**len(axis)) | (len(coefs) < len(U)) | (2**max(axis) > len(coefs))):
        print('Некорректно задано преобразование')
        return 0
    a = coefs.copy()
    b = [0]*len(a)
    
    ax = axis[0]
    
    N = len(U)  # размерность преобразования 
    K = len(coefs) # размерность состояния системы
    
    axis_comb = get_combinations(axis)
    old_ind_set = {K+1}
    
    for i in range(K):
        
        zero_index = (K-1)
            
        for ax in axis:
            zero_index = zero_index^(1<<ax)
            
        zero_index = i&zero_index
        
        if zero_index in old_ind_set:
            continue
            
        old_ind_set.add(zero_index)

                
        index = zero_index
        
        b[index] += a[zero_index]*U[0][0]
        
        
        u1 = 0
        for ax_list in axis_comb:
            u1 += 1
            m = 0
            for ax in ax_list:
                m += 1<<ax

            b[index] += a[zero_index^m]*U[0][u1]
            

        u0 = 0
        for ax_list in axis_comb:
            u0 += 1
            r = 0
            for ax in ax_list:
                r += 1<<ax
            
            index = zero_index^r
            old_ind_set.add(index)
            
            b[index] += a[zero_index]*U[u0][0]
            
            
            
            u1 = 0
            for ax1_list in axis_comb:
                u1 += 1
                m = 0
                for ax1 in ax1_list:
                    m += 1<<ax1
                b[index] += a[zero_index^m]*U[u0][u1]
            
        
    return b


def apply_Toffoli(coefs, axis):
    coefs = apply_U(U = H(), coefs = coefs, axis = [axis[2]])
    coefs = apply_U(U = CX(), coefs = coefs, axis = [axis[1],axis[2]])
    
    coefs = apply_U(U = T().conj(), coefs = coefs, axis = [axis[2]])
    coefs = apply_U(U = CX(), coefs = coefs, axis = [axis[0],axis[2]])
    
    coefs = apply_U(U = T(), coefs = coefs, axis = [axis[2]])
    coefs = apply_U(U = CX(), coefs = coefs, axis = [axis[1],axis[2]])
    
    coefs = apply_U(U = T().conj(), coefs = coefs, axis = [axis[2]])
    coefs = apply_U(U = CX(), coefs = coefs, axis = [axis[0],axis[2]])
    
    coefs = apply_U(U = T(), coefs = coefs, axis = [axis[1]])
    coefs = apply_U(U = T(), coefs = coefs, axis = [axis[2]])
    
    coefs = apply_U(U = CX(), coefs = coefs, axis = [axis[0],axis[1]])
    coefs = apply_U(U = H(), coefs = coefs, axis = [axis[2]])
    
    coefs = apply_U(U = T(), coefs = coefs, axis = [axis[0]])
    coefs = apply_U(U = T().conj(), coefs = coefs, axis = [axis[1]])
    
    coefs = apply_U(U = CX(), coefs = coefs, axis = [axis[0],axis[1]])
    
    return coefs

    
    
def expand_state(coefs, n_qubits):
    
    c = coefs.copy()
    exp_coefs = np.zeros(len(c)*2**n_qubits)
    exp_coefs[:len(c)] = c
    
    return exp_coefs

def reduce_state(coefs, N):
    c = coefs.copy()
    return c[:2**N]   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def apply_long_Toffoli(coefs = None, axis = None, N = 2):
    
    if N == 1:
        new_coefs = apply_U(CX(), coefs.copy(), axis)
    
    elif N == 2:
        new_coefs = apply_Toffoli(coefs.copy(), axis)
        
    else:
    
        anc_num = N-1
        tq_num = len(axis) - 1
        vq_num = int(np.log2(len(coefs)))
        
        
        
        
        c = expand_state(coefs, anc_num )


        c = apply_Toffoli(coefs = c, axis = [axis[0], axis[1], vq_num])

        
        for i in range(tq_num-2):
            c = apply_Toffoli(coefs = c, axis = [axis[i + 2], vq_num + i, vq_num + i + 1])


        c = apply_U(CX(), coefs = c, axis = [vq_num + anc_num - 1, axis[-1]])     
                  

        for i in reversed(range(tq_num-2)):
            c = apply_Toffoli(coefs = c, axis = [axis[i + 2], vq_num + i, vq_num + i + 1])

        c = apply_Toffoli(coefs = c, axis = [axis[0], axis[1], vq_num])

        new_coefs = reduce_state(coefs = c, N = vq_num)

    
    
    return new_coefs
    
   
        
    


    
                         
    
    
    
        
