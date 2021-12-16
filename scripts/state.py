import numpy as np
import copy
import itertools
from scipy import linalg as sLA
from numpy import linalg as LA

from Gates import *
from functions import *

class Qpsi():
    def __init__(self, N = 2):
        self.N = N
        self.coefs = None
    def build_random_state(self):
        self.coefs = np.random.rand(2**self.N) + 1j*np.random.rand(2**self.N)
        self.coefs = self.coefs/np.sqrt(np.sum(self.coefs*self.coefs.conj()))
        
    def build_zero_state(self):
        coefs = np.zeros(2**self.N) + 1j*np.zeros(2**self.N)
        coefs[0] = 1
        self.coefs = coefs
        
    def set_coefs(self, coefs):
        if len(coefs) == 2**self.N:
            self.coefs = np.array(coefs)
            
    def get_coefs(self):
        return self.coefs
    
    def apply_U(self, U, axis):
        self.coefs = apply_U(U, self.coefs.copy(), axis)
        
    def apply_long_Toffoli(self, axis):
        self.coefs = apply_long_Toffoli(self.coefs.copy(), axis = axis, N = len(axis)-1)
        
        
        
class Qstate():
    def __init__(self, N = 2):
        self.A_dims = []
        self.B_dims = []
        self.tensor = np.zeros(tuple([2]*N))
        self.vec = np.zeros(2**N)
        self.N = N
        self.rho = None
        
        
    def build_pure_random_state(self): # строит случайное состояние заданной размерности
#         vec = np.random.randn(*([2]*self.N)) + 1j*np.random.randn(*([2]*self.N))
#         self.tensor = vec

        vec = np.random.randn(2**self.N) + 1j*np.random.randn(2**self.N)
        self.vec = vec/np.sqrt(np.sum(vec*vec.conj()))
#         self.tensor = np.reshape(self.vec, [2]*self.N)
        self.rho = np.kron(np.reshape(self.vec, (self.vec.shape[0], 1)), self.vec.conj())
    
    def build_state(self, vec):
        self.vec = vec
        #         self.tensor = np.reshape(self.vec, [2]*self.N)
        self.rho = np.kron(np.reshape(self.vec, (self.vec.shape[0], 1)), self.vec.conj())
        
        
    def apply_U(self, U):
        self.vec = np.dot(U,self.vec)
#         self.tensor = np.reshape(self.vec, [2]*self.N)
        self.rho = np.kron(np.reshape(self.vec, (self.vec.shape[0], 1)), self.vec.conj())
        
    def get_tensor(self): # возвращаем вектор
        return self.tensor
    
    def get_vec(self): # возвращаем вектор
        return self.vec
    
    def get_rho(self): # возвращаем матрицу плотности
        return self.rho
      
        
        
