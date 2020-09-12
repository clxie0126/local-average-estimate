#__author__ = 'Chuanlong'

import numpy as np
from scipy.linalg import block_diag
from sklearn.linear_model import LassoLarsCV
from numpy import array, pi, sin, cos, sqrt, var, log, exp, dot

"""
Code for local average estimator
response (n by 1) is the response of a semevarying coefficient models. 
linear (n by p) is the predictors of the linear part of a semivarying coefficent model.
varying (n by d) is the predictors of the varying linear part of a semivarying coefficient model.
index (n by 1) is the index variable.
the semi-varying coefficient model is  response = a(index).T@varying + beta.T@linear + vep
the varying coefficient model is  response = a(index).T@varying + vep
a(.) is a function of index with d outputs
# b is a p-length vector, which is of interest here.
the response, linear and varying are sorted according to index.
"""

class semivaring_coefficient_data:
    def package(self, response, linear, varying, index):
        if response.shape[0] == linear.shape[0] and response.shape[0] == varying.shape[0] and response.shape[0] == index.shape[0]:
            self.num_sample = response.shape[0]
            self.linear_dim = linear.shape[1]
            self.varying_dim = varying.shape[1]
            self.response = response
            self.linear = linear
            self.varying = varying
            self.index = index
        else:
            print('The sample sizes are not consistent')            
    def sort(self):
        self.rank = np.argsort(self.index)
        self.response = self.response[self.rank]
        self.linear = self.linear[self.rank]
        self.varying = self.varying[self.rank]
        self.index = self.index[self.rank]

class varing_coefficient_data:
    def package(self, response, varying, index):
        if  response.shape[0] == varying.shape[0] and response.shape[0] == index.shape[0]:
            self.num_sample = response.shape[0]
            self.varying_dim = varying.shape[1]
            self.response = response
            self.varying = varying
            self.index = index
        else:
            print('The sample sizes are not consistent')
    def sort(self):
        self.rank = np.argsort(self.index)
        self.response = self.response[self.rank]
        self.varying = self.varying[self.rank]
        self.index = self.index[self.rank]
                
def la_design(varying, k):
    local_average_design = block_diag(*np.vsplit(varying, k))
    return(local_average_design)

def hat_matrix(design):
    hat_matrix = dot(dot(design, np.linalg.pinv(dot(design.T, design))), design.T )
    return(hat_matrix)

def full_estimate(semi_data, I):
    """estimate the varying and linear coefficients of a semi-varying coefficient model"""
    n = semi_data.num_sample
    k = int(np.floor(n/I))
    varying_design = la_design(semi_data.varying, k)
    full_design = np.hstack((varying_design, semi_data.linear))
    full_hat = np.linalg.pinv(full_design.T @ full_design) @ (full_design.T @ semi_data.response)
    varying_hat = full_hat[:(k*semi_data.varying_dim)]
    varying_hat = np.reshape(varying_hat, (k,semi_data.varying_dim))
    linear_hat = full_hat[(k*semi_data.varying_dim):]
    index_design = np.mean( np.reshape(semi_data.index, (k,I)), axis=1)
    return(varying_hat, linear_hat, index_design)
        

def varying_estimate(varying_data, I):
    """estimate the varying coefficients of a varying coefficient model"""
    n = varying_data.num_sample
    k = int(np.floor(n/I))
    design = la_design(varying_data.varying, k)
    varying_hat = np.linalg.pinv(design.T @ design) @ (design.T @ varying_data.response)
    varying_hat = np.reshape(varying_hat, (k,varying_data.varying_dim))
    index_design = np.mean( np.reshape(varying_data.index, (k,I)), axis=1)
    return(varying_hat,index_design)
           