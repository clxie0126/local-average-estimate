"""example for semi-varying coefficient model"""

import numpy as np
from numpy import pi, sin, cos, exp, sqrt, std
from laepy import lae
from localreg import *
import matplotlib.pyplot as plt

n, p, d, snr = 500, 10, 2, 20
beta = [3,1.5,2,0,0,0,0,0,0,0]
z = np.random.normal(0, 1,[n,p])
x = np.random.normal(0, 1,[n,d])
u = np.random.uniform(0,1,n)
coef = np.array([sin(60*u), 4*u*np.subtract(1,u)])
coef = coef.T
myrate = np.sum(x*coef, axis=1) + np.dot(z,beta)
sigmae = std(myrate/sqrt(snr))
y = myrate + np.random.normal(0, sigmae,n)


""" 
assemble data and estimate coefficients via local averaging
I is the tuning parameter of local average estimate
"""

semi_data = lae.semivaring_coefficient_data()
semi_data.package(y, z, x, u)
semi_data.sort()

print('linear dim:' + str(semi_data.linear_dim))
print('varying dim:' + str(semi_data.varying_dim))

I = 5
la_design = lae.la_design(semi_data.varying, I)
varying_hat, linear_hat, index_design = lae.full_estimate(semi_data, I)

print('constant coefficients:')
print(linear_hat)

"""
further estimate via local polynomial regression 
"""
target = varying_hat[:,1]
index = index_design
a_hat = localreg(index, target, degree=1, kernel=epanechnikov, width=0.3)

plt.style.use('ggplot')
plt.plot(index, a_hat, label='Local average + Local linear')




"""example for varying coefficient model"""

import numpy as np
from numpy import pi, sin, cos, exp, sqrt, std
from laepy import lae
from localreg import *
import matplotlib.pyplot as plt

n, d, snr = 500, 2, 5
x = np.random.normal(0, 1,[n,d])
u = np.random.uniform(0,1,n)
coef = np.array([sin(60*u), 4*u*np.subtract(1,u)])
coef = coef.T
myrate = np.sum(x*coef, axis=1)
sigmae = std(myrate/sqrt(snr))
y = myrate + np.random.normal(0, sigmae,n)

""" 
assemble data and estimate coefficients
I is the tuning parameter of local average estimate
"""
vary_data = lae.varing_coefficient_data()
vary_data.package(y, x, u)
vary_data.sort()

I = 5
la_design = lae.la_design(vary_data.varying, I)
varying_hat, index_design = lae.varying_estimate(vary_data, I)

"""
further estimate via local polynomial regression 
"""
target = varying_hat[:,1]
index = index_design
a_hat = localreg(index, target, degree=1, kernel=epanechnikov, width=0.3)

plt.style.use('ggplot')
plt.plot(index, a_hat, label='Local average + Local linear')






