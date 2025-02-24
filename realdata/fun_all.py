
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from numpy import random
from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF

import torch


"""###########  radamacher  ###############"""
#def generate_one_minus_one_distribution(size):
#    return torch.tensor([random.choice([-1.0, 1.0]) for _ in range(size)])

"""####################################"""


"""###########  radamacher  ###############"""
"""
def generate_mammen(pv):
    #pv = torch.rand(size)
    reject = (np.sqrt(5)+1)/(2* np.sqrt(5));
    choice1 = -(np.sqrt(5)-1)/2;
    choice2 = (np.sqrt(5)+1)/2;
    M_mm = pv
    for i in range(0,5*100):
        if pv[i] > reject:
            M_mm[i] = choice2
        else:
            M_mm[i] = choice1
    return M_mm
"""
"""####################################"""


def generate_one_minus_one_distribution(a, b):
    return 2 * torch.randint(0, 2, (a, b)) - 1
    # pv = torch.rand(a, b)
    # reject = 0.5
    # choice1 = -1
    # choice2 = 1

    # mask = pv > reject
    # M_red = np.where(mask, choice2, choice1)

    # return torch.tensor(M_red)



def generate_mammen(a, b):
    pv = torch.rand(a, b)
    reject = (np.sqrt(5) + 1) / (2 * np.sqrt(5))
    choice1 = -(np.sqrt(5) - 1) / 2
    choice2 = (np.sqrt(5) + 1) / 2

    mask = pv > reject
    M_mm = np.where(mask, choice2, choice1)

    return torch.tensor(M_mm)



###########  my-empirical  ############### 
"""
def my_empirical1(x_1, x_2):
    n1 = x_1.shape[0]
    n2 = x_2.shape[0]
    y_1 = torch.zeros((n2,1))
    for i in range(0,n2):
        y_2 = torch.zeros((n1,1))
        for j in range(0,n1):
            y_2[j,0] = x_1[j] <= x_2[i]
        if torch.mean(y_2) <= 1/n1:
            y_1[i,0] = 1/n1
        elif torch.mean(y_2) >= (1-1/n1):
            y_1[i,0] = (1-1/n1)
        else:
            y_1[i,0] = torch.mean(y_2)
    return y_1

"""




def my_empirical(x_1, x_2):
    n1 = x_1.shape[0]
    n2 = x_2.shape[0]

    y_1 = torch.empty(n2, 1, dtype=x_1.dtype, device=x_1.device)

    for i in range(n2):
        y_2 = x_1.view(-1, 1) <= x_2[i]  # Vectorized comparison
        mean_y_2 = y_2.float().mean()

        y_1[i, 0] = torch.clamp(mean_y_2, 1/n1, (1-1/n1))

    return y_1


"""
import time


x_1  = torch.rand(50,100)
x_2 = torch.rand(50,100)


st =time.time()
a= data_transform(x_1, x_2)
et = time.time()
print(et-st)
st =time.time()
b= my_empirical1(x_1, x_2)
et = time.time()
print(et-st)
a-b
"""


"""###########  data-transform  ###############"""
def data_transform(x_1,x_2):
    n_1 = x_1.shape[0]
    dim_x = x_2.shape
    pp = dim_x[1]
    y_1 = torch.zeros(dim_x)
    ecdf = [ECDF(x_1[:, i]) for i in range(x_1.shape[1])]
    for j in range(0,pp):
        y_1[:, j] = torch.tensor(ecdf[j](x_2[:, j]), dtype=torch.float32)
        y_1[:, j] = torch.clamp(y_1[:, j], 1/n_1, (1-1/n_1))
        # y_1[:, j]= my_empirical(x_1[:, j], x_2[:, j]).T
        y_1[:, j] = torch.tensor(norm.ppf(y_1[:,j], 0, 1))
    return y_1
    
"""####################################"""
"""
def data_transform1(x_1,x_2):
    dim_x = x_2.shape
    nn = dim_x[0]
    pp = dim_x[1]
    y_1 = torch.zeros(dim_x)
    y_2 = torch.zeros(dim_x)
    for j in range(0,pp):
        y_1[:, j]= my_empirical(x_1[:, j], x_2[:, j]).T
        for i in range(0, nn):
            y_2[i,j] = norm.ppf(y_1[i,j], 0, 1)
    return y_2

"""













