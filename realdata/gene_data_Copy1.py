# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 16:20:43 2024

@author: duyue
"""

import torch
import numpy as np
import itertools


"""#########################  ex6 ######################"""

def my_data_ex6(n,p,q,m,dep):
    # 创建StudentT分布实例
    t_distribution = torch.distributions.StudentT(2)
    x = t_distribution.sample((n,p))
    y = t_distribution.sample((n,q))
    z = t_distribution.sample((n,m))
    w = t_distribution.sample((n,p))
    v = t_distribution.sample((n,q))

    

    # 获取所有可能的两两列组合的索引
    combinations = list(itertools.combinations(range(m), 2))

    # 计算组合的数量
    num_combinations = len(combinations)

    # 对每一对组合的对应元素进行相乘，并将结果存储在x矩阵中
    for i, (col1, col2) in enumerate(combinations):
        x[:, i] = z[:, col1] * z[:, col2]
        # x[:, (i+num_combinations)] = torch.sin(z[:,col1]) + torch.cos(z[:,col2]) + z[:,col1]**2 + z[:,col2]**2
        y[:, i] = z[:, col1] + z[:, col2] 
    if dep!=0:
        for i in range(dep):
            t_distribution2 = torch.distributions.StudentT(1)
            eps = t_distribution2.sample((1,n))
            x[:, i] = x[:, i] + eps + 3* eps**3  
            y[:, i] = y[:, i] + eps + 3* eps**3  
         
    #save as an array
    dat_my = [
        ("x", x),
        ("y", y),
        ("z", z)
    ]       

    data_label = {label: matrix for label, matrix in dat_my}
    return  data_label




# def my_data_ex6(n,p,q,m,dep):
#     # 创建StudentT分布实例
#     t_distribution = torch.distributions.StudentT(1)
#     x = t_distribution.sample((n,p))
#     y = t_distribution.sample((n,q))
#     zz = t_distribution.sample((n,m))
#     mm = int(p/20)

#     x[:, 0:4] = zz[:, 0:4]
#     x[:, 4] = zz[:,3] + zz[:,4]
#     y[:, 0] = torch.sum(zz[:, 0:4], dim=1)

#     # for i in range(1,mm+1):
#     #     index1 = (i-1)*6 + 1
#     #     index2 = i*6
#     #     x[:, (index1-1):(index2-2)] = zz[:, (index1-1):(index2-2)]
#     #     x[:, (index2-2)] = zz[:,(index2-2)] + zz[:,(index2-1)]
#     #     y[:, (i-1)] = torch.sum(zz[:, (index1-1):index2], dim=1)
            
#     if dep == 0:
#         z = zz
#     elif dep == 1:
#         z= t_distribution.sample((n,m))
#         # indexm =  np.arange(mm, mm*6, 6)
#         z[:, 0] = zz[:, 1]
#     else:
#          z= t_distribution.sample((n,m))
#          z[:, 0:4] = zz[:, 0:4]
         
#     #save as an array
#     dat_my = [
#         ("x", x),
#         ("y", y),
#         ("z", z)
#     ]       

#     data_label = {label: matrix for label, matrix in dat_my}
#     return  data_label





    


"""######################  ex7 #####################"""
def my_data_ex7(n,p,q,m,dep):
     x = torch.randn((n,p))
     y = torch.randn((n,q))
     z = torch.rand((n,m)) * 2 - 1
     rho =dep * 0.1

     for i in range(m):
         eps_y1 = torch.rand((n,48)) * 0.5 - 0.25
         eps_y = torch.sum(eps_y1, dim =1)
         y[:,i] = z[:,i] + 0.25 * z[:, i]**2 + eps_y
         beta = rho / (2*np.sqrt(1 - rho**2))
         eps_x = torch.randn(1,n)
         x[:, i] = 5*beta * y[:, i] + z[:, i] + eps_x
     for i in range(int(p-m)):
         x[:, (i+m)] = x[:, (i+m)] + 5*beta * y[:, (i+m)] + torch.randn((1,n))
     #save as an array  
     dat_my = [
         ("x", x),
         ("y", y),
         ("z", z)
     ]       

     data_label = {label: matrix for label, matrix in dat_my}
     return  data_label

 
 
    
 
"""######################  ex8 #####################"""
def my_data_ex8(n,p,q,m,dep):
     t_distribution = torch.distributions.StudentT(1)
     w = torch.randn((n,p))
     v = torch.randn((n,q))
     z = torch.randn((n,m)) 
     x = torch.randn(n,p)
     y = torch.randn(n,q)

     for i in range(m):
         a = 0.7 *( z[: , i]** 3 / 5 + z[:, i]/2) + torch.tanh(w[:,i])
         b = (z[:, i]**3 /4 + z[:, i]) / 3 + v[:, i]
         x[:,i] = a + a**3/3 + torch.tanh(a/3)/2
         y[:, i] = (b + torch.tanh(b/3)) **3
     if dep!=0:
         for i in range(dep):
             eps = 3*t_distribution.sample((1,n))
             x[:, i] = x[:, i] + eps
             y[:, i] = y[:, i] + eps
     #save as an array 
     dat_my = [
         ("x", x),
         ("y", y),
         ("z", z)
     ]       

     data_label = {label: matrix for label, matrix in dat_my}
     return  data_label

 



"""######################  ex9 #####################"""
def my_data_ex9(n,p,q,m,dep):
     t_distribution = torch.distributions.StudentT(1)
     w = torch.randn((n,p))
     v = torch.randn((n,q))
     z = torch.randn((n,m)) 
     x = torch.randn(n,p)
     y = torch.randn(n,q)

     for i in range(m):
         a = 0.5 *( z[: , i]** 3 / 7 + z[:, i]/2)
         b = (z[:, i]**3 /2 + z[:, i]) / 3 
         a1 = a + torch.tanh(w[:, i])
         x[:,i] = a1 + a1**3/3  
         b1 = b + v[:, i]
         y[:, i] = b1 + torch.tanh(b1/3)
     if dep!=0:
         for i in range(dep):
             eps = t_distribution.sample((1,n))
             x[:, i] = 0.5 * x[:, i] + 3 * torch.cosh(eps) 
             y[:, i] = 0.5 * y[:, i] + 3 *  torch.cosh(eps**2) 
     "save as an array"   
     dat_my = [
         ("x", x),
         ("y", y),
         ("z", z)
     ]       

     data_label = {label: matrix for label, matrix in dat_my}
     return  data_label




def my_data_ex10(n,p,q,m,dep):
     w = torch.randn((n,p))
     v = torch.randn((n,q))
     z = torch.randn((n,m)) 
     x = torch.randn(n,p)
     y = torch.randn(n,q)
     
     if dep ==0:
         for i in range(int(p/4)):
             zz = torch.mean(z[:, 0:m], dim=1)
             x[:,i] = torch.tanh(zz + w[:, i])
             y[:, i] = (zz+v[:,i])**3
     else:
         for i in range(int(p/4)):
             zz = torch.mean(z[:, 0:m], dim=1)
             x[:,i] = torch.tanh(zz + w[:, i])
             y[:, i] = (zz+v[:,i])**3
         for i in range(dep):
             eps = 3* torch.randn(1,n)
             zz = torch.mean(z[:, 0:m], dim=1)
             x[:,i] = torch.tanh(zz + w[:, i] + eps)
             y[:, i] = (zz+v[:,i] + eps)**3
     "save as an array"   
     dat_my = [
         ("x", x),
         ("y", y),
         ("z", z)
     ]       

     data_label = {label: matrix for label, matrix in dat_my}
     return  data_label




 
















