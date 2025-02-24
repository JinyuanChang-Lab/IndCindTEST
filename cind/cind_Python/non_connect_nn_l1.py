# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 08:51:21 2024

@author: duyue
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

from torch.utils.data import DataLoader, TensorDataset
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


# if torch.cuda.is_available():
#     device = torch.device("cuda:0")  # you can continue going on here, like cuda:1 cuda:2....etc. 
#     print("Running on the GPU")
# else:
#     device = torch.device("cpu")
#     print("Running on the CPU")
    
 
class CustomPartiallyConnectedLayer(nn.Module):
    def __init__(self, hidden_features1, hidden_features2, device):
        super(CustomPartiallyConnectedLayer, self).__init__()
        self.device = device
        self.hidden_features1 = hidden_features1
        self.hidden_features2 = hidden_features2
        
        self.weights = nn.Parameter(torch.randn(self.hidden_features1))
        self.bias = nn.Parameter(torch.randn(self.hidden_features2))
        self.colsize = int(self.hidden_features1 / self.hidden_features2)
        # (p, hid2, hid1)
        self.mask = torch.zeros(self.hidden_features1, self.hidden_features2).to(self.device)
        for j in range (self.hidden_features2):
            self.mask[(j*self.colsize):((j+1)*self.colsize), j] = 1

    def forward(self, x):
        mask1 = self.weights.unsqueeze(-1).expand(-1, self.hidden_features2)
        mask = mask1 * self.mask
        # print(mask.shape, x.shape ,self.bias.shape)
        output = torch.matmul(x, mask) + self.bias
        #print(output)
        return output



 


import torch
import torch.nn as nn
import torch.nn.functional as F

 


class AdaptiveNet(nn.Module):
    def __init__(self, 
                 input_features: int, 
                 hidden_features1: int, 
                 hidden_features2: int, 
                 output_features: int):
        """
        参数:
            input_features:  输入层大小
            hidden_features1: 隐藏层1大小 (被分成 hidden_features2 组)
            hidden_features2: 隐藏层2的神经元数(也是分组数)
            output_features: 输出层大小
        """
        super(AdaptiveNet, self).__init__()
        
        # 1) 输入层 -> 隐藏层1
        self.fc1 = nn.Linear(input_features, hidden_features1)
        
        # 2) 按照 "隐藏层1分成 hidden_features2 组" 的需求:
        #    每组大小 = hidden_features1 // hidden_features2
        #    每一组映射到隐藏层2中的 1 个神经元
        self.hidden_features2 = hidden_features2
        self.group_size = hidden_features1 // hidden_features2  # 要确保能整除
        
        # 用 ModuleList 存储这一批小的全连接
        # 每个 linear 是 (group_size -> 1)
        self.fc2_list = nn.ModuleList([
            nn.Linear(self.group_size, 1) for _ in range(hidden_features2)
        ])
        
        # 3) 隐藏层2 (合并成 hidden_features2 维) -> 输出层
        self.fc3 = nn.Linear(hidden_features2, output_features)

    def forward(self, x):
        """
        x 的形状: [batch_size, input_features]
        """
        # 第一层，简单用 ReLU 作为激活示例
        x = F.sigmoid(self.fc1(x))  # 得到 [batch_size, hidden_features1]
        
        # 在特征维度上分 chunk
        # 比如 hidden_features1=12, hidden_features2=3, 就会分成 3 份, 每份 4 维
        x_chunks = torch.chunk(x, self.hidden_features2, dim=1)
        
        # 分别通过每个小 Linear (group_size -> 1)
        outs = []
        for i, chunk in enumerate(x_chunks):
            # 注意激活函数也可以自己选，或去掉
            out_i = F.sigmoid(self.fc2_list[i](chunk))  # [batch_size, 1]
            outs.append(out_i)
        
        # 把这 hidden_features2 个输出拼起来: [batch_size, hidden_features2]
        x2 = torch.cat(outs, dim=1)
        
        # 最后一层: (hidden_features2 -> output_features)
        out = self.fc3(x2)  # [batch_size, output_features]
        return out

class ThreeLayerPartiallyConnectedNet_l1(nn.Module):
    def __init__(self, input_features, hidden_features1, hidden_features2, output_features, device):
        super(ThreeLayerPartiallyConnectedNet_l1, self).__init__()
        # 第一层是标准的全连接层
        self.hidden_features1 = hidden_features1
        self.hidden_features2 = hidden_features2
        # self.pdim = pdim
        self.input_features = input_features
        self.output_features = output_features

        self.fc1 = nn.Linear(self.input_features, self.hidden_features1)
        # 第二层是自定义的部分连接层
        self.fc2 = CustomPartiallyConnectedLayer(self.hidden_features1, self.hidden_features2, device)
        # 第三层是标准的全连接层
        self.fc3 = nn.Linear(self.hidden_features2, self.output_features)


        
    def forward(self, x):
        # 通过第一层全连接层
        #print(self.fc1weights.shape,x.shape,self.fc1bias.shape)
        # (p, hid1, d) * (B, p, d, 1) -> (B, p, hid1, 1)
        x = torch.sigmoid(self.fc1(x))
        x = torch.sigmoid(self.fc2(x))
        x = self.fc3(x)
        # print(x.shape)
 
        return x
    
"""    
n_batch = 30
input_features = 100
hidden_features1 = 8
hidden_features2 = 2
output_features = 1
pdim =  10
x = torch.randn(n_batch,pdim, input_features,1).to(device)
test1 = ThreeLayerPartiallyConnectedNet_l1(pdim, input_features, hidden_features1, hidden_features2, output_features)
b = test1.forward(x)
"""