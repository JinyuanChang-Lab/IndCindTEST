# -*- coding: utf-8 -*-
"""
The hierarchical neural networks in Section 6.2
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

from torch.utils.data import DataLoader, TensorDataset
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)



  


class AdaptiveNet(nn.Module):
    def __init__(self, 
                 input_features: int, 
                 hidden_features1: int, 
                 hidden_features2: int, 
                 output_features: int):
        """
        Args:
            input_features (int): Size of the input layer (feature dimension).
            hidden_features1 (int): Size of hidden layer 1 (will be partitioned into groups).
            hidden_features2 (int): Size of hidden layer 2 (also the number of groups).
            output_features (int): Size of the output layer.
        """
        super(AdaptiveNet, self).__init__()
        
        # 1) Input layer - Hidden layer 1
        self.fc1 = nn.Linear(input_features, hidden_features1)
        
        # 2) Partition hidden layer 1 into `hidden_features2` groups; each group maps to one unit
        self.hidden_features2 = hidden_features2            # Number of groups / size of hidden layer 2
        self.group_size = hidden_features1 // hidden_features2  # Size of each group (assumes divisibility)
        

        # every linear is (group_size -> 1)
        self.fc2_list = nn.ModuleList([
            nn.Linear(self.group_size, 1) for _ in range(hidden_features2)   # Per-group linear reduction
        ])
        
        # 3) Hidden layer 2 (concatenated group outputs) - Output layer
        self.fc3 = nn.Linear(hidden_features2, output_features)         # Final dense map: hidden_features2 - output_features

    def forward(self, x):
        """
        size of x: [batch_size, input_features]
        """
         
        x = F.sigmoid(self.fc1(x))    # Apply first linear + sigmoid, size [batch_size, hidden_features1] 
        
        # Split along feature dimension into `hidden_features2` chunks  
        # Example: hidden_features1=12, hidden_features2=3 - 3 chunks of size 4
        x_chunks = torch.chunk(x, self.hidden_features2, dim=1)
        
        # Map each chunk through its corresponding small linear layer (group_size - 1)
        outs = []                                                          
        for i, chunk in enumerate(x_chunks):                               
            out_i = F.sigmoid(self.fc2_list[i](chunk))                    # Apply group linear + sigmoid, size [batch_size, hidden_features1] 
            outs.append(out_i)                                            # Append group output
        
        x2 = torch.cat(outs, dim=1)                                       # Concatenate to form hidden layer 2, size [batch_size, hidden_features2] 
        
        out = self.fc3(x2)                                                # Final linear projection, size [batch_size, output_features]  
        return out    



