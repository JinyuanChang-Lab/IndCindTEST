"""
Notes:
    - The number of cores for parallelization can be adjusted by the "cores" parameter.
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"  
os.environ["MKL_NUM_THREADS"] = "1"   
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_WAIT_POLICY"] = "PASSIVE"
os.environ["KMP_INIT_AT_FORK"] = "FALSE"

os.environ["TORCH_NUM_THREADS"] = "1"            
os.environ["TORCH_NUM_INTEROP_THREADS"] = "1" 

import logging

import multiprocessing as mp
try:
    mp.set_start_method("spawn", force=True)   #  fork  
except RuntimeError:
    pass
import numpy as np
from numpy import random
from joblib import Parallel, delayed
PARALLEL_KW = dict(
    backend="loky",
    temp_folder="./tmp/joblib",
    max_nbytes=None,
    mmap_mode=None,
    # inner_max_num_threads=1,   #   joblib>=1.3
)
os.makedirs("./tmp/joblib", exist_ok=True)

import torch
torch.set_num_threads(1)  #  

from tqdm import tqdm
import gc
import math

import pandas as pd
import traceback

""" the proposed functions """
import gene_data_fun
import Cind_gaussian_fun



if torch.cuda.is_available():
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

"""Initialize the random number generator with a fixed seed."""
def seed_torch(seed=42):
    """For reproducibility"""
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)


""" Select GPU if available; otherwise use CPU and cap BLAS/Torch threads to 1 to avoid oversubscription. """
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

print(device)


"""
Define the log function
"""
def init_logging(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
     
    file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    file_handler.setFormatter(formatter)
    
     
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    # 
    
    if not logger.handlers:
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)
    
    return logger



"""
Apply the same random permutation to shuffle x, y, z along the first dimension.
Returns: x_shuf, y_shuf, z_shuf, perm
"""
def same_perm_shuffle(x, y, z, generator=None):

    assert x.size(0) == y.size(0) == z.size(0), "x,y,z must be same"
    perm = torch.randperm(x.size(0), generator=generator, device=x.device)
    return x[perm], y[perm], z[perm], perm



"""
Main function: 
    - perform the proposed conditional independence test based on nonparametric regressions with sample splitting and without sample splitting
    - We perform ten different sample permutations to assess whether randomness in the splitting procedure affects performance. 
"""
def run_simulation(i, n, p, q, m, alpha, gene_data, para_time, t=9, seed_base=12345):
    seed_torch((i + 1000) * seed_base)
 
    data_my = gene_data(n, p, q, m, dep)
    x = data_my["x"]; y = data_my["y"]; z = data_my["z"]

    """ 
    # Configure parallelization parameters (e.g., use 8 GPUs and 100 CPU cores).
    para_time = 100
    tmp_device = 'cuda'
    tasks_per_machine = para_time // 8
    tmp = tasks_per_machine * 8 
    adjusted_index = (i - tmp if i > (tmp-1) else i) // tasks_per_machine
    machine_index = adjusted_index % 8
    device = tmp_device +':'+str(machine_index)
    """
   
    order = ['Gaussian', 'Mammen', 'Rademacher']
    #  A) Full-sample result (setting = "full-sample")
    batchsize = max(1, math.ceil((n / 100) * 32))
    res_all_1 = Cind_gaussian_fun.Cind_Gtest_py(device, "full-sample", 'all', x, y, z, alpha, batchsize, 
                                              hidden_features1 = 128, hidden_features2 = 32,  lr=0.01,n_epochs= 400, patience=30) 
    res_all = [int(res_all_1[k]['reject']) for k in order]
    
    #  B) Original order result (setting = "S-selectn3")
    batchsize = max(1, math.ceil((n / 100) * 16))
    res_orig_1 = Cind_gaussian_fun.Cind_Gtest_py(device, "S-selectn3", 'all', x, y, z, alpha, batchsize, 
                                              hidden_features1 = 128, hidden_features2 = 32,  lr=0.01,n_epochs= 400, patience=30) 
    res_orig = [int(res_orig_1[k]['reject']) for k in order]
    
    #  C) Shuffle with the same Generator for t rounds (different perm each time), save separately (setting = "S-selectn3")
    g = torch.Generator(device=x.device).manual_seed((1 + 200) * seed_base)  # set seed
    shuffle_results = []   
    perms = []            
    for kk in range(t):
        xb, yb, zb, perm_k = same_perm_shuffle(x, y, z, generator=g)  # apply the same perm to x,y,z
        rb_1 = Cind_gaussian_fun.Cind_Gtest_py(device, "S-selectn3", 'all', xb, yb, zb, alpha, batchsize, 
                                              hidden_features1 = 128, hidden_features2 = 32,  lr=0.01,n_epochs= 400, patience=30)
        rb =  [int(rb_1[k]['reject']) for k in order]
        shuffle_results.append(rb)  # collect result of this shuffle
        perms.append(perm_k)        # collect permutation used in this shuffle


    col1 = torch.tensor(res_all)                      # full-sample results: 3-by-1 vector
    col2 = torch.tensor(res_orig)                     # sample-splitting results: 3-by-1 vector
    cols_shuf = torch.tensor(shuffle_results).T       # Sample-splitting results over 9 permutations: 3-by-t matrix with t=9

    M_res = torch.column_stack([col1, col2, cols_shuf]) # concatenate all results
     
    return M_res



"""
Convert to a unified NumPy array
"""
def to_numpy(x):
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().numpy()
    return np.asarray(x)


#--------------------- simulation setting --------------------- 
alpha=0.05
p=100
q=100
m=5
n = 100 
n_my = 1000
#### random splitting number 
ave = 9
sim_list = {'my_data_ex6': 0}   

#-----------------------------------------------------------------

cores = 170  # the number of CPU cores

# conduct experiments using parallel processing
for sim, dep in sim_list.items():
    logger = init_logging(f"log_{sim}_nn.log")
    try:
        with Parallel(n_jobs=cores, **PARALLEL_KW) as parallel:
            tmp_results = parallel(
                delayed(run_simulation)(
                    i, n, p, q, m, alpha, getattr(gene_data_fun, sim), cores, ave, seed_base=12345
                )
                for i in tqdm(range(n_my))
            )

        # save the results
        # tmp_results: a list of 1000 3-by-K matrices (K=2+ave)
        mats = [to_numpy(m) for m in tmp_results]     # convert to a NumPy array
        n_my = len(mats) 
        K = mats[0].shape[1]
        assert all(M.shape == (3, K) for M in mats)
        
        # row names
        row_names = ['Gaussian', 'Mammen', 'Rademacher']
        # column names
        col_names = ['Full','Orig'] + [f'Shuf{i+1}' for i in range(K-2)]
        
        out_dir = os.path.join('.', f'results_{sim}')    # create a folder
        os.makedirs(out_dir, exist_ok=True)
        
        means_rows = []
        
        for r, rname in enumerate(row_names):
            data = np.vstack([M[r, :] for M in mats])  # n_my-by-K matrix
            df = pd.DataFrame(data, columns=col_names, index=[f'run_{i+1}' for i in range(n_my)])
            df.to_csv(os.path.join(out_dir, f"{sim}_{rname}.csv"))
        
            row_mean = data.mean(axis=0, dtype=float)  # calculate the means of each column of all_runs
            means_rows.append(row_mean)                # 1-by-K vector    
        
        means_mat = np.vstack(means_rows)               # concatenate all results: 3-by-K matrix
        # save to a csv file
        df_means  = pd.DataFrame(means_mat, index=row_names, columns=col_names)
        print(df_means.to_string())
        logger.info("Column means:\n%s", df_means.to_string())
        df_means.to_csv(os.path.join(out_dir, f"{sim}_ALL_col_means.csv"))
        
        # cleanup
        del tmp_results, mats, means_rows, means_mat, df_means
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            torch.cuda.synchronize()

    except Exception:
        error_message = f"code is wrong:\n{traceback.format_exc()}"
        print(error_message)
         
