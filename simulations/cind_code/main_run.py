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
from pcind_functions_dep.Cind_gaussian_fun import *



if torch.cuda.is_available():
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

"""Initialize the random number generator with a fixed seed."""
def seed_torch(seed=42):
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
def init_logging_for_sim(sim_name):

    logger_name = f"runner.{sim_name}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.propagate = False   

    if not logger.handlers:
        fmt = logging.Formatter('%(asctime)s - %(message)s')
        fh = logging.FileHandler(f"{sim_name}_nn.log", mode='a', encoding='utf-8')
        fh.setLevel(logging.INFO); fh.setFormatter(fmt)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO); ch.setFormatter(fmt)
        logger.addHandler(fh); logger.addHandler(ch)

    return logger

def close_logger_handlers(logger: logging.Logger):
    for h in list(logger.handlers):
        try:
            h.flush(); h.close()
        except Exception:
            pass
        logger.removeHandler(h)





# --- Main function ---
def run_simulation(i, n, p, q, m, alpha, gene_data, dep, para_time, seed_base=12345):
    seed_torch((i + 1000) * seed_base)
 
    data_my = gene_data(n, p, q, m, dep)
    x = data_my["x"]; y = data_my["y"]; z = data_my["z"]


    """ Configure parallelization parameters (e.g., use 8 GPUs and 100 CPU cores)."""
    para_time = 100
    tmp_device = 'cuda'
    tasks_per_machine = para_time // 8
    tmp = tasks_per_machine * 8 
    adjusted_index = (i - tmp if i > (tmp-1) else i) // tasks_per_machine
    machine_index = adjusted_index % 8
    device = tmp_device +':'+str(machine_index)


    # inplement the proposed method with selecting n_3^{opt}
    order = ['Gaussian', 'Mammen', 'Rademacher']
    batchsize = max(1, math.ceil((n / 100) * 16))
    res_orig_1 =  Cind_Gtest_py(device, "S-selectn3", 'all', x, y, z, alpha, batchsize, 
                                              hidden_features1 = 128, hidden_features2 = 32,  lr=0.01,n_epochs= 400, patience=30) 
    res_orig = [int(res_orig_1[k]['reject']) for k in order]
    

    col1 = torch.tensor(res_orig)                      # 3-by-1 vector
    
    return col1

  

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

####

pairs = [('my_data_ex6', int(0)), ('my_data_ex6', int(p/10)), ('my_data_ex6', int(p/5)),
         ('my_data_ex7', int(0)), ('my_data_ex7', int(7)), ('my_data_ex7', int(8)),
         ('my_data_ex8', int(0)), ('my_data_ex8', int(p/10)), ('my_data_ex8', int(p/5)),
         ('my_data_ex9', int(0)), ('my_data_ex9', int(p/10)), ('my_data_ex9', int(p/5)),
         ('my_data_ex10', int(0)), ('my_data_ex10', int(p/10)), ('my_data_ex10', int(p/5))]
#-----------------------------------------------------------------

cores = 100  # the number of CPU cores

# conduct experiments using parallel processing
loggers_by_sim = {}         
try:
    for sim, dep in pairs:
        # create log file
        logger = loggers_by_sim.get(sim)
        if logger is None:
            logger = init_logging_for_sim(sim)
            loggers_by_sim[sim] = logger

        try:
            with Parallel(n_jobs=cores, **PARALLEL_KW) as parallel:
                tmp_results = parallel(
                    delayed(run_simulation)(
                        i, n, p, q, m, alpha, getattr(gene_data_fun, sim), dep, cores, seed_base=12345
                    )
                    for i in tqdm(range(n_my), desc=f"{sim}-dep{dep}")
                )
            # save the results
            out_dir = os.path.join('.', f'results_{sim}')  # create a folder
            os.makedirs(out_dir, exist_ok=True)

            mats = [to_numpy(m) for m in tmp_results]
            all_runs = np.vstack([to_numpy(M).reshape(3,) for M in mats])   
            means = np.nanmean(all_runs, axis=0)

            row_names = ['Gaussian', 'Mammen', 'Rademacher']
            df_means  = pd.DataFrame(means, index=row_names)
            print(df_means.to_string())
            logger.info("dep=%s Column means:\n%s", dep, df_means.to_string())

            # save to a csv file
            df_means.to_csv(os.path.join(out_dir, f"{sim}_dep_{dep}.csv"))

            # cleanup
            del tmp_results, mats, row_names, means, df_means
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
                torch.cuda.synchronize()

        except Exception:
            error_message = f"code is wrong:\n{traceback.format_exc()}"
            print(error_message)
            logger.error(error_message)

finally:
    for lg in loggers_by_sim.values():
        close_logger_handlers(lg)




         
