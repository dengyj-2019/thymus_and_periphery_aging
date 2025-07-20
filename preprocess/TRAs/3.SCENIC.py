#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import glob


f_db_glob = "/data1/00.software/01.database/SCENIC/databases/hg19*species*feather"
f_db_names = ' '.join( glob.glob(f_db_glob) )
f_motif_path = "/data1/00.software/01.database/SCENIC/resources/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
f_tfs = "/data1/00.software/01.database/SCENIC/resources/hs_hgnc_tfs.txt"

import sys


# In[2]:


get_ipython().system(" /home/dengyj/.conda/envs/pyscenic/bin/pyscenic grn {'filtered_TRA_TEC.loom'} {f_tfs} -o filtered_TRA_TEC_adj.csv --num_workers 16 --seed 1234")


# In[3]:


get_ipython().system("/home/dengyj/.conda/envs/pyscenic/bin/pyscenic ctx filtered_TRA_TEC_adj.csv     {f_db_names}     --annotations_fname {f_motif_path}     --expression_mtx_fname {'filtered_TRA_TEC.loom'}     --output filtered_TRA_TEC_reg.csv     --mask_dropouts     --num_workers 16     --mode custom_multiprocessing")


# In[4]:


get_ipython().system("/home/dengyj/.conda/envs/pyscenic/bin/pyscenic aucell     {'filtered_TRA_TEC.loom'}     filtered_TRA_TEC_reg.csv     --output {'filtered_TRA_TEC_auc.loom'}     --num_workers 16 --seed 1234")


# In[ ]:




