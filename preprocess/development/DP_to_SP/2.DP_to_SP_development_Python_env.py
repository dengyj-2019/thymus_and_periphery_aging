#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os,sys
sys.path.remove('/home/dengyj/.local/lib/python3.8/site-packages')
import numpy as np 
import scanpy as sc
import pandas as pd
import loompy as lp


# In[ ]:





# In[3]:


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300, facecolor='white')


# In[4]:


Dev_1 = sc.read_loom('Dev_1.loom')
Dev_1_cp = Dev_1
Dev_1_metadata = pd.read_csv('Dev_1_metadata.csv', index_col= 0)
Dev_1_harmony = np.array(pd.read_csv('Dev_1_harmony.csv', index_col= 0))
Dev_1.obs = Dev_1_metadata
Dev_1.obsm['X_pca'] = Dev_1_harmony

sc.pp.neighbors(Dev_1, n_pcs=20)
sc.tl.diffmap(Dev_1)
Dev_1.uns['iroot'] = np.flatnonzero(Dev_1.obs['identity']  == 'DP_P')[0]
sc.tl.dpt(Dev_1)
sc.tl.paga(Dev_1, groups='identity_2')


# In[5]:


sc.pl.paga(Dev_1, color = ['identity_2'])
sc.pl.paga(Dev_1, color=['dpt_pseudotime'])

sc.pl.diffmap(Dev_1, color = 'identity_2')
sc.pl.diffmap(Dev_1, color = 'dpt_pseudotime')


# In[6]:


sc.tl.draw_graph(Dev_1, maxiter = 100, layout = 'fa',init_pos='paga')
sc.pl.draw_graph(Dev_1, color = ['identity_2'])
sc.pl.draw_graph(Dev_1, color = ['dpt_pseudotime'])
sc.pl.paga_compare(Dev_1)


# In[8]:


sc.pl.violin(Dev_1, ['dpt_pseudotime'], groupby='identity_2', rotation = 90)


# In[9]:


###save data for visualization
Dev_1_FDG_df = pd.DataFrame(
    Dev_1.obsm['X_draw_graph_fa'], 
    index = Dev_1.obs.index, columns = ['x', 'y']).join(Dev_1.obs[['identity_2', 'dpt_pseudotime']])
Dev_1_FDG_df.to_csv('Dev_1_FDG_df.csv')

Dev_1_paga_pos_df = pd.DataFrame(
    Dev_1.uns['paga']['pos'], index = Dev_1.obs['identity_2'].cat.categories, columns = ['x', 'y'])
Dev_1_paga_pos_df.to_csv('Dev_1_paga_pos_df.csv')

Dev_1_paga_con_df = pd.DataFrame(
    Dev_1.uns['paga']['connectivities'].toarray(), 
    index = Dev_1.obs['identity_2'].cat.categories, 
    columns = Dev_1.obs['identity_2'].cat.categories)
Dev_1_paga_con_df.to_csv('Dev_1_paga_con_df.csv')

Dev_1.write_h5ad('Dev_1.h5ad')


# In[ ]:





# In[ ]:





# In[10]:


Dev_2 = sc.read_loom('Dev_2.loom')
Dev_2_cp = Dev_2
Dev_2_metadata = pd.read_csv('Dev_2_metadata.csv', index_col= 0)
Dev_2_harmony = np.array(pd.read_csv('Dev_2_harmony.csv', index_col= 0))
Dev_2.obs = Dev_2_metadata
Dev_2.obsm['X_pca'] = Dev_2_harmony
sc.pp.neighbors(Dev_2, n_pcs=20)
sc.tl.diffmap(Dev_2)
Dev_2.uns['iroot'] = np.flatnonzero(Dev_2.obs['identity_2']  == 'DP_Q')[0]
sc.tl.dpt(Dev_2)
sc.tl.paga(Dev_2, groups='identity_2')


# In[11]:


sc.pl.paga(Dev_2, color = ['identity_2'])
sc.pl.paga(Dev_2, color=['dpt_pseudotime'])

sc.pl.diffmap(Dev_2, color = 'identity_2')
sc.pl.diffmap(Dev_2, color = 'dpt_pseudotime')


# In[12]:


sc.tl.draw_graph(Dev_2, maxiter = 100, layout = 'fa',init_pos='paga')
sc.pl.draw_graph(Dev_2, color = ['identity_2'])
sc.pl.draw_graph(Dev_2, color = ['dpt_pseudotime'])
sc.pl.paga_compare(Dev_2)


# In[14]:


sc.pl.violin(Dev_2, ['dpt_pseudotime'], groupby='identity_2', rotation = 90)


# In[15]:


###save data for visualization
Dev_2_FDG_df = pd.DataFrame(
    Dev_2.obsm['X_draw_graph_fa'], 
    index = Dev_2.obs.index, columns = ['x', 'y']).join(Dev_2.obs[['identity_2', 'dpt_pseudotime']])
Dev_2_FDG_df.to_csv('Dev_2_FDG_df.csv')

Dev_2_paga_pos_df = pd.DataFrame(
    Dev_2.uns['paga']['pos'], index = Dev_2.obs['identity_2'].cat.categories, columns = ['x', 'y'])
Dev_2_paga_pos_df.to_csv('Dev_2_paga_pos_df.csv')

Dev_2_paga_con_df = pd.DataFrame(
    Dev_2.uns['paga']['connectivities'].toarray(), 
    index = Dev_2.obs['identity_2'].cat.categories, 
    columns = Dev_2.obs['identity_2'].cat.categories)
Dev_2_paga_con_df.to_csv('Dev_2_paga_con_df.csv')

Dev_2.write_h5ad('Dev_2.h5ad')

