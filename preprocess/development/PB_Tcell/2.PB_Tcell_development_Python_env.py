#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
import scanpy as sc
import pandas as pd
import loompy as lp


# In[2]:


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300, facecolor='white')


# In[3]:


PB_CD4T = sc.read_loom('PB_CD4T.loom')
PB_CD4T_cp = PB_CD4T
PB_CD4T_metadata = pd.read_csv('PB_CD4T_metadata.csv', index_col= 0)
PB_CD4T_harmony = np.array(pd.read_csv('PB_CD4T_harmony.csv', index_col= 0))
PB_CD4T.obs = PB_CD4T_metadata
PB_CD4T.obsm['X_pca'] = PB_CD4T_harmony

sc.pp.neighbors(PB_CD4T, n_pcs=30)
sc.tl.diffmap(PB_CD4T)
PB_CD4T.uns['iroot'] = np.flatnonzero(PB_CD4T.obs['identity']  == 'CD4+RTE')[0]
sc.tl.dpt(PB_CD4T)
sc.tl.paga(PB_CD4T, groups='identity')


# In[4]:


sc.pl.paga(PB_CD4T, color = ['identity'])
sc.pl.paga(PB_CD4T, color=['dpt_pseudotime'])
sc.pl.diffmap(PB_CD4T, color = 'identity')
sc.pl.diffmap(PB_CD4T, color = 'dpt_pseudotime')


# In[5]:


sc.tl.draw_graph(PB_CD4T, maxiter = 100, layout = 'fa',init_pos='paga')
sc.pl.draw_graph(PB_CD4T, color = ['identity'])
sc.pl.draw_graph(PB_CD4T, color = ['dpt_pseudotime'])
sc.pl.paga_compare(PB_CD4T)


# In[7]:


sc.pl.violin(PB_CD4T, ['dpt_pseudotime'], groupby='identity', rotation = 90)


# In[8]:


###save data for visualization
PB_CD4T.write_h5ad('PB_CD4T.h5ad')

PB_CD4T_FDG_df = pd.DataFrame(
    PB_CD4T.obsm['X_draw_graph_fa'], 
    index = PB_CD4T.obs.index, columns = ['x', 'y']).join(PB_CD4T.obs[['identity', 'dpt_pseudotime']])
PB_CD4T_FDG_df.to_csv('PB_CD4T_FDG_df.csv')

PB_CD4T_paga_pos_df = pd.DataFrame(
    PB_CD4T.uns['paga']['pos'], index = PB_CD4T.obs['identity'].cat.categories, columns = ['x', 'y'])
PB_CD4T_paga_pos_df.to_csv('PB_CD4T_paga_pos_df.csv')

PB_CD4T_paga_con_df = pd.DataFrame(
    PB_CD4T.uns['paga']['connectivities'].toarray(), 
    index = PB_CD4T.obs['identity'].cat.categories, 
    columns = PB_CD4T.obs['identity'].cat.categories)
PB_CD4T_paga_con_df.to_csv('PB_CD4T_paga_con_df.csv')


# In[ ]:





# In[ ]:





# In[22]:


PB_CD8T = sc.read_loom('PB_CD8T.loom')
PB_CD8T_cp = PB_CD8T
PB_CD8T_metadata = pd.read_csv('PB_CD8T_metadata.csv', index_col= 0)
PB_CD8T_harmony = np.array(pd.read_csv('PB_CD8T_harmony.csv', index_col= 0))
PB_CD8T.obs = PB_CD8T_metadata
PB_CD8T.obsm['X_pca'] = PB_CD8T_harmony
sc.pp.neighbors(PB_CD8T, n_pcs=15)
sc.tl.diffmap(PB_CD8T)
PB_CD8T.uns['iroot'] = np.flatnonzero(PB_CD8T.obs['identity']  == 'CD8+RTE')[0]
sc.tl.dpt(PB_CD8T)
sc.tl.paga(PB_CD8T, groups='identity')


# In[23]:


sc.pl.paga(PB_CD8T, color = ['identity'])
sc.pl.paga(PB_CD8T, color=['dpt_pseudotime'])
sc.pl.diffmap(PB_CD8T, color = 'identity')
sc.pl.diffmap(PB_CD8T, color = 'dpt_pseudotime')


# In[24]:


sc.tl.draw_graph(PB_CD8T, maxiter = 100, layout = 'fa',init_pos='paga')
sc.pl.draw_graph(PB_CD8T, color = ['identity'])
sc.pl.draw_graph(PB_CD8T, color = ['dpt_pseudotime'])
sc.pl.paga_compare(PB_CD8T)


# In[25]:


sc.pl.paga(PB_CD8T, color = ['identity'], pos = PB_CD8T.uns['paga']['pos'])
sc.pl.paga(PB_CD8T, color=['dpt_pseudotime'], pos = PB_CD8T.uns['paga']['pos'])


# In[27]:


sc.pl.violin(PB_CD8T, ['dpt_pseudotime'], groupby='identity', rotation = 90)


# In[ ]:





# In[ ]:





# In[28]:


PB_CD8T_FDG_df = pd.DataFrame(
    PB_CD8T.obsm['X_draw_graph_fa'], 
    index = PB_CD8T.obs.index, columns = ['x', 'y']).join(PB_CD8T.obs[['identity', 'dpt_pseudotime']])
PB_CD8T_FDG_df.to_csv('PB_CD8T_FDG_df.csv')

PB_CD8T_paga_pos_df = pd.DataFrame(
    PB_CD8T.uns['paga']['pos'], index = PB_CD8T.obs['identity'].cat.categories, columns = ['x', 'y'])
PB_CD8T_paga_pos_df.to_csv('PB_CD8T_paga_pos_df.csv')

PB_CD8T_paga_con_df = pd.DataFrame(
    PB_CD8T.uns['paga']['connectivities'].toarray(), 
    index = PB_CD8T.obs['identity'].cat.categories, 
    columns = PB_CD8T.obs['identity'].cat.categories)
PB_CD8T_paga_con_df.to_csv('PB_CD8T_paga_con_df.csv')


# In[ ]:





# In[ ]:





# In[ ]:


PB_CD8T.write_h5ad('PB_CD8T.h5ad')


# In[ ]:




