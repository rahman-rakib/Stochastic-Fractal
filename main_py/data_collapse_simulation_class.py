#!/usr/bin/env python
# coding: utf-8

# In[1]:


import random, numpy as np, pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import os


# ## Loading class module

# In[2]:


get_ipython().run_line_magic('run', '../stochastic_fractal_m.py')


# In[3]:


alphap = 1  # any real positive value larger than 0.5, but anaylytical plot exists only for alpha=1,2,3
probp = 0.75 # any value between 0 and 1

ensemble_sizep = 20_000
min_iterationp = 100_000
bin_sizep = 0.0001


# In[4]:


expon = 2 * alphap - 1


# ## Defining directory

# In[5]:


dir_data = "../data/alpha{}/".format(alphap)
dir_fig  ="../figures/alpha{}/".format(alphap)


## While testing
dir_data = "../data/temp/alpha{}/".format(alphap)
dir_fig  ="../figures/temp/alpha{}/".format(alphap)

# create directory if it does not exists
for dir_name in [dir_data, dir_fig]:
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
        print("directory ", dir_name, " is created")
    else:
        print("directory ", dir_name, " exists")


# In[6]:


figure_name_collapse = dir_fig + "pn_data_coll_alpha_{}.pdf".format(alphap)
figure_name = dir_fig + "pn_alpha_{}.pdf".format(alphap)


# In[ ]:





# ## Simulation using class

# In[7]:


stochastic_frag = TrueLengths(alpha=alphap, probability=probp)


# In[8]:


# total_iteration = 10000
# lengths_enselbme = stochastic_frag.run_ensemble(ensemble_sizep, total_iteration)


# In[9]:


ana_soln = AnalyticSoln()


# ## Other functions

# In[10]:


def bin_maker(max_value, bin_number):
    bin_size = max_value/bin_number
    bins = [0.]
    for i in range(bin_number):
        bin_edge = (i+1)*bin_size
        bins.append(bin_edge)
        pass
    return bins


# In[11]:


def histrogram_data(data, bin_size = bin_sizep):
    bin_number = int((np.max(data) - np.min(data))//bin_size) + 1
    y, x = np.histogram(data, bins = bin_number)
    return x[1:], y


# ## fitting data
# 
# $\phi \equiv \text{simulation}$
# $\phi^{(th)} \equiv \text{theoritical or analytical}$
# 
# $$cost = \sum_t \sum_i (b*\phi_{i t} - \phi_{i t}^{(th)})^2 \$$
# 
# we need to minimize $cost$. We get the value of $b$ as
# 
# $$b = \frac{\sum_t \sum_i \phi_{i t} \phi_{i t}^{(th)}}{\sum_t \sum_i \phi_{i t}^2}$$

# In[12]:


def fitting_parameter(given_array,ref_array):
    a1 = np.multiply(given_array,ref_array)
    a2 = np.multiply(given_array,given_array)
    a1_sum = np.sum(a1)
    a2_sum = np.sum(a2)
    return a1_sum, a2_sum


# ## Simulation  and plot

# In[13]:


def plot_data_simulation(total_iteration, ens_data):
    """
    data collapse points from simulation
    """
    df = ana_soln.find_df(alphap,probp)
#     ens_data = stochastic_frag.run_ensemble(ensemble_sizep, total_iteration)  # from class
    x,y = histrogram_data(ens_data)
    xi = x*total_iteration**(1/(2*alphap-1))
    phi = y/(total_iteration**((1+df)/(2*alphap-1)))

    return xi, phi 


# In[14]:


def plot_data_fitted(total_iteration, ens_data):
    xi, phi  = plot_data_simulation(total_iteration, ens_data)
    xi, phi_analytic = ana_soln.phi_list(alphap,probp,xi)
    a1_sum, a2_sum = fitting_parameter(phi,phi_analytic)
    return xi, phi, a1_sum, a2_sum


# In[15]:


# def plot_data(total_iteration):
#     df = ana_soln.find_df(alphap,probp)
    
#     # doing the simulation here
#     ens_data = stochastic_frag.run_ensemble(ensemble_sizep, total_iteration)  # from class
#     x,y = histrogram_data(ens_data)
#     xi = x*total_iteration**(1/(2*alphap-1))
#     phi = y/(total_iteration**((1+df)/(2*alphap-1)))
#     xi, phi_analytic = ana_soln.phi_list(alphap,probp,xi)  # from class
# #     print(len(phi_analytic))
#     a1_sum, a2_sum = fitting_parameter(phi, phi_analytic)
#     print(a1_sum, ", ", a2_sum)
#     return xi, phi, a1_sum, a2_sum


# ## Simulation

# In[16]:


ens_data_dct = dict()
# for i in range(3):
#     time_t = min_iterationp*(i+1)
#     ens_data_dct[time_t] = stochastic_frag.run_ensemble(ensemble_sizep, time_t)  # from class
#     print(time_t)


# In[18]:


time_t = min_iterationp*1
ens_data_dct[time_t] = stochastic_frag.run_ensemble_parallel(ensemble_sizep, time_t, 10)  # from class
# print(time_t)
filename = "alpha_1_p_75_true_lengths_t_{}".format(time_t)
np.savetxt(filename, ens_data_dct[time_t])


# In[19]:


time_t = min_iterationp*2
ens_data_dct[time_t] = stochastic_frag.run_ensemble_parallel(ensemble_sizep, time_t, 10)  # from class
# print(time_t)
filename = "alpha_1_p_75_true_lengths_t_{}".format(time_t)
np.savetxt(filename, ens_data_dct[time_t])


# In[ ]:


time_t = min_iterationp*3
ens_data_dct[time_t] = stochastic_frag.run_ensemble_parallel(ensemble_sizep, time_t, 10)  # from class
# print(time_t)
filename = "alpha_1_p_75_true_lengths_t_{}".format(time_t)
np.savetxt(filename, ens_data_dct[time_t])


# In[ ]:


# ens_data_dct[time_t]


# In[ ]:


data_dict = {}
a1_sum, a2_sum = 0,0
for time_t in ens_data_dct.keys():
   
    if alphap in [1,2,3]:
        print("alpha is within known values")
        xi, phi, a1, a2 = plot_data_fitted(time_t, ens_data_dct[time_t])
    else:
        xi, phi = plot_data_simulation(time_t,  ens_data_dct[time_t])
        a1, a2 = 1, 1
        pass
    data_dict[time_t]=[xi,phi]
    a1_sum += a1
    a2_sum += a2
    
for i in range(3):
    time_t = min_iterationp*(i+1)
    data_dict[time_t][1]*=(a1_sum/a2_sum)


# In[ ]:


# print(data_dict.keys())
# print(a1_sum)
# print(a2_sum)


# ### Saving data for later use

# In[ ]:


signature = 'data_collapse_alpha_{}_t_{}k'
for time_t in data_dict.keys():
    x,y = data_dict[time_t]
    file_name = signature.format(alphap, int(time_t//1000))
#     np.savetxt(dir_data+file_name,np.c_[x,y])


# ## Plotting

# In[ ]:


fig, axes = plt.subplots(1,1,figsize = (5,3.5),dpi = 300)

for i in range(3):
    time_t = min_iterationp*(i+1)
    x,y = data_dict[time_t]
    plt.plot(x, y,"o", markersize=1,label=r"$t={}k$".format(str(int(time_t//1000))))

if alphap in [1,2,3]:
    xi_th = np.linspace(0,2,1000)
    xis, phi_th = ana_soln.phi_list(alphap,probp,xi_th)  # from class
    plt.plot(xi_th, phi_th,color='black',linewidth=0.5)

plt.legend(loc=1, prop={'size': 12})
plt.xlabel(r"$x/t^{1/z}$",  fontsize=16)
plt.ylabel(r"$c(x,t)/t^\theta$",  fontsize=16)
plt.xlim([0, 2.1])
plt.text(0.8, 0.55, r"$\beta$={}".format(alphap-1), transform = axes.transAxes)
plt.text(0.8, 0.5, r"$p$={}".format(probp), transform = axes.transAxes)

axes.set_position([0.15, 0.15, 0.8, 0.8])

print(figure_name_collapse)
plt.savefig(figure_name_collapse)


# In[ ]:


fig, axes = plt.subplots(1,1,figsize = (5,3.5),dpi = 300)

for i in range(3):
    time_t = min_iterationp*(i+1)
    x,y = histrogram_data(ens_data_dct[time_t])
    plt.plot(x, y,"o", markersize=1,label=r"$t={}k$".format(str(int(time_t//1000))))


plt.legend(loc=1, prop={'size': 12})
plt.xlabel(r"$x$", fontsize=16)
plt.ylabel(r"$c(x,t)$", fontsize=16)
plt.text(0.8, 0.55, r"$\beta$={}".format(alphap-1), transform = axes.transAxes)
plt.text(0.8, 0.5, r"$p$={}".format(probp), transform = axes.transAxes)
# plt.xlim([0, 2.1])
axes.set_position([0.15, 0.15, 0.8, 0.8])


print(figure_name)
plt.savefig(figure_name)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


import numpy as np
from time import time

# Prepare data
np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[200000, 5])
data = arr.tolist()
data[:5]


# In[ ]:


# Solution Without Paralleization

def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

results = []
for row in data:
    results.append(howmany_within_range(row, minimum=4, maximum=8))

print(results[:10])
#> [3, 1, 4, 4, 4, 2, 1, 1, 3, 3]


# In[ ]:


# Parallelizing using Pool.apply()

import multiprocessing as mp

# Step 1: Init multiprocessing.Pool()
pool = mp.Pool(mp.cpu_count())

# Step 2: `pool.apply` the `howmany_within_range()`
results = [pool.apply(howmany_within_range, args=(row, 4, 8)) for row in data]

# Step 3: Don't forget to close
pool.close()    

print(results[:10])
#> [3, 1, 4, 4, 4, 2, 1, 1, 3, 3]


# In[ ]:




