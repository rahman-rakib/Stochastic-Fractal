import numpy as np
from scipy.special import gamma

n_star = np.linspace(0,1,100_001)

def func_trans(n_star,alpha,p):
    """left-hand side of transcendental equation defined as a numpy array"""
    return (gamma(n_star+alpha)*gamma(2*alpha))/(gamma(n_star+2*alpha)*gamma(alpha))-1/(1+p)

def find_df(alpha,p):
    """root finding by minimizing the cost function: square of left-hanad side"""
    cost_array = func_trans(n_star,alpha,p)**2
    min_index = np.argmin(cost_array)
    df_value = n_star[min_index]
    df_rounded = float("{:.6f}".format(df_value))
    return df_rounded

def test():
    print(find_df(2, 0.75))
    print(find_df(3, 0.75))