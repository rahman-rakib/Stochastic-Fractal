# from main_py import df_determination
import numpy as np
import scipy.special as sc
import pandas as pd
import decimal as dc
import os
from pathlib import Path
from scipy.special import gamma

dc.getcontext().rounding = dc.ROUND_DOWN


# find_df = df_determination.find_df
# analytic_data_dir = "../data/alpha3/analytic/"  # data directory
temp_fig_dir = "../figures/temp/"
# dir_data = "../../data/alpha3/"
# data_filename = dir_data+"function_phi_alpha3_dataframe.csv"


class AnalyticSoln:
    """
    Analytic solution of the Binary Fragmentation Equation.
    for alpha=1,2,3 only
    """
    def __init__(self):
        filename = "function_phi_alpha3_dataframe.csv"
        root_dir = Path(os.path.dirname(os.path.abspath(__file__))).parent
        analytic_data_dir = str(root_dir) + "/data/alpha3/"  # data directory
        self.filename = os.path.abspath(analytic_data_dir + filename)
        # print("analytic function data file ", self.filename)
        # self.data_dct = None
        self.df_data = None
        self.n_star = np.linspace(0, 1, 1_000_001)
        # read and load data
        # self.read_and_load(probability)
        pass

    def func_trans(self, alpha, probability_p, n_star=None):
        """
        left-hand side of transcendental equation defined as a numpy array
        :param alpha:
        :param probability_p:
        :param n_star:
        :return:
        """
        if n_star is not None:
            self.n_star = n_star
        return (gamma(self.n_star + alpha) * gamma(2 * alpha)) / (gamma(self.n_star + 2 * alpha) * gamma(alpha)) - 1 / (1 + probability_p)

    def find_df(self, alpha, p, n_star=None):
        """root finding by minimizing the cost function: square of left-hanad side"""
        cost_array = self.func_trans(alpha, p, n_star) ** 2
        min_index = np.argmin(cost_array)
        df_value = self.n_star[min_index]
        df_rounded = float("{:.6f}".format(df_value))
        return df_rounded

    def phi_value_alpha_3(self, xi, probability_p):
        """
        df_data : pandas DataFrame
        xi      :
        probability_p       : probablity
        """

        p_0 = float(round(dc.Decimal(probability_p), 2))
        p_1 = float(round(dc.Decimal(p_0 + 0.01), 2))
        del_p = abs(probability_p - p_0)

        if del_p <= 1e-9:
            col_name = 'p={:.2f}'.format(p_0)
            phi_list = self.df_data[col_name].values

        else:
            col_name0 = 'p={:.2f}'.format(p_0)
            phi_list0 = self.df_data[col_name0].values

            col_name1 = 'p={:.2f}'.format(p_1)
            phi_list1 = self.df_data[col_name1].values

            del_phi_list = (phi_list1 - phi_list0) * del_p / 0.01
            phi_list = phi_list0 + del_phi_list

        xi_list = self.df_data['xi'].values
        if xi > 2:
            phi_value = 0
        elif xi in xi_list:
            index = list(xi_list).index(xi)
            phi_value = phi_list[index]
        else:
            index0 = int(xi // 0.0002)
            index1 = int(index0 + 1)
            del_xi = xi % 0.0002
            phi_0 = phi_list[index0]
            phi_1 = phi_list[index1]
            del_phi = (phi_1 - phi_0) * del_xi / 0.0002
            phi_value = phi_0 + del_phi

        return phi_value

    def phi_list_alpha_3(self, xi_list, p):
        phi_list = []
        for xi in xi_list:
            phi_value = self.phi_value_alpha_3(xi, p)
            phi_list.append(phi_value)

        return np.array(phi_list)

    def phi_list(self, alpha, probability_p, xi_list=None):
        """

        :param alpha:           value of alpha
        :param probability_p:   value of probability p
        :param xi_list:         independent variable xi. If None then a predefined xi_list is used and retured
        :return:
        """
        if xi_list is None:
            xi_list = np.linspace(0, 2, 10_001)
            pass

        df = self.find_df(alpha, probability_p)

        if alpha == 1:
            density = np.exp(-xi_list)

        elif alpha == 2:
            c2 = -(sc.gamma(1 / 3) / sc.gamma(5 / 3)) * (sc.gamma((df + 5) / 3) / sc.gamma((df + 3) / 3))
            density1 = -3 * (df + 2) * (xi_list ** 2) * sc.hyp1f1(-(df - 1) / 3, 4 / 3, -xi_list ** 3)
            density2 = -(3 / 5) * c2 * df * (xi_list ** 4) * sc.hyp1f1(-(df - 3) / 3, 8 / 3, -xi_list ** 3)
            density3 = -2 * c2 * (xi_list) * sc.hyp1f1(-df / 3, 5 / 3, -xi_list ** 3)
            density = density1 + density2 + density3

        elif alpha == 3:
            if self.df_data is None:
                # Read data only once.
                self.df_data = pd.read_csv(self.filename)
            density = self.phi_list_alpha_3(xi_list, probability_p)

        else:
            print('the analytical solution for alpha = {} value is unknown'.format(alpha))
            return

        return xi_list, density
    pass

def test():
    alpha_array = np.linspace(2, 3, 1_000)
    soln = AnalyticSoln()
    p = 0.5
    for i in alpha_array:
        df_value = soln.find_df(i, p)
        print("alpha = ", i, " df ", df_value)

def test_2():
    import matplotlib.pyplot as plt
    xi_list = np.linspace(0, 2, 5000)
    sonl = AnalyticSoln()
    for alphap in range(2, 4):
        print("doing alpha = ", alphap)
        # if alphap == 1:
        #     xi_list = np.linspace(-3, 10, 5000)
        xi_list, phi = sonl.phi_list(alphap, 0.75, xi_list)
        plt.plot(xi_list, phi, label="alpha={}".format(alphap))
        pass
    # alphap = 1
    # xi_list = np.linspace(-3, 10, 5000)
    # phi = phi_list(alphap, 0.75, xi_list)
    # plt.plot(xi_list, phi, label="alpha={}".format(alphap))
    plt.legend()
    plt.savefig(temp_fig_dir + "testing-analytic-plot2")
