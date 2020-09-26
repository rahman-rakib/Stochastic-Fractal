import numpy as np
import scipy.special as sc
from main_py import df_determination

find_df = df_determination.find_df
analytic_data_dir = "../data/alpha3/analytic/"  # data directory
temp_fig_dir = "../figures/temp/"


def phi_value_alpha_3(xi, p):
    xi_list = np.linspace(0, 2, 10_001)
    if abs(p - 0.1) <= 1e-6:
        filename = 'phi_0_10_list.csv'
    elif abs(p - 0.5) <= 1e-6:
        filename = 'phi_0_50_list.csv'
    elif abs(p - 0.9) <= 1e-6:
        filename = 'phi_0_90_list.csv'
    elif abs(p - 0.75) <= 1e-6:
        filename = 'phi_0_75_list.csv'
    else:
        print('analytical values available only for p equal 0.10, 0.50, 0.75 and 0.90')
        return None

    phi_list = np.loadtxt(analytic_data_dir + filename)

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


class AnalyticSoln:
    def __init__(self, alpha, probability):
        self.xi_list = None
        self.phi_list = None
        filename = "phi_list_analytic_alpha_3.txt"
        analytic_data_dir = "../data/alpha3/analytic/"  # data directory
        self.filename = analytic_data_dir + filename
        self.data_dct = None
        self.defined_p = [0.1, 0.5, 0.75, 0.9]

        # read and load data
        self.read_and_load(probability)
        pass

    def read_and_load(self, probability):
        ret = self.check_if_in_list(probability)
        if ret is None:
            print('analytical values available only for p equal ', self.defined_p)
            return None

        p = ret  # so that the values is exactly what it is supposed to be
        if self.data_dct is None:
            # load data only once
            phi_list_tmp = np.loadtxt(self.filename)
            self.xi_list = phi_list_tmp[:, 0]
            self.data_dct = dict()
            self.data_dct[0.1] = phi_list_tmp[:, 1]
            self.data_dct[0.5] = phi_list_tmp[:, 2]
            self.data_dct[0.75] = phi_list_tmp[:, 3]
            self.data_dct[0.9] = phi_list_tmp[:, 4]
            pass

        self.phi_list = self.data_dct[p]
        pass

    def check_if_in_list(self, p):
        flag = False
        for pp in self.defined_p:
            flag = abs(p - pp) <= 1e-6
            if flag:
                return pp
            pass
        return None

    def phi_value_alpha_3_v2(self, xi):

        if xi > 2:
            phi_value = 0
        elif xi in self.xi_list:
            index = list(self.xi_list).index(xi)
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
    pass






def phi_list_alpha_3(xi_list, p):
    phi_list = []
    for xi in xi_list:
        phi_value = phi_value_alpha_3(xi, p)
        phi_list.append(phi_value)

    return np.array(phi_list)


def phi_list(alpha, p, xi_list):
    df = find_df(alpha, p)

    if alpha == 1:
        density = np.exp(-xi_list)

    elif alpha == 2:
        c2 = -(sc.gamma(1 / 3) / sc.gamma(5 / 3)) * (sc.gamma((df + 5) / 3) / sc.gamma((df + 3) / 3))
        density1 = -3 * (df + 2) * (xi_list ** 2) * sc.hyp1f1(-(df - 1) / 3, 4 / 3, -xi_list ** 3)
        density2 = -(3 / 5) * c2 * df * (xi_list ** 4) * sc.hyp1f1(-(df - 3) / 3, 8 / 3, -xi_list ** 3)
        density3 = -2 * c2 * (xi_list) * sc.hyp1f1(-df / 3, 5 / 3, -xi_list ** 3)
        density = density1 + density2 + density3

    elif alpha == 3:
        density = phi_list_alpha_3(xi_list, p)
        # density = phi_value_alpha_3_v2(xi_list, p)
    else:
        density = None
        print('the analytical solution for alpha = {} value is unknown'.format(alpha))

    return density


def test():
    import matplotlib.pyplot as plt
    xi_list = np.linspace(0, 2, 5000)
    for alphap in range(2, 4):
        print("doing alpha = ", alphap)
        # if alphap == 1:
        #     xi_list = np.linspace(-3, 10, 5000)
        phi = phi_list(alphap, 0.75, xi_list)

        plt.plot(xi_list, phi, label="alpha={}".format(alphap))
        pass
    alphap = 1
    xi_list = np.linspace(-3, 10, 5000)
    phi = phi_list(alphap, 0.75, xi_list)
    plt.plot(xi_list, phi, label="alpha={}".format(alphap))
    plt.legend()
    plt.savefig(temp_fig_dir + "testing-analytic-plot2")
