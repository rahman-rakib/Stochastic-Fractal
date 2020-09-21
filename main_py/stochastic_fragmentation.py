import random
import numpy as np
import datetime
import json
import time


class StochasticFragmentation:
    """
    Simulatin for Stochastic Binary Fragmentation
    """
    def __init__(self, **kwargs):
        """

        :param kwargs: following keywords are defined
            alpha :
            beta  :
            probability :
        """
        print("kwargs ", kwargs)
        self.alpha = kwargs['alpha']
        self.prob = kwargs['probability']

        # variables
        self.expon = 2 * self.alpha - 1
        self.normC = 1.0  # normalization constant
        self.probability_list = [1.0]
        self.length_list = [1.0]
        self.flag_list = [True]
        self.choose_pivot = self.betadist

        # how many times
        self.ensemble_size = 1 # at least once
        self.time_iteration = 0 # it will be set when run() executes

        # loggin
        self.logging = False
        pass

    def log(self, flag=False):
        self.logging = flag

    def get_signature(self):
        return "StochasticFragmentation"

    def reset(self):
        self.normC = 1.0  # normalization constant
        self.probability_list = [1.0]
        self.length_list = [1.0]
        self.flag_list = [True]

    def set_pivot_choosing_method(self, method='logistic'):
        """
        method  : way to choose the pivot point of a given segment
        returns : random value in [0, 1]
        """
        if method is 'logistic':
            self.choose_pivot = self.logistic_choice
            pass
        self.choose_pivot = self.betadist
        pass

    def betadist(self):
        """gives a random number from beta distribution"""
        #         print("betadist")
        return random.betavariate(self.alpha, self.alpha)

    def logistic_xn(self, n, a=2):
        x0 = random.random()
        x1 = 0
        if n < 0:
            print("n cannot be negative")
            n = 1
            pass
        for i in range(n):
            x1 = a * x0 * (1 - x0)
            x0 = x1
            pass
        return x0

    def logistic_choice(self):
        #         print("logistic_choice")
        beta = self.alpha - 1
        RR = self.logistic_xn(beta)
        RRp = 1 - RR
        r = random.random()
        if r < 0.5:
            return RR
        return RRp

    def decision(self):
        """
        decides with a given probability whether to keep the right part
        """
        if self.prob > random.random():
            return True
        else:
            return False

    def splitting(self, segment):
        """
        splits a given segment. left and right are endpoints of the segment
        segment : length of a segment
        returns :
            xL -> length of the left segment
            xR -> length of the right segment
            flag -> keeping the right segment
            xLp, xRp -> probability(unnormalized) for being selected
            change -> change of normalization const
        """
        xL = segment * self.choose_pivot()
        xR = segment - xL
        flag = self.decision()
        xLp = xL ** self.expon
        xRp = xR ** self.expon
        change = xLp + xRp - segment ** self.expon
        return xL, xR, flag, xLp, xRp, change

    def pickindex(self):
        """
        picks up a segment to be subsequently split
        """
        # return self.pickindex_v1()
        return self.pickindex_v2()


    def pickindex_v1(self):
        """
        from old algorighm
        picks up a segment to be subsequently split
        """
        r = random.uniform(0, 1)
        sum_ = 0
        for index in range(len(self.probability_list)):
            sum_ += self.probability_list[index] / self.normC
            if sum_ < r:
                continue
            else:
                return index
            pass
        print("out of range. return None")
        pass

    def pickindex_v2(self):
        """
        NEW algorithm. Seems to be accurate one. agrees with Theoritical funciton.
        picks up a segment to be subsequently split
        """
        # print("pickindex_v2")
        # print("self.probability_list ", self.probability_list)
        r = random.uniform(0, 1.0)  # 1 is the initial particle size
        if r > self.normC:
            # print("r > self.normC => return None")
            return None
        sum_ = 0
        for index in range(len(self.probability_list)):
            sum_ += self.probability_list[index]
            # print("[", index, "] => ", self.probability_list[index], " .cumsum = ", sum_)
            if sum_ < r:
                continue
            else:
                return index
            pass
        print("out of range. return None")
        return None

    def view(self):
        print("viewing status")
        print("alpha = ", self.alpha)
        print("keeping probability ", self.prob)
        print("<length>  <probability>  <flag>")
        for i in range(len(self.length_list)):
            print("{:.5e}, {:.5e}, {:5}".format(self.length_list[i],
                                                self.probability_list[i] / self.normC,
                                                self.flag_list[i]
                                                ))
            pass
        print(np.sum(self.length_list))
        print(np.sum(self.length_list))
        pass

    def one_time_step(self):
        # print("StochasticFragmentation.one_time_step")
        index = self.pickindex()
        # print("index ", index)
        if (type(index) == int) and self.flag_list[index]:
            xL, xR, flag, xLp, xRp, change = self.splitting(self.length_list[index])

            self.length_list[index] = xL
            self.length_list.append(xR)
            self.flag_list.append(flag)
            self.probability_list[index] = xLp
            self.probability_list.append(xRp)
            self.normC += change
            pass

        pass

    # def run(self, time_iteration):
    #
    #     for i in range(time_iteration + 1):
    #         #             print("time step ", i)
    #         self.one_time_step()
    #
    #         lengths = np.array(self.length_list)
    #         lengths = lengths[self.flag_list]
    #
    #
    #     pass

    def get_header(self, header_description="NA"):
        x = datetime.datetime.now()
        date_time = x.strftime("%Y%m%d_%H%M%S")
        j_header = dict()
        j_header['time_iteration'] = self.time_iteration
        j_header['date_time'] = date_time
        j_header['ensemble_size'] = self.ensemble_size
        j_header['alpha'] = self.alpha
        j_header['probability'] = self.prob
        j_header['desc'] = header_description
        self.header_str = json.dumps(j_header)
        return j_header

    def get_filename(self):
        signature = self.get_signature()
        filename = signature
        filename += "_time-iteration_{}".format(self.time_iteration)
        filename += "_ensemble_{}_".format(self.ensemble_size)
        header = self.get_header()
        date_time = header['date_time']
        filename += date_time
        filename += ".txt"
        return filename


# df
class NumberLength(StochasticFragmentation):
    """

    """
    def get_signature(self):
        return super(NumberLength, self).get_signature() + "NumberLength"

    def number_length(self):
        lengths = np.array(self.length_list)
        lengths = lengths[self.flag_list]
        segment_count = lengths.shape[0]
        surviving_length_sum = np.sum(lengths)

        #         surviving_length_sum, N = 0, 0
        #         for i in range(len(self.flag_list)):
        #             if self.flag_list[i]:
        #                 N += 1
        #                 surviving_length_sum += self.length_list[i]
        #                 pass
        #             pass
        #         if abs(M1-surviving_length_sum) > 1e-10:
        #             print("not equal")
        return segment_count, surviving_length_sum

    def run(self, time_iteration, min_iteration, number_of_points):
        """
        we run the `one_time_step()` method `time_iteration` times. We record some information
        at `iteration_step` step interval starting from `min_iteration`
        :param time_iteration: maximum time step
        :param min_iteration: starting point to record data
        :param number_of_points: number of data points
        :return: a numpy array with two columns,
                1st column is the particle counts and
                2nd column is the sum of surviving particle sizes
        """
        # iteration_list = list(range(min_iteration, time_iteration + 1, iteration_step))
        N_realization = []
        M_realization = []
        step_size = int((time_iteration - min_iteration)/number_of_points)
        for i in range(time_iteration + 1):
            #             print("time step ", i)
            self.one_time_step()
            if (i > min_iteration) and (i % step_size == 0):

            # if i + 1 in iteration_list:
                segment_count, surviving_length_sum = self.number_length()
                N_realization.append(segment_count)
                M_realization.append(surviving_length_sum)
            pass

        N_list = np.array(N_realization)
        M_list = np.array(M_realization)
        self.time_iteration = time_iteration
        return np.c_[N_list, M_list]

    def run_ensemble(self, ensemble_size, time_iteration, start_at, number_of_data_points):
        """

        :param ensemble_size: ensemble size
        :param time_iteration: total number of time steps
        :param start_at:       minimum number of step before we start recording data
        :param number_of_data_points: number of data poitns. starting from `start_at` it will take `number_of_data_points`
                data points upto time `time_iteration`
        :return: [N, M] average value
                where, N = number of particles at particular time steps
                       M = sum of surviving length at those time steps
        """

        ensemble_data = None
        step_temp=int(ensemble_size / 100)
        start_time = time.time()
        for i in range(ensemble_size):
            if i % step_temp == 0 and self.logging:
                print("realization ", i, " . Time spent ", (time.time()-start_time), " sec")
                start_time = time.time()
                pass
            self.reset()
            out_data_N_M = self.run(time_iteration, start_at, number_of_data_points)
            if ensemble_data is None:
                ensemble_data = out_data_N_M
            else:
                ensemble_data += out_data_N_M
            pass

        data_average = ensemble_data / ensemble_size

        self.time_iteration = time_iteration
        self.ensemble_size = ensemble_size

        return data_average
    pass


# xdf
class Moment(StochasticFragmentation):
    """

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        key = "fractal_dim"
        self.exponent = None
        if key in kwargs.keys():
            self.exponent = kwargs[key]
            pass
        key = "exponent"
        if key in kwargs.keys():
            self.exponent = kwargs[key]
            pass
        if self.exponent is None:
            print("Key 'fractal_dim' or 'exponent' not found!")
        pass

    def k_th_moment(self):
        M_frac = 0

        # segment_lengths = self.length_list[self.flag_list]
        for i in range(len(self.flag_list)):
            if self.flag_list[i]:
                M_frac += self.length_list[i] ** self.exponent
        # for ll in segment_lengths:
        #     M_frac += ll**self.fractal_dim

        return M_frac

    def run(self, time_iteration, min_iteration, iteration_step):
        """

        :param time_iteration: Total time iteration
        :param min_iteration:  when to start recording
        :param iteration_step: interval after which it will be recorded
        :return:
        """
        # lengths = [1.]
        # flags = [True]
        # frag_prob = [1.]  # raw probability, not normalized
        # frag_prob_sum = 1.0  # normalization const
        #
        # iteration_list = list(range(min_iteration, time_iteration + 1, iteration_step))

        M_realization = []
        for i in range(time_iteration + 1):
            self.one_time_step()
            if (i > min_iteration) and (i % iteration_step == 0):
                M_frac = self.k_th_moment()
            # if i + 1 in iteration_list:
            #     M_frac = fractal_length(lengths, flags)
                M_realization.append(M_frac)
            pass

        return np.array(M_realization)


    def run_ensemble(self, ensemble_size, time_iteration, start_at, step_interval):
        """

        :param ensemble_size: ensemble size
        :param time_iteration: total number of time steps
        :param start_at:       minimum number of step before we start recording data
        :param step_interval: number of step between successive data record
        """
        M_ensemble = None

        step=int(ensemble_size/100)
        start_time = time.time()
        for i in range(ensemble_size):
            if i % step == 0 and self.logging:
                print("realization ", i, " . Time spent ", (start_time - time.time()), " sec")
                start_time = time.time()
                pass
            self.reset()
            M_list = self.run(time_iteration, start_at, step_interval)
            if M_ensemble is None:
                M_ensemble = M_list
            else:
                M_ensemble += M_list
            pass

        M_average = M_ensemble / ensemble_size

        print(M_average)
        return M_average

    pass


# to get the lengths after n iteration
class TrueLengths(StochasticFragmentation):
    """

    """
    def __init__(self, **kwargs):
        super(TrueLengths, self).__init__(**kwargs)
        print("Turning on logging")
        self.log(True)

    def get_signature(self):
        a = super().get_signature()
        return a + "_Lengths_"

    def save_to_file(self, directory, header_description):
        filename = self.get_filename()
        # header = self.get_header(header_description)
        # date_time = header['date_time']
        # filename += date_time
        # filename += ".txt"
        np.savetxt(directory + filename, self.lengths_ensemble, header=self.header_str)

    def run(self, time_iteration):
        self.time_iteration = time_iteration
        for i in range(time_iteration + 1):
            self.one_time_step()
            pass

        lengths = np.array(self.length_list)
        true_lengths = lengths[self.flag_list]

        return true_lengths

    def run_ensemble(self, ensemble_size, time_iteration):
        """

        :param ensemble_size:
        :param time_iteration:
        :return:
        """
        self.ensemble_size = ensemble_size
        self.lengths_ensemble = np.array([])
        step=int(ensemble_size/100)
        start_time = time.time()
        for i in range(ensemble_size):
            if i % step == 0 and self.logging:
                print("realization ", i, " . Time spent ", (start_time - time.time()), " sec")
                start_time = time.time()
                pass
            self.reset()
            length = self.run(time_iteration)
            self.lengths_ensemble = np.append(self.lengths_ensemble, length)
            pass

        return self.lengths_ensemble
    pass