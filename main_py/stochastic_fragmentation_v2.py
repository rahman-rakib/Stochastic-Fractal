import random
import numpy as np
import datetime
import json
import time


class StochasticFragmentation_v2:
    """
    Simulation for Stochastic Binary Fragmentation
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
        self.frac_sum = 1.0  # normalization constant
        self.frac_prob = [1.0]
        self.length_list = [1.0]
        self.choose_pivot = self.betadist

        # how many times
        self.ensemble_size = 1 # at least once
        self.time_iteration = 0 # it will be set when run() executes

        # loggin
        self.logging = False
        pass

    def get_alpha(self):
        return self.alpha

    def get_probability(self):
        return self.prob

    def get_df(self):
        from main_py import analytic_solution
        ana_soln = analytic_solution.AnalyticSoln()
        return ana_soln.find_df(self.alpha, self.prob)

    def get_time_iteration(self):
        return self.time_iteration

    def log(self, flag=False):
        """
        to enable logging during simulation. Useful to monitor the progress of the simulation
        :param flag:
        :return:
        """
        self.logging = flag

    def get_signature(self):
        """
        signature for automatically generated data file.
        :return:
        """
        sig = "StochasticFragmentation"
        sig += "_alpha_{}".format(self.alpha)
        sig += "_p_{}".format(self.prob)
        return sig

    def reset(self):
        """
        After each realization the state of the system must be reset
        :return:
        """
        self.frac_sum = 1.0  # normalization constant TODO
        self.frac_prob = [1.0]
        self.length_list = [1.0]

    def betadist(self):
        """gives a random number from beta distribution"""
        #         print("betadist")
        return random.betavariate(self.alpha, self.alpha)

    def decision(self):
        """
        decides with a given probability whether to keep the right part
        """
        if self.prob > random.random():
            return 1.
        else:
            return 0.

    def splitting(self, segment):
        """
        splits a given segment. left and right are endpoints of the segment
        segment : length of a segment
        returns :
            xL -> length of the left segment
            xR -> length of the right segment
            flag -> keeping the right segment
            xLp, xRp -> probability of being selected
            change -> change of normalization const
        """
        xL = segment * self.choose_pivot()

        flag = self.decision()
        xR = (segment - xL) * flag
        xLp = xL ** self.expon
        xRp = xR ** self.expon
        change = xLp + xRp - segment ** self.expon
        return xL, xR, flag, xLp, xRp, change

    def pickindex(self):
        """
        picks up a segment to be subsequently split
        """
        # print("pickindex")
        # print("self.probability_list ", self.probability_list)
        r = random.uniform(0, 1.0)  # 1 is the initial particle size
        if r > self.frac_sum:
            # print("r > self.normC => return None")
            return None
        sum_ = 0
        for index in range(len(self.frac_prob)):
            sum_ += self.frac_prob[index]
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
            print("{:.5e}, {:.5e}".format(self.length_list[i],
                                                self.frac_prob[i] / self.frac_sum
                                                ))
            pass
        print("sum of all lengths : ", np.sum(self.length_list))
        pass

    def one_time_step(self):
        """
        One step of each time iteration
        :return:
        """
        # print("StochasticFragmentation.one_time_step")
        index = self.pickindex()
        # print("index ", index)
        if index is not None:
            xL, xR, flag, xLp, xRp, change = self.splitting(self.length_list[index])

            self.length_list[index] = xL
            self.frac_prob[index] = xLp
            if flag:
                # self.flag_list.append(flag)
                self.length_list.append(xR)
                self.frac_prob.append(xRp)
                pass
            self.frac_sum += change
            pass

        pass

    def get_header(self, header_description="NA"):
        """
        basic header information to write to the data file.
        First line of the data fill generated and saved from this program will
        contatin a json header.
        :param header_description:
        :return:
        """
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
        """
        Multiple run of smaller ensemble to get a large ensemble is a good way of saving time,
        but the data file name of each independent run must be different for this to work.
        :return:
        """
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
class NumberLength_v2(StochasticFragmentation_v2):
    """
    In order to find the fractal dimension
    """
    def get_signature(self):
        return super(NumberLength_v2, self).get_signature() + "_NumberLength"

    def number_length(self):
        """

        :return: N, M
            N -> 0-th moment or number of surviving segments
            M -> 1st moment or sum of sizes of `N` segments
        """
        # lengths = np.array(self.length_list)
        # lengths = lengths[self.flag_list]
        # segment_count = lengths.shape[0]
        segment_count = len(self.length_list)
        surviving_length_sum = np.sum(self.length_list)

        return segment_count, surviving_length_sum

    def run(self, time_iteration, min_iteration, iteration_step):
        """
        One realization.
        we run the `one_time_step()` method `time_iteration` times. We record some information
        at each `iteration_step` step interval starting from `min_iteration`
        :param time_iteration: maximum time step
        :param min_iteration: starting point to record data
        :param iteration_step: number of data points we want to generate.
        :return: a numpy array with two columns,
                1st column is the particle counts and
                2nd column is the sum of surviving particle sizes
        """
        # iteration_list = list(range(min_iteration, time_iteration + 1, iteration_step))
        N_realization = []
        M_realization = []
        # step_size = self.get_step_size(time_iteration, min_iteration, number_of_points)
        # print("step size = ", step_size)
        for i in range(1, time_iteration + 1):
            #             print("time step ", i)
            self.one_time_step()
            if (i >= min_iteration) and (i % iteration_step == 0):
                # print("t = ", i)
                segment_count, surviving_length_sum = self.number_length()
                N_realization.append(segment_count)
                M_realization.append(surviving_length_sum)
            pass

        N_list = np.array(N_realization)
        M_list = np.array(M_realization)
        self.time_iteration = time_iteration
        return np.c_[N_list, M_list]

    def run_ensemble(self, ensemble_size, time_iteration, start_at, iteration_step):
        """

        :param ensemble_size: ensemble size
        :param time_iteration: total number of time steps
        :param start_at:       minimum number of step before we start recording data
        :param iteration_step: starting from `start_at` it will take `number_of_data_points`
                data points upto time `time_iteration`
        :return: [N, M] average value
                where, N = number of particles at particular time steps
                       M = sum of surviving length at those time steps
        """

        ensemble_data = None
        step_temp = int(ensemble_size / 100)
        if step_temp == 0:
            step_temp += 1
            pass
        start_time = time.time()
        interval_time = time.time()
        for i in range(1, ensemble_size + 1):

            self.reset()
            out_data_N_M = self.run(time_iteration, start_at, iteration_step)
            if ensemble_data is None:
                ensemble_data = out_data_N_M
            else:
                ensemble_data += out_data_N_M
                pass
            if i % step_temp == 0 and self.logging:
                print("realization ", i, " . Time spent ", (time.time() - interval_time), " sec")
                interval_time = time.time()
                pass
            pass

        data_average = ensemble_data / ensemble_size

        self.time_iteration = time_iteration
        self.ensemble_size = ensemble_size
        print("Total time spent ", (time.time() - start_time), " sec")
        return data_average

    pass


# xdf
class Moment_v2(StochasticFragmentation_v2):
    """
    finding n-th moment.
    df-th moment is always conserved.
    """
    def __init__(self, **kwargs):
        super(Moment_v2, self).__init__(**kwargs)
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
            self.exponent = self.get_df()
            print("Key 'fractal_dim' or 'exponent' not found!")
            print("Theoretical value is used. df = ", self.exponent)
        pass

    def get_signature(self):
        sig = super(Moment_v2, self).get_signature()
        sig += "_Moment"
        return sig

    def k_th_moment(self):
        """
        returns the k-th moment for the current state of the system. Here k=`self.exponent`
        :return:
        """
        M_frac = 0
        for ll in self.length_list:
                M_frac += ll ** self.exponent
                pass
        return M_frac

    def run(self, time_iteration, start_at, iteration_step):
        """
        one realization.
        :param time_iteration: Total time iteration
        :param start_at:  when to start recording
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
            if (i >= start_at) and (i % iteration_step == 0):
                M_frac = self.k_th_moment()
            # if i + 1 in iteration_list:
            #     M_frac = fractal_length(lengths, flags)
                M_realization.append(M_frac)
            pass

        return np.array(M_realization)

    def run_ensemble(self, ensemble_size, time_iteration, start_at, step_interval):
        """
        A number of realization
        :param ensemble_size: ensemble size
        :param time_iteration: total number of time steps
        :param start_at:       minimum number of step before we start recording data
        :param step_interval: number of step between successive data record
        """
        M_ensemble = None

        step_temp=int(ensemble_size/100)
        if step_temp == 0:
            step_temp += 1
            pass
        start_time = time.time()
        interval_time = time.time()
        for i in range(1, ensemble_size+1):

            self.reset()
            M_list = self.run(time_iteration, start_at, step_interval)
            if M_ensemble is None:
                M_ensemble = M_list
            else:
                M_ensemble += M_list
                pass
            if i % step_temp == 0 and self.logging:
                print("realization ", i, " . Time spent ", (time.time() - interval_time), " sec")
                interval_time = time.time()
                pass
            pass

        M_average = M_ensemble / ensemble_size
        print("Total time spent ", (time.time() - start_time), " sec")
        # print(M_average)
        return M_average

    pass


# to get the lengths after n iteration
class TrueLengths_v2(StochasticFragmentation_v2):
    """
    Segment lengths are appended to a list after `T` iteration in each realization.
    Histogram is generated from that list of segment lengths.
    """
    def __init__(self, **kwargs):
        super(TrueLengths_v2, self).__init__(**kwargs)
        print("Turning on logging")
        self.log(True)
        self.lengths_ensemble = None

    def get_signature(self):
        a = super().get_signature()
        return a + "_Lengths"

    def save_to_file(self, directory, header_description='NA'):
        """
        To save the simulation data.
        :param directory:
        :param header_description:
        :return:
        """
        if self.lengths_ensemble is None:
            print("No data to write to a file")
            return None
        filename = self.get_filename()
        header = self.get_header(header_description)
        # date_time = header['date_time']
        # filename += date_time
        # filename += ".txt"
        np.savetxt(directory + filename, self.lengths_ensemble, header=self.header_str)

    def run(self, time_iteration):
        """
        One realization.
        :param time_iteration: maximum time iteration
        :return:
        """
        self.time_iteration = time_iteration
        for i in range(time_iteration + 1):
            self.one_time_step()
            pass

        lengths = np.array(self.length_list)
        return lengths

    def run_ensemble(self, ensemble_size, time_iteration, thread_id=-1, return_dict=None):
        """
        :param ensemble_size: number of independent realization
        :param time_iteration: maximum time iteration
        :return:
        """
        self.ensemble_size = ensemble_size
        self.lengths_ensemble = np.array([])
        step_temp=int(ensemble_size/100)
        if step_temp == 0:
            step_temp = 1
            pass
        start_time = time.time()
        interval_time = time.time()
        for i in range(1, ensemble_size+1):

            self.reset()
            length = self.run(time_iteration)
            self.lengths_ensemble = np.append(self.lengths_ensemble, length)
            if i % step_temp == 0 and self.logging:

                if thread_id >= 0:
                    # print("t[", thread_id, "]",  "realization ", i, " . Time spent ", (time.time() - interval_time), " sec")
                    pass
                else:
                    print("realization ", i, " . Time spent ", (time.time() - interval_time), " sec")
                interval_time = time.time()
                pass
            pass
        print("Total time spent ", (time.time() - start_time), " sec")
        if thread_id >=0:
            return_dict[thread_id] = self.lengths_ensemble
            pass
        else:
            return self.lengths_ensemble

    def run_ensemble_parallel(self, ensemble_size, time_iteration, threads=2):
        import multiprocessing as mp

        thread_counts = mp.cpu_count()
        if threads > thread_counts:
            thread_counts = mp.cpu_count()
        else:
            thread_counts = threads
            pass

        if (ensemble_size < 1000) and (time_iteration < 10000):
            print("Maybe you don't need multi threading")
            pass
        EN_per_thread = ensemble_size // thread_counts

        manager = mp.Manager()
        return_dict = manager.dict()

        all_processes = []
        print("thread_counts ", thread_counts)
        for i in range(thread_counts):
            # inargs = {'length':LL, "ensembleSize":En, "thread_count":i}
            process = mp.Process(target=self.run_ensemble,
                                              args=(EN_per_thread, time_iteration, i, return_dict))
            all_processes.append(process)
            process.start()

            pass
        for process in all_processes:
            process.join()
            pass

        length_list = []
        for key in return_dict.keys():
            tolistsss = return_dict[key].tolist()
            # print(len(tolistsss))
            length_list += tolistsss
            # print(length_list)
            # print("length_list ", len(length_list))
            pass
        # length_list = np.array(return_dict.values())
        return length_list

    pass


if __name__ == '__main__':
    stochastic_frag = TrueLengths_v2(alpha=3, probability=0.75)
    ens_data_dct = stochastic_frag.run_ensemble_parallel(100, 1000, 4)  # from class
    print(len(ens_data_dct))
