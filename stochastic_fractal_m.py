
#from main_py import stochastic_fragmentation_v1
from main_py import stochastic_fragmentation_v2

#from main_py.stochastic_fragmentation_v1 import *
from main_py import analytic_solution
from main_py.analytic_solution import AnalyticSoln
from main_py.stochastic_fragmentation_v2 import  *

StochasticFragmentation = StochasticFragmentation_v2
Moment = Moment_v2
NumberLength = NumberLength_v2
TrueLengths = TrueLengths_v2

FractalLength = Moment
# find_df = df_determination.find_df


## Creating temp directory
import os
repo_name="Stochastic-Fractal"
directory=os.getcwd().split(repo_name)
ROOT_directory = directory[0] + repo_name

dir_data = ROOT_directory + "/data/temp"
dir_fig  = ROOT_directory + "/figures/temp"

for dir_name in [dir_data, dir_fig]:
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
        print("directory ", dir_name, " is created")
        pass
    pass
pass

