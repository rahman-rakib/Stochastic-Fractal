from main_py.stochastic_fragmentation import *
from main_py import df_determination
from main_py import analytic_solution
from main_py.analytic_solution import AnalyticSoln
import os

FractalLength = Moment
find_df = df_determination.find_df


from pathlib import Path
file = "stochastic_fractal_m.py"
# ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = Path(os.path.dirname(os.path.abspath(file))).parent
print("ROOT_DIR ", ROOT_DIR)