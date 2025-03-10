# import gplugins.luminescent as gl
from gdsfactory.generic_tech import *
from .runs_utils import *
from .utils import *
from .constants import *
from .materials import *
from .sol import *

from .pic.gcells import *
from .pic.setup import *
from .pic.inverse_design import *
from .pic.sparams import *

from .gen.setup import *
print("imported luminescent")


# __all__ = [
#     "pic_design_problem",
#     "sparams_problem",
#     "solve",
#     # "apply_design",
#     "load_solution",
#     "finetune",
#     "load_problem",
#     # "make_training_movie",
#     # "make_simulation_movie",
#     "mimo",
# ]
