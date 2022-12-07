import os, sys
repo_path = f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}"
sys.path.append(repo_path)
from enhanced_sampling.schema import SimParams

for file in os.listdir("../sim_params"):
    sim_params = SimParams(os.path.join("../sim_params", file))
    print(file)
    print(sim_params.get_time_scales(), '\n')