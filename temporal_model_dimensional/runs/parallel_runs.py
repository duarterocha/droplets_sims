from stationary_model.utils import *
from pyoomph.expressions.units import *
from pyoomph.utils.paramscan import *

if __name__=="__main__":

    #Create a parameter scanner, give the script to run and the max number of processes to run simultaneously
    scanner=ParallelParameterScan("single_run.py", max_procs=6)
    Strength_param_range = parameter_range(numperoder=20, log_scale=False, Strength=1)["Strength"]

    for strength_value in Strength_param_range: #Scan the spring constant
        sim=scanner.new_sim("Ma_G_"+str(round(strength_value)))

        default_surfactants_diffusivity = 1e-9 * meter ** 2 / second

        #Modify the parameters
        sim.Ma_G = strength_value

    #Run all (and rerun also already finished sims)
    scanner.run_all(skip_done=True)

