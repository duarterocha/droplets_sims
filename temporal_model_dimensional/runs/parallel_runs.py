from stationary_model.utils import *
from pyoomph.expressions.units import *
from pyoomph.utils.paramscan import *

if __name__=="__main__":

    #Create a parameter scanner, give the script to run and the max number of processes to run simultaneously
    scanner=ParallelParameterScan("single_run.py", max_procs=6)
    Strength_param_range = parameter_range(numperoder=20, log_scale=False, Strength=1)["Strength"]

    for strength_value in Strength_param_range: #Scan the spring constant
        sim=scanner.new_sim("diffusivity_"+str(strength_value))

        default_surfactant_diffusivity = 1e-10

        #Modify the parameters (strength_value ranges from 0 to 0.02; with this, the surfactants diffusivity will range
        # from 1e-10 to 1e-6)
        sim.surfactants_diffusivity = default_surfactant_diffusivity + strength_value/0.02*(1e-6-default_surfactant_diffusivity)
        sim.surfactants_diffusivity *= meter**2 / second

    sim.surfactants_diffusivity = 1e-9 * meter ** 2 / second
    for strength_value in Strength_param_range:  # Scan the spring constant
        sim = scanner.new_sim("Ma_G_" + str(strength_value))

        # Modify the parameters (Ma_G ranging from 0 to 0.02)
        sim.Ma_G = strength_value

    #Run all (and rerun also already finished sims)
    scanner.run_all(skip_done=True)

