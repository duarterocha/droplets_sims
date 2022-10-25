from stationary_model.utils import *
from pyoomph.utils.paramscan import *

if __name__=="__main__":

    #Create a parameter scanner, give the script to run and the max number of processes to run simultaneously
    scanner=ParallelParameterScan("param_scan.py",max_procs=6)
    Strength_param_range = parameter_range(Strength=4)["Strength"]

    for strength_value in Strength_param_range: #Scan the spring constant
        sim=scanner.new_sim("param_scan_run_"+str(round(strength_value)))

        #Modify the parameters
        sim.param_Strength=strength_value

    #Run all (and rerun also already finished sims)
    scanner.run_all(skip_done=True)

