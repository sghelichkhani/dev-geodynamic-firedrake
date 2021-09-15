import os
import numpy as np
import subprocess
import itertools as it

# number of processes to be used
num_proc = 16


#selected_prms = {
#	'SUFFICIENTDECREASETOLERANCE': [0.1, 0.01, 0.001, 0.0001],
#    'STORAGE_BFGS' : [5, 10, 20],
#    'USEPREVIOSBOOL' : [True, False],
#    'LINESEARCHMETHOD' : ["\"Brent's\"", "\"Cubic Interpolation\"", "\"Backtracking\"", "\"Bisection\""],
#        }
#selected_prms = {
#	'SUFFICIENTDECREASETOLERANCE': [0.1, 0.01, 0.001],
#    'STORAGE_BFGS' : [5, 10],
#    'USEPREVIOSBOOL' : [True, False],
#    'LINESEARCHMETHOD' : ["\"Brent's\""],
#        }
selected_prms = {
        'STORAGE_BFGS' : [5, 10, 20, 50, 100],
        }
name_pattern="test_TR"

template_name = "template_prms_TR.py"
all_prms = sorted(selected_prms)
combinations = it.product(*(selected_prms[Name] for Name in all_prms))
for k , combi in enumerate(combinations):
    dir_name = str("%s_%3.3i" %(name_pattern, k))
    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)

    # Writing out the python scripts
    with open(os.path.join(dir_name, "./optimisation.py"), mode='w') as fi_out:
        with open(template_name, mode='r') as fi_in:
            for line in fi_in:
                temp_line = line
                for i, key in enumerate(all_prms):
                    if key in temp_line:
                        temp_line = temp_line.replace(key, str(combi[i]))
                fi_out.write(temp_line)

    # Writing out qsub file 
    with open(os.path.join(dir_name, "qsub.sh"), mode='w') as q_fi_out:
        with open('./qsub.sh', mode='r') as q_fi_in:
            for line in q_fi_in:
                if "SIMUNAME" in line:
                    q_fi_out.write(line.replace("SIMUNAME", dir_name))
                elif "ALL_PRMS" in line:
                    q_fi_out.write(line.replace("ALL_PRMS", " ".join(np.array(all_prms + list(combi)))))
                else: 
                    q_fi_out.write(line)

