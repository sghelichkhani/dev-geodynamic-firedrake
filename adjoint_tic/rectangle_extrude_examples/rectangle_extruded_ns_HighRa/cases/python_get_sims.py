#!/bin/python3
import numpy as np
import subprocess
import os 
qstat_out = np.array(subprocess.run("qstat", capture_output=True).stdout.decode('utf-8').split('\n'))
qstat_out = qstat_out[['gadi-pbs' in i for i in qstat_out]]
JobIds   = [i.split()[0] for i in qstat_out] 


# Getting all the names
JobSpits = [subprocess.run(["qstat","-f", i.split()[0]], capture_output=True).stdout.decode('utf-8') for i in qstat_out]
pths = [JobSpit[JobSpit.find("Output_Path ")+len("Output_Path"): JobSpit.find("Priority")].replace(" ","").replace('\t','').replace("\n","") for JobSpit in JobSpits] 
pths = [pth[pth.find('/'):] for pth in pths]
for i, jid in enumerate(JobIds):
    with open(pths[i], mode='w') as fout: 
        subprocess.run(["qcat", jid], stdout=fout)
    
