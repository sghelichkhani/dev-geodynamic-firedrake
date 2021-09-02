#!/bin/python3
import numpy as np
import subprocess
import os 
qstat_out = np.array(subprocess.run("qstat", capture_output=True).stdout.decode('utf-8').split('\n'))
qstat_out = qstat_out[['gadi-pbs' in i for i in qstat_out]]
JobIds   = [i.split()[0] for i in qstat_out] 


# Getting all the names
JobNames = [subprocess.run(["qstat","-f", i.split()[0]], capture_output=True).stdout.decode('utf-8').split('\n')[1].split('=')[-1].replace(' ','') for i in qstat_out]


for i, jid in enumerate(JobIds):
    os.chdir(JobNames[i])
    with open(JobNames[i]+'.o'+jid[0:jid.rfind('.')], mode='w') as fout: 
        subprocess.run(["qcat", jid], stdout=fout)
    os.chdir('..')
    
