#!/bin/bash

# names for simulations
SuffVals=(1e-1 1e-2 1e-3 1e-4)
dirname=("Brent" )
cur_pth=$(pwd);
source ~/firedrake/bin/activate
export PYTHONPATH='/home/sia/local_src/lib/python3.8/site-packages/':${PYTHONPATH}

for ((k=0; k<"${#dirname[@]}"; k++)); do
   for i in ${SuffVals[@]};do
      dir=$(printf "%s" BFGS_${dirname[$k]}_Suff_${i});
      mkdir -pv ${dir}; 
      cp -v adj_reconst_ROL_fromprms_cache.py ${dir};
      echo ${mymethod}
      cd ${dir}
      #sed -i '265s/.*/                                "Type":  "${mymethod}",/' adj_reconst_ROL_fromprms_cache.py ;
      sed -i '284s/.*/                "Sufficient Decrease Tolerance": '$(printf "%4.4f" ${i})',/' adj_reconst_ROL_fromprms_cache.py ;
      mpiexec -np 16 python3 adj_reconst_ROL_fromprms_cache.py &> ${dir}.out &
      cd ..;
   done
done
