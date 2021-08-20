#!/bin/bash

methods=("Brent's" "Cubic Interpolation")
methods_name=("Brent" "Cubic")
sims=(0.1 0.01 0.001 0.0001)

for ((i=0; i<${#methods[@]}; i++));do
   for j in ${sims[@]};do
      dir=$(printf "BFGS_%s_%s" "${methods_name[i]}" ${j})
      mkdir -vp ${dir};
      cp -v template.py ${dir}/adj_reconst_ROL_fromprms_cache.py;
      cd ${dir}
      echo "${methods[$i]}"
      sed -i'.original' -e '268s|.*|                                "Type":  "'"${methods[$i]}"'",|' adj_reconst_ROL_fromprms_cache.py
      sed -i'.original' -e "265s/.*/                'Sufficient Decrease Tolerance': $(printf "%s" ${j}),/" adj_reconst_ROL_fromprms_cache.py
      rm *original
      echo ${dir}
      cd .. 
   done
done

