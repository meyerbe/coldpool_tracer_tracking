#!bin/bash/

# set range of parameters for z*, r*, th' (3 values per index)

# read in parameters
read -p "dTh: " dTh; 
#read -p "tmin: " tmin; 
#read -p "tmax: " tmax; 


path="/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/"
casename="ColdPoolDry_single_3D"

# set geometry parameters
z_params=( 1600 )
r_params=( 700 )
nx=( 400 )
dx=100

n_geom=${#r_params[@]}
#n_therm=${#th_params[@]}
#n_tot=$(( $n_geom*$n_therm ))
echo "dTh:" $dTh
echo "z-parameters:" ${z_params[@]} 
echo "r-parameters:" ${r_params[@]}
echo "#geometry parameters:" $n_geom
echo " "

count_geom=0
k=0

  while [ $count_geom -lt $n_geom ]
  do
    zstar=${z_params[$count_geom]}
    rstar=${r_params[$count_geom]}
    echo "parameters:" $zstar $rstar
  
    id="dTh"$dTh"_z"$zstar"_r"$rstar
    echo $id
  
    fullpath=$path$id
    echo $fullpath
    echo " "

    echo "run tracers (k="$k")"
    ./run_raintrack_loop.job "$fullpath" $k $rstar ${nx[$count_geom]} $dx $interpol 0 

    #echo "run tracers"
    #./run_raintrack_bet.job
    #./run_raintrack_bet.job "$fullpath"
  
    echo " "
    ((count_geom++))
  done


echo "finished bash script"

