#!bin/bash/

# set range of parameters for z*, r*, th' (3 values per index)

# read in parameters
read -p "dTh: " dTh;
read -p "nx: " nx;
read -p "dx: " dx;


#path="/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m/"
path="/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run3_dx50m/"
casename="ColdPoolDry_single_3D"

# set geometry parameters
if [ $dTh -eq 1 ]
then
  z_params=( 3465 1730 1155 )  #run1
  r_params=( 1155 1730 3465 )  #run1
elif [ $dTh -eq 2 ]
then
  #z_params=( 2450 1225 815 ) #run1
  #r_params=( 815 1225 2450 ) #run1 
  z_params=( 500 900 1600 1900 2500 ) #run2
  r_params=( 1900 1300 900 800 600 )  #run2
  z_params=( 1600 )
  r_params=( 900 )
elif [ $dTh -eq 3 ]
then
  #z_params=( 2000 250 500 670 1000 1500 2000 4000 ) #run1
  #r_params=( 2000 4000 2000 1500 1000 670 500 250 ) #run1
  #z_params=( 500 1000 1600 2000 2500 ) #run2
  #r_params=( 1500 1000 700 600 500 )   #run2
  z_params=( 500 1000 2000 ) #run3
  r_params=( 1500 1000 600 ) #run3
  z_params=( 1000 ) #run4
  r_params=( 1000 ) #run4
elif [ $dTh -eq 4 ]
then
  #z_params=( 430 870 1730 ) #run1
  #r_params=( 1730 870 430 ) #run1
  z_params=( 500 900 1600 2000 2500 ) #run2
  r_params=( 1300 900 600 500 400 )   #run2
  #z_params=( 2500 )
  #r_params=( 400 )
fi





n_geom=${#z_params[@]}
n_therm=${#th_params[@]}
n_tot=$(( $n_geom*$n_therm ))
#echo "dTh:" $dTh
echo "z-parameters:" ${z_params[@]} 
echo "r-parameters:" ${r_params[@]}

echo "#geometry parameters:" $n_geom




echo " "
echo "---- start loop ----"
echo " "

count_geom=0
k=0

  while [ $count_geom -lt $n_geom ]
  do
    zstar=${z_params[$count_geom]}
    rstar=${r_params[$count_geom]}
    echo "-- parameters:" $zstar $rstar
  
    id="dTh"$dTh"_z"$zstar"_r"$rstar
    echo $id
  
    fullpath=$path$id
    echo $fullpath

    echo "--- run tracers (k="$k") ---"
    ./run_raintrack_loop.job "$fullpath" $k $rstar $nx $dx

    #echo "run tracers"
    #./run_raintrack_bet.job
    #./run_raintrack_bet.job "$fullpath"
  
    echo " "
    ((count_geom++))
  done


echo "finished bash script"

