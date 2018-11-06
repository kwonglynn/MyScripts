#eight replicas
nrep=10
# "effective" temperature range
tmin=310
tmax=600

# build geometric progression
list=$(
awk -v n=$nrep \
    -v tmin=$tmin \
    -v tmax=$tmax \
  'BEGIN{for(i=0;i<n;i++){
    t=tmin*exp(i*log(tmax/tmin)/(n-1));
    printf(t); if(i<n-1)printf(",");
  }
}'
)

# clean directory
#rm -fr \#*
#rm -fr topol*

for((i=0;i<nrep;i++))
do

# choose lambda as T[0]/T[i]
# remember that high temperature is equivalent to low lambda
  lambda=$(echo $list | awk 'BEGIN{FS=",";}{print $1/$'$((i+1))';}') 
# process topology # (if you are curious, try "diff topol0.top topol1.top" to see the changes)
  plumed partial_tempering $lambda < processed.top > topol$i.top 
# prepare tpr file # -maxwarn is often needed because box could be charged  
  gmx_mpi grompp -c npt.gro -n index.ndx  -maxwarn 1 -o topol$i.tpr -f rest2_membrane.mdp -p topol$i.top
done
