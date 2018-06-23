#1. Judge if a file exists and do the processing, used in umbrella sampling simulations:
for i in {0..29}: do
    if [ -f complex_$i.gro ]; then
        mkdir Umbrella$i
        cp complex${i}.gro common/* Umbrella$i
        cd Umbrella$i
        cp complex${i}.gro complex.gro
        sh gromacs_common.sh
        cd ..
    fi
done

#2. 