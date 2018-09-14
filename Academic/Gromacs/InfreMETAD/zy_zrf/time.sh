time.bash
################
filename=$1

kb=2.476;  #kt at room temperature

dtmd=1000; # #time in femtosec between two lines of COLVAR

temperature=300; # temperature in PLUMED


dir=$PWD

awk '(!/FIELDS/){print}' $filename | awk -v dt=$dtmd -v kbt=$kb  -v T=$temperature '{

      count=count+1;

        sum=sum+dt*exp(300/kbt/T*($4));

        sum2=sum2+exp(300/kbt/T*($4));

        alpha=sum/(dt*count)

        de_alpha=( sum2 - sum/(dt*count) ) / (dt*count)

      printf("Reweighted Unbinding time = %8.3f ms | Metadynamics time = %6.2f ns | alpha = %7.2e | %s \n", sum/1e12, dt*count/1e6, alpha, dir)

     }' | tail -1