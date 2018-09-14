for i in {2..13}
do 
cd trp_${i}
#cd leu

kb=2.496;  #kt at room temperature

dtmd=1000; # #time in femtosec between two lines of COLVAR 

temperature=300; # temperature in PLUMED


dir=$PWD
awk '(!/FIELDS/){print}' COLVAR | awk -v dt=$dtmd -v kbt=$kb  -v T=$temperature '{

      count=count+1;

        sum=sum+dt*exp(300/kbt/T*($5));  

        sum2=sum2+exp(300/kbt/T*($5));  

        alpha=sum/(dt*count)

        de_alpha=( sum2 - sum/(dt*count) ) / (dt*count)

      printf("Reweighted Unbinding time = %8.3f ms | Metadynamics time = %6.2f ns | alpha = %7.2e | %s \n", sum/1e12, dt*count/1e6, alpha, dir)

     }' | tail -1
cd ../
done

