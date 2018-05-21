for i in {4..20}; do echo Node${i} >> performance.dat; grep 'Performance:' Node${i}/md1.log >> performance.dat; done
