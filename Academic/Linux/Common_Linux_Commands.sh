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

#2. Change the permissions of a folder and file:
# Reference: https://stackoverflow.com/questions/3740152/how-do-i-set-chmod-for-a-folder-and-all-of-its-subfolders-and-files
# To change all the directories to 755 (drwxr-xr-x):
find /opt/lampp/htdocs -type d -exec chmod 755 {} \;
# To change all the files to 644 (-rw-r--r--):
find /opt/lampp/htdocs -type f -exec chmod 644 {} \;