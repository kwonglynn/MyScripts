#Make COLVAR have the same length as the trajectory.

# Step 1, remove the comment lines after the third line:
fi = open('COLVAR', 'r')
fo = open('COLVAR-no-overlap', 'w')

i = 0
time0 = -1
for line in fi:
    # Only write the first comment lines, omit the commet lines for restarting runs.
    if line.startswith('#') and i < 3:
        fo.write(line)
        i += 1
    elif not line.startswith('#'):
        time = int(line.strip().split()[0].split('.')[0])
        # Coninuous output
        if time > time0:
            time0 = time
            fo.write(line)
        # Redundant output, ignore it.
        elif time < time0:
            continue
        
fi.close()
fo.close()

# Step 2, Extract the data every 10th line:
fi = open('COLVAR-no-overlap', 'r')
fo = open('COLVAR-100', 'w')

lines = fi.readlines()
i = 0
while i < len(lines):
    if lines[i].startswith('#'):
        fo.write(lines[i])
        i += 1
    else:
        try:
            fo.write(lines[i])
            i += 100
        except IndexError:
            break

fi.close()
fo.close()
