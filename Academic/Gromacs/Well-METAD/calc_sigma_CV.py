import numpy as np

fi = open('COLVAR', 'r')

CV1 = []
CV2 = []

for line in fi:
    if not line.startswith('#'):
        terms = line.split()
        CV1.append(float(terms[1]))
        CV2.append(float(terms[2]))

print("The standard deviation of CV1 is:\t%5.2f" % np.std(CV1))
print("The standard deviation of CV2 is:\t%5.2f" % np.std(CV2))

fi.close()
