# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

d = pd.read_csv('runs.csv')
score_india = d['score_india']
legend = ['India', 'Pakistan']
score_pk = d['score_pk']
plt.hist([score_india, score_pk], color=['orange', 'green'])
plt.xlabel("Runs/Delivery")
plt.ylabel("Frequency")
plt.legend(legend)
plt.xticks(range(0, 7))
plt.yticks(range(1, 20))
plt.title('Champions Trophy 2017 Final\n Runs scored in 3 overs')
plt.show()