# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
x = np.random.random_integers(1, 100, 5)
plt.hist(x, bins=20)
plt.ylabel('No of times')
plt.show()