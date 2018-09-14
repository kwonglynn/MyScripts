# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 12:00:08 2018

@author: Guanglin Kuang
"""

def main():
    KbT0 = 2.479            # KbT at 298 K, in kJ/mol
    KbT = T / 298 * KbT0    # Calculate the KbT at the simulation temperature
    print KbT
 
if __name__ == "__main__":
    N = 20                  # Number of independent simulations.
    T = 300                 # The temperature of the simulation in K
    main()