# control parameters defined in functions, taken from either figures or equations in
# https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.16.024060

import numpy as np

def Omega(t, eta=1, tau=500):
    if eta == 0:
        sigma = 2000
        return np.exp(-1 * (t - 125)**2 / sigma) * 12 +  np.exp(-1 * (t - 375)**2 / sigma) * 12
    else:
        print(f'Eta = {eta} not implemented yet')
        
def phase1(t, tau=500):
    return np.pi / 2