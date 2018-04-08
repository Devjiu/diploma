import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt

def suppose_parabola(parabola_arg, candidates):
    if len(parabola_arg) != 3:
        raise ValueError("Input error: ", parabola_arg)
    # print("Initial arguments: x_1 = %f, x_2 = %f, x_3 = %f" % (parabola_arg[0], parabola_arg[1], parabola_arg[2]))
    
    parabola_coefficients = {'a': 0.0, 'b': 0.0, 'c': 0.0}
    denominator = \
        parabola_arg[0] * np.power(parabola_arg[1], 2) - np.power(parabola_arg[0], 2) * parabola_arg[1] \
        - parabola_arg[0] * np.power(parabola_arg[2], 2) + parabola_arg[1] * np.power(parabola_arg[2], 2) \
        + np.power(parabola_arg[0], 2) * parabola_arg[2] - np.power(parabola_arg[1], 2) * parabola_arg[2] 

    parabola_coefficients['a'] = \
        - parabola_arg[1] * candidates[0] + parabola_arg[2] * candidates[0] \
        + parabola_arg[0] * candidates[1] - parabola_arg[2] * candidates[1] \
        - parabola_arg[0] * candidates[2] + parabola_arg[1] * candidates[2]

    parabola_coefficients['b'] = \
        np.power(parabola_arg[1], 2) * candidates[0] - np.power(parabola_arg[2], 2) * candidates[0] \
        - np.power(parabola_arg[0], 2) * candidates[1] + np.power(parabola_arg[2], 2) * candidates[1] \
        + np.power(parabola_arg[0], 2) * candidates[2] - np.power(parabola_arg[1], 2) * candidates[2]

    parabola_coefficients['c'] = \
        parabola_arg[1] * np.power(parabola_arg[2], 2) * candidates[0] - \
            np.power(parabola_arg[1], 2) * parabola_arg[2] * candidates[0] \
        - parabola_arg[0] * np.power(parabola_arg[2], 2) * candidates[1] + \
            np.power(parabola_arg[0], 2) * parabola_arg[2] * candidates[1] \
        + parabola_arg[0] * np.power(parabola_arg[1], 2) * candidates[2] - \
            np.power(parabola_arg[0], 2) * parabola_arg[1] * candidates[2]
    
    skip = False
    for coef in parabola_coefficients:
        parabola_coefficients[coef] = parabola_coefficients[coef] / denominator
        # if np.isnan(parabola_coefficients[coef]):
        #    skip = True
        #    print("Strange coefs: ", coef, " - ", parabola_coefficients[coef])
        # print("Resulting coefficients ", coef, ": ", parabola_coefficients[coef])
        # else:
        #    skip = True

    # TODO handle NaNs correctly
    if np.NaN in parabola_arg:
        print("nan here")
        ind = np.where(parabola_arg != np.NaN)[0]
        if len(ind) != 1:
            if np.isnan(ind[0]):
                ind = 0
            else:
                ind = ind[0]
        parabola_min_arg = parabola_arg[ind]
    else:
        parabola_min_arg = np.float64( (- parabola_coefficients['b']) / (2 * parabola_coefficients['a']))
        # atoms[num][1] = atoms_old[num][1] + parabola_min_arg
    
    E_new = parabola_coefficients['a'] * np.power(parabola_min_arg, 2) + \
            parabola_coefficients['b'] * parabola_min_arg + parabola_coefficients['c']
    
    # print("Minimum of parabola: ", parabola_min_arg, " Energy: ", E_new)
    parab_minimum = namedtuple("parab_minimum", ["E", "arg", "grad_0", "grad_1", "grad_2"])
    
    print("Calculated coefs: \n\ta = ", parabola_coefficients['a'], "\n\tb = ", parabola_coefficients['b'], "\n\tc = ", parabola_coefficients['c'])
    grad_0 = 2 * parabola_coefficients['a'] * parabola_arg[0] + parabola_coefficients['b']
    grad_1 = 2 * parabola_coefficients['a'] * parabola_arg[1] + parabola_coefficients['b']
    grad_2 = 2 * parabola_coefficients['a'] * parabola_arg[2] + parabola_coefficients['b']
    
    
    #x = np.array(range(-100, 100, 10))
    #y = parabola_coefficients['a'] * np.power(x, 2) + parabola_coefficients['b'] * x + parabola_coefficients['c']
    #plt.plot(x,y)  
    #plt.grid(True)
    #plt.show() 
    
    return parab_minimum(E_new, parabola_min_arg, grad_0, grad_1, grad_2)
