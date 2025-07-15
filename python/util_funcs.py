import numpy as np

def all_max_ind(a, eps = 1E-9):
    """
    returns indices of all maxes within a tolerances
    """
    if len(a) == 0:
        return []
    all_ = [0]
    max_ = a[0]
    for i in range(1, len(a)):
        if a[i] > max_:
            all_ = [i]
            max_ = a[i]
        elif np.abs(a[i] - max_) < eps:
            all_.append(i)
    return all_
    