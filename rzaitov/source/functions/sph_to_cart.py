import numpy as np

# Proposed API. There is no k parameter.
# Here we follow python's name conventions for functions (all lower case with underscore)
def sph_to_cart(trd, plg):
    cn = np.cos(trd) * np.cos(plg)
    ce = np.sin(trd) * np.cos(plg)
    cd = np.sin(plg)
    return cn, ce, cd