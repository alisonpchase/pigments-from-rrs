# define input arrays
#WAVEL_WATER = importlib.resources('water_iops', 'data/water_abs.csv')

# Matlab function definition line:
#    `function [a_sw,bb_sw] = get_water_iops(wave,T,S)`

def get_water_iops(wavel, T, S):
    """Description of what this function does.

    Parameters
    ----------
    wavel: (array)
        Wavelenghts (nm) for iops to be returned at.
    T: (float or int)
        Temperature (degC) of the water.
    S: (float or int)
        Salinity (psu) of the water.

    Returns
    -------
    a_sw: (array)
        Absorption of seawater at each wavelength in wavel for the specified temperature and salinity values.
    bb_sw: (array)
        Backscattering of seawater at each wavelenght in wavel for the specified temperature and salinity values.
    """
    pass   # this function doesn't do anything yet, so we're just going to `pass` for now


# Matlab function definition line:
#   `function [betasw124, bsw, beta90sw, theta] = betasw124_ZHH2009(lambda, S, Tc, delta)`

def betasw124_ZHH2009(wavel, S, Tc, delta):
    """
    
    Parameters
    ----------
    wavel: (array)
        Wavelenghts (nm) 
    S: (float or int)
        Salinity of water
    Tc: (float or int)
        Temperature (degC) of the water
    delta: (float)
        Depolarization ratio, default = 0.039

    Returns
    -------
    
    
    """
    pass