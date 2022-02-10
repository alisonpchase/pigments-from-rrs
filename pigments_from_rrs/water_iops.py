"""TODO: Add link to Matlab source.
"""
# define input arrays
#WAVEL_WATER = importlib.resources('water_iops', 'data/water_abs.csv')

# Matlab function definition line:
#    `function [a_sw,bb_sw] = get_water_iops(wave,T,S)`

def RInw(lambda_, Tc, S):
    """Refractive index of air is from Ciddor (1996,Applied Optics)
    Refractive index of seawater is from Quan and Fry (1994, Applied Optics)
    """
    n_air = (
        1.0 + (5792105.0 / (238.0185 - 1 / (lambda_ / 1e3) ** 2)
        + 167917.0 / (57.362 - 1 / (lambda_ / 1e3) ** 2)) / 1e8
    )
    n0 = 1.31405
    n1 = 1.779e-4
    n2 = -1.05e-6
    n3 = 1.6e-8
    n4 = -2.02e-6
    n5 = 15.868
    n6 = 0.01155
    n7 = -0.00423
    n8 = -4382
    n9 = 1.1455e6
    nsw = (
        n0 + (n1 + n2 * Tc + n3 * Tc ** 2) * S + n4 * Tc ** 2 
        + (n5 + n6 * S + n7 * Tc) / lambda_ + n8 / lambda_ ** 2 + n9
        / lambda_ ** 3
    )
    nsw = nsw * n_air
    dnswds = (n1 + n2 * Tc + n3 * Tc ** 2 + n6 / lambda_) * n_air
    return nsw, dnswds


def BetaT(Tc, S):
    """
    function IsoComp = BetaT(Tc, S)
    % pure water secant bulk Millero (1980, Deep-sea Research)
    kw = 19652.21+148.4206*Tc-2.327105*Tc.^2+1.360477e-2*Tc.^3-5.155288e-5*Tc.^4;
    Btw_cal = 1./kw;

    % isothermal compressibility from Kell sound measurement in pure water
    % Btw = (50.88630+0.717582*Tc+0.7819867e-3*Tc.^2+31.62214e-6*Tc.^3-0.1323594e-6*Tc.^4+0.634575e-9*Tc.^5)./(1+21.65928e-3*Tc)*1e-6;

    % seawater secant bulk
    a0 = 54.6746-0.603459*Tc+1.09987e-2*Tc.^2-6.167e-5*Tc.^3;
    b0 = 7.944e-2+1.6483e-2*Tc-5.3009e-4*Tc.^2;

    Ks =kw + a0.*S + b0.*S.^1.5;

    % calculate seawater isothermal compressibility from the secant bulk
    IsoComp = 1./Ks*1e-5; % unit is pa
    end
    """
    kw = (
        19652.21 + 148.4206 * Tc - 2.327105 * Tc ** 2 + 1.360477e-2 
        * Tc ** 3 - 5.155288e-5 * Tc ** 4
    )
    Btw_cal = 1 / kw
    a0 = 54.6746 - 0.603459 * Tc + 1.09987e-2 * Tc ** 2-6.167e-5 * Tc ** 3
    b0 = 7.944e-2 + 1.6483e-2 * Tc - 5.3009e-4 * Tc ** 2

    Ks = kw + a0 * S + b0 * S ** 1.5
    IsoComp = 1 / Ks * 1e0 - 5

    return IsoComp


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