import numpy as np


def calculate_spectra(use_mock_data):
    """Calculate and return spectra

    Parameters
    ----------
    use_mock_data: bool
        Whether or not to use mock data.

    Returns
    -------
    dict
    """

    if use_mock_data:
        wavelengths = np.arange(400, 750)  # units: nm

        spectra = {
            "wavelengths": wavelengths,
            "Chla": ((np.sin((wavelengths/10)) + 1) / 60),
            "Chlb": (np.cos(wavelengths/15) + 1) / 60,
            "Chlc": (np.cos(wavelengths/20) + 1) / 60,
            "PSC": (np.cos(wavelengths/25) + 1) / 60,
            "PPC": (np.cos(wavelengths/30) + 1) / 60,
            "PE PEB": (np.cos(wavelengths/35) + 1) / 60,
            "PE PUB": (np.cos(wavelengths/40) + 1) / 60,
        }
        return spectra
    else:
        # Here's where real spectra would be returned
        raise ValueError("`use_mock_data` must be set to `True`")
