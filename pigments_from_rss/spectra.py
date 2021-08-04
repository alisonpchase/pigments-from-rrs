import numpy as np


def _read_raw_data_from_source(filepath=None, use_mock_data=False):
    """Read data from source file or location

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

        raw = {
            "wavelengths": wavelengths,
            "Chla": (np.sin((wavelengths/10)) + 1),
            "Chlb": (np.cos(wavelengths/15) + 1),
            "Chlc": (np.cos(wavelengths/20) + 1),
            "PSC": (np.cos(wavelengths/25) + 1),
            "PPC": (np.cos(wavelengths/30) + 1),
            "PE PEB": (np.cos(wavelengths/35) + 1),
            "PE PUB": (np.cos(wavelengths/40) + 1),
        }
    elif not use_mock_data:
        # Here's where real spectra would be read from a file
        # and turned into a Python dictionary
        raise ValueError(f"Reading data from {filepath} not yet implemented")
    else:
        raise ValueError(f"got `use_mock_data={use_mock_data}`; must be bool")

    return raw


def _divide_raw_curves(use_mock_data, divisor):
    """Calculate and return spectra

    Parameters
    ----------
    use_mock_data: bool
        Whether or not to use mock data.

    divisor: int
        An integer by which to divide dict values

    Returns
    -------
    dict
    """
    raw = _read_raw_data_from_source(use_mock_data=use_mock_data)

    return {k: (v / divisor) for k, v in raw.items()}


def calculate_spectra(use_mock_data, divisor):
    """Calculate and return spectra

    Parameters
    ----------
    use_mock_data: bool
        Whether or not to use mock data.

    divisor: int
        An integer by which to divide dict values

    Returns
    -------
    dict
    """

    return _divide_raw_curves(use_mock_data=use_mock_data, divisor=divisor)
