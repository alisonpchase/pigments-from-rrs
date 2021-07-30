import numpy as np
import matplotlib.pyplot as plt


def plot_spectra(
    spectra,
    title,
    figsize=(10, 8),
    style='seaborn-whitegrid',
    fontsize=14,
    labelpad=6.0,
):
    """Plot pigment curves

    Parameters
    ----------
    spectra: dict
        Keys are label names, values are 1-dimensional `numpy.ndarray`s.
        The dictionary must contain a key `"wavelengths"` correseponding to a
        `numpy.ndarray` of x-axis values for the plot.
    title: str
        A title for the plot
    figsize: tuple
        A size for the figure
    style: str
        One of the styles listed at
        https://matplotlib.org/stable/tutorials/introductory/customizing.html
    fontsize: int
        A fontsize for the labels
    labelpad: float
        Spacing between axis and labels

    Notes
    -----
    Plots a figure via the `matplotlib.pyplot` interface

    Returns
    -------
    None
    """
    plt.style.use(style)
    plt.rcParams["figure.figsize"] = figsize
    plt.rcParams["axes.labelpad"] = labelpad

    for key, value in spectra.items():
        if key != "wavelengths":
            plt.plot(spectra["wavelengths"], value, label=key)

    yticks = np.arange(start=0, stop=0.08, step=0.01)
    plt.yticks(yticks, fontsize=fontsize-2)
    plt.xticks(fontsize=fontsize-2)
    plt.ylim(top=0.08)

    plt.ylabel("$a * (m^2 mg^{-1})$", fontsize=fontsize)
    plt.xlabel("wavelength (nm)", fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(loc='upper right')

    plt.margins(0)
    plt.show()
