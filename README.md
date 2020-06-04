
# Python wrapper for the Bayesian Distance Calculator (BD_wrapper)

## About
The ``BD_wrapper`` is a Python wrapper for the Bayesian Distance Calculator (BDC) tool presented in [Reid et al. 2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...823...77R/abstract) and [Reid et al. 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...885..131R/abstract). The BDC uses a Bayesian approach to establish distances to sources within the Milky Way given Galactic coordinates (GLON, GLAT) and a line of sight velocity (VLSR). The BDC is a data product from the [The Bar And Spiral Structure Legacy (BeSSeL) Survey](http://bessel.vlbi-astrometry.org/) and can be accessed from its webpage. Currently, the ``BD_wrapper`` supports both v1 and v2.4 of the BDC.

The ``BD_wrapper`` provides a Python implementation for bulk distance assignments for a large number of input sources. Moreover, it contains two new priors that help to solve for the kinematic distance ambiguity (KDA). The first new prior considers literature solutions for the KDA by automatically matching the source coordinates with literature sources. The second new prior requires additional information about the velocity dispersion of the source and uses a size-linewidth relationship to solve for the KDA.

For a detailed description about the new priors and results of tests performed with the ``BD_wrapper`` on decomposition results from the Galactic Ring Survey ([Riener et al. 2020a](https://ui.adsabs.harvard.edu/abs/2020A%26A...633A..14R/abstract)) please see:

Riener et al. 2020b (in prep.)

All credit for the (BDC) and its original Fortran implementation is due to [Reid et al. 2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...823...77R/abstract) and [Reid et al. 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...885..131R/abstract) and should be acknowledged as such.

For tips on how to get started with the ``BD_wrapper`` see the section [Getting started](#gettingstarted) further below. See also the notebook [Parameter_settings.ipynb](tutorials/Parameter_settings.ipynb) for an overview about the parameter settings.

### Version

The currently recommended version of the ``BD_wrapper`` is v0.1. See the [BD_wrapper Changelog](CHANGES.md) for an overview of the major changes and improvements introduced by newer versions currently in development.

New updates to the code are first tested and developed in the ``dev`` branch. Users cloning the ``dev`` branch should beware that these versions are not guaranteed to be stable.

## Installation

### Dependencies

The ``BD_wrapper`` already includes the necessary files to run the BDC tool. However, you need to be able to execute Fortran (.f) files on your operating system. Moreover, you will need the following packages to run the ``BD_wrapper``. We list the version of each package which we know to be compatible with the ``BD_wrapper``.

* [python 3.5](https://www.python.org/)
* [astropy (v3.0.4)](http://www.astropy.org/)
* [matplotlib (v2.2.2)](http://matplotlib.org/)
* [numpy (v1.14.2)](http://www.numpy.org/)
* [tqdm (v4.19.4)](https://tqdm.github.io/)

If you do not already have Python 3.5, you can install the [Anaconda Scientific Python distribution](https://store.continuum.io/cshop/anaconda/), which comes pre-loaded with numpy.

### Download the BD_wrapper

Download the BD_wrapper using git `$ git clone https://github.com/mriener/BD_wrapper.git`


### Installing Dependencies on Linux

Install pip for easy installation of python packages:

```bash
sudo apt-get install python-pip
```

Then install the required python packages:

```bash
sudo pip install astropy matplotlib numpy tqdm
```

### Installing Dependencies on OSX

Install pip for easy installation of python packages:

```bash
sudo easy_install pip
```

Then install the required python packages:

```bash
sudo pip install astropy matplotlib numpy tqdm
```

<a id="gettingstarted"></a>
## Getting started

The notebook [Parameter_settings.ipynb](tutorials/Parameter_settings.ipynb) gives an overview about parameter settings for the ``BD_wrapper``.

The [Tutorial-batch_distance_estimation.ipynb](tutorials/Tutorial-batch_distance_estimation.ipynb) notebook guides users through the process of using the `BD_wrapper` to estimate distances for a large table of input sources.

The [Tutorial-plot_distance_pdf.ipynb](tutorials/Tutorial-plot_distance_pdf.ipynb) notebook shows how to plot the distance probability density results obtained by the BDC and how to retain temporary files that can be important for debugging and obtaining diagnostics of the distance calculation.

The ``BD_wrapper`` uses literature distance results to inform the P_far prior that helps resolve the kinematic distance ambiguity (KDA). The ``BD_wrapper`` currently contains twelve catalogues (called KDA info tables) that mostly cover regions in the first Galactic quadrant. In the [Example-KDA_info_table.ipynb](tutorials/Example-KDA_info_table.ipynb) we show how to easily create a new KDA info table for a catalogue that is not yet included in the `KDA_info` directory.

## Citing the BD_wrapper

If you make use of this package in a publication, please consider the following acknowledgements:

```
@ARTICLE{2019ApJ...885..131R,
       author = {{Reid}, M.~J. and {Menten}, K.~M. and {Brunthaler}, A. and
         {Zheng}, X.~W. and {Dame}, T.~M. and {Xu}, Y. and {Li}, J. and
         {Sakai}, N. and {Wu}, Y. and {Immer}, K. and {Zhang}, B. and
         {Sanna}, A. and {Moscadelli}, L. and {Rygl}, K.~L.~J. and
         {Bartkiewicz}, A. and {Hu}, B. and {Quiroga-Nu{\~n}ez}, L.~H. and
         {van Langevelde}, H.~J.},
        title = "{Trigonometric Parallaxes of High-mass Star-forming Regions: Our View of the Milky Way}",
      journal = {\apj},
     keywords = {Milky Way, Milky Way dynamics, Milky Way rotation, Trigonometric parallax, Star formation, Gravitational wave sources, Astrophysics - Astrophysics of Galaxies},
         year = 2019,
        month = nov,
       volume = {885},
       number = {2},
          eid = {131},
        pages = {131},
          doi = {10.3847/1538-4357/ab4a11},
archivePrefix = {arXiv},
       eprint = {1910.03357},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019ApJ...885..131R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{2016ApJ...823...77R,
       author = {{Reid}, M.~J. and {Dame}, T.~M. and {Menten}, K.~M. and {Brunthaler}, A.},
        title = "{A Parallax-based Distance Estimator for Spiral Arm Sources}",
      journal = {\apj},
     keywords = {Galaxy: structure, parallaxes, stars: formation, Astrophysics - Astrophysics of Galaxies},
         year = 2016,
        month = jun,
       volume = {823},
       number = {2},
          eid = {77},
        pages = {77},
          doi = {10.3847/0004-637X/823/2/77},
archivePrefix = {arXiv},
       eprint = {1604.02433},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2016ApJ...823...77R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
Citation courtesy of [ADS](https://ui.adsabs.harvard.edu/#).

Please also consider acknowledgements to the required dependencies in your work.


## Feedback

If you should find that the ``BD_wrapper`` does not perform as intended for your dataset or if you should come across bugs or have suggestions for improvement, please get into contact with us or open a new Issue or Pull request.

## Contributing to GaussPy+

To contribute to the ``BD_wrapper``, see [Contributing to the BD_wrapper](CONTRIBUTING.md)
