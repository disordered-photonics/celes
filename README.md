<p align="center">
<img src="https://disordered-photonics.github.io/celes/readme_logo.svg", height="150">
</p>

## CELES
> CUDA-accelerated Electromagnetic scattering for Large Ensembles of Spheres

CELES (latin for 'fast ship') is a MATLAB/CUDA MEX implementation of the multi-sphere T-matrix method (also known as Generalized Multiparticle Mie method).

The main aim of the software is to rigorously solve electrodynamic problems comprising extremely large numbers of spherical scatterers. As such, it can be applied to study light propagation in macroscopic aggregates of particles in order to derive their bulk transport properties.

![coverimage](doc/images/coverimage.png)

If you use CELES, please cite it as follows:

> Egel A, Pattelli L, Mazzamuto G, Wiersma DS, and Lemmer U.
_CELES: CUDA-accelerated simulation of electromagnetic scattering by large ensembles of spheres_,
Journal of Quantitative Spectroscopy and Radiative Transfer 199C (2017) pp. 103-110. [[link](https://doi.org/10.1016/j.jqsrt.2017.05.010)] [[bibtex](doc/celes.bib)]

### Features
CELES is written in MATLAB in order to provide a user-friendly, fully scriptable interface to configure and run simulations. Its prominent features are

* massively parallel execution on CUDA-capable NVIDIA GPU hardware
* block-diagonal preconditioning for faster convergence of iterative solvers
* lookup-table approach to evaluate spherical Hankel functions
* rich output (power flux, near- and far-field distributions)
* Gaussian beam excitation
* GUI (experimental)

What's new in version 2.0:
* support for polydisperse samples of spheres (thanks to Alan Zhan for his valuable contributions)
* patching of Matlab's GMRES solver for a live convergence progress monitor

### Requirements
In order to run CELES, the following software (in addition to MATLAB) should be installed on your system:
* the [CUDA toolkit](https://developer.nvidia.com/cuda-downloads).
* a [C++ compiler](https://it.mathworks.com/support/compilers.html) which is supported by MATLAB in combination with the given CUDA version.

You can check which CUDA version you need to install by running the command `gpuDevice` in MATLAB and look for `ToolkitVersion` in the output.

CELES has been successfully tested on Linux using the built-in gcc compiler and on Windows with the following software combinations: MATLAB R2016b + CUDA 7.5 + MS Visual Studio 2013 and MATLAB R2017a + CUDA 8 + MS Visual Studio 2013. 
MS Visual Studio 2013 can be downloaded [here](https://www.microsoft.com/en-us/download/details.aspx?id=48138).

In order to fully take advantage of preconditioned iterative solvers we recommend running CELES on a workstation with sufficient RAM (~several 10GB for 10000+ scattering particles).

### Getting started
CELES can be installed via cloning the GitHub repository with
```bash
git clone git@github.com:disordered-photonics/celes.git
```
or by downloading and extracting one of the [releases](https://github.com/disordered-photonics/celes/releases). Please note that the releases do not always represent the most up to date version.

In MATLAB, remember to add CELES to your search path with
```matlab
addpath(genpath('path/to/celes'));
```

As an example input you can execute the `CELES_MAIN` script. Comments in the script explain how the simulation parameters are specified. Alternatively, you can use the `CELES_model_wizard` app, a GUI that helps in the specification of the simulation parameters.

For more information, please refer to the
[documentation](https://disordered-photonics.github.io/celes/).

### Contributing
If you add any improvement or implement new features to the software please consider contributing them following the [GitHub flow](https://guides.github.com/introduction/flow/)

1. Fork the project
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push the branch: `git push origin my-new-feature`
5. Submit a [Pull request](https://github.com/disordered-photonics/celes/pulls)

If you have questions, bug reports or feature requests, please use the [Issues](https://github.com/disordered-photonics/celes/issues) section to report them.

### License
This software is published under the BSD 3-clause license, please read the [LICENCE](LICENCE) file for more information.

### Credits
CELES was initiated by Amos Egel, Lorenzo Pattelli and Giacomo Mazzamuto. In addition, Alan Zhan has contributed code to the project. We thank Yasuhiko Okada for valuable comments and feedback.

CELES uses the following codes from other programmers:
* [polarplot3d](https://it.mathworks.com/matlabcentral/fileexchange/13200-3d-polar-plot/content/polarplot3d.m) from Kenn Gerard
* [wigner3j](https://it.mathworks.com/matlabcentral/fileexchange/20619-wigner3j-symbol) from Kobi Kraus
* Iterative solvers based on the [Templates for the Solution of Linear Systems](http://it.mathworks.com/matlabcentral/fileexchange/2158-templates-for-the-solution-of-linear-systems)
