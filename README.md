# DynaProg
[![View DynaProg on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/84260-dynaprog)

Solve multi-stage deterministic decision problems.

## Purpose
DynaProg is a MATLAB toolbox to solve a finite horizon multi-stage deterministic decision problem, which is a problem where a decision must be made at each stage for a system that evolves through a finite number of stages, minimizing the total cost incurred.

## Installing DynaProg
There are two ways to use DynaProg:
- Download the toolbox from MATLAB File Exchange,
- Download the source code.

### Install the toolbox
The most straightforward way to install DynaProg is to directly install it from MATLAB's Add-on explorer or from the [File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/84260-dynaprog).
Doing this also installs the documentation in MATLAB's Help Browser.

### Download the source code
Alternitavely, you can simply clone the git repository. Then, add the `src` folder to your MATLAB's search path. This method is recommended if you want to inspect/modify the source code.

Doing this will not install the documentation in MATLAB's Help Browser. You can still visualize it with any web browser: browse to the `src/html` folder and open `index.html`.

## Documentation
There are three ways to access the documentation:
- **From the MATLAB Help Browser.** Open the Help Browser. You can find _DynaProg Toolbox_ in the Contents menu (on the left), under _Supplemental Software_.
This option is only available if you installed DynaProg from MATLAB's Add-on explorer or from the [File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/84260-dynaprog).
- **Online.** You can access the online documentation at [this link](https://fmiretti.github.io/DynaProg/).
- **From the html files.** If you downloaded DynaProg as a .zip folder, you can find the docs in the _html_ folder. You should start by opening _index.html_ in your browser.

## Licensing
DynaProg is available under the [MIT license](LICENSE.md).

Some of the examples packaged with the toolbox make use of data derived from the [ADVISOR(R) Software](https://sourceforge.net/p/adv-vehicle-sim/code/HEAD/tree/).

If you use this software for research purposes, please cite the accompanying paper:
Federico Miretti, Daniela Misul, Ezio Spessa, DynaProg: Deterministic Dynamic Programming solver for finite horizon multi-stage decision problems, SoftwareX, Volume 14, 2021, 100690, ISSN 2352-7110, https://doi.org/10.1016/j.softx.2021.100690.


## Contact
For all questions or suggestions, feel free to contact me at federico.miretti@polito.it.
