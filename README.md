# MScThesisAdRC_AnnexI

The scripts are part of the Master's Thesis "Microelectronic design of a pseudo voltage clamp for a 0.18um CMOS based mesoscale neural interface for intracellular in-vitro recording and stimulation", available at http://repository.tudelft.nl

by

Aitor del Rivero Cortázar

at

Section Bioelectronics, Department of Microelectronics, EEMCS, Delft University of Technology, the Netherlands.

Bio Engineering Lab, D-BSSE, ETH Zürich, Switzerland.

The scripts run on SLiCAP (Symbolic Linear Circuit Analysis Program) for Python, by Anton Montagne.
SLiCAP is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
More information available at: http://www.analog-electronics.eu/slicap/slicap.html

------


PLEASE, BEFORE RUNNING ANY OF THE SCRIPTS INSTALL THE REQUIRED SOFTWARE:


1. Install Python 3 (https://www.microsoft.com/en-us/p/python-37/9nj46sx7x90p or equivalent)
2. Install the following packages for Pyhton 3 (python -m pip install ............):
    - docutils==0.16
    - numpy>=1.18.1
    - sympy>=1.5.1
    - scipy>=1.4.1
    - ply>=3.11
    - matplotlib>=3.1.3
    - in-place==0.4.0
    - sphinx-rtd-theme==0.5.0
3. Install MaximaCAS (http://maxima.sourceforge.net/download.html)
4. Install LTSpice: (https://www.analog.com/en/design-center/design-tools-and-calculators.html)
5. Install SliCAP:
    git clone https://github.com/aitordelrivero/SLiCAP_python.git
    python setup.py install --user


------

PLEASE NOTE THAT

1) The scripts might not run for different SLiCAP versions. The scripts have been tested with the version at https://github.com/aitordelrivero/SLiCAP_python/tree/e3f4d5fddc5466398a813fc3eebde7475d7c4795.

2) The scripts can run with default EKV parameters. Process parameters are available in the internal Wiki only for authorized users.

------

Last update: 15/08/2021
Aitor del Rivero Cortázar. Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
