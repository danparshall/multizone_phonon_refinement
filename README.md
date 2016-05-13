# Multizone Phonon Refinement

Traditionally, phonon measurements have been made at specific wave vectors in one or two BZs. This works well when phonons are well separated in energy, but in cases of multiple overlapping peaks it is often impossible to distinguish a single broad peak from two or more narrow peaks. TOF inelastic neutron-scattering instruments capture data over many BZs.  Thus phonon peaks belonging to a multiplet of overlapping peaks at a specific reduced wave vector will have different intensities in different BZs, but they will have the same energies and intrinsic linewidths. To perform the fitting, we first start with predictions of DFT calculations, which give approximate phonon frequencies and structure factors. Some of these are very accurate and agree with the data without any fitting, whereas others are significantly off (see the left panel in Fig. 3). By simultaneously fitting the experimental phonon spectrum at the same reduced wave vector in multiple BZs, we can identify the experimental positions of individual phonons quite well if we constrain their widths and energies to be BZ independent while letting their intensities vary from zone to zone (see the right panel in Fig. 3).

http://dx.doi.org/10.1103/PhysRevB.89.064310

"van_var.m" is from the Horace package ( http://horace.isis.rl.ac.uk/Main_Page )
"merchop.m" is an brutal hack of MCHOP, also part of Horace.
Both of these functions will go away in future releases, to be replaced with a
resolution-width calculator.

This refinement program uses the "lsqnonlin.m" command from the "optim" package (and hence is compatible with the same command from Matlab).
You can install optim by typing this on the Octave command line :
```
pkg install -forge struct
pkg install -forge optim
```
You may have to install dependencies first.


Octave has a well-known bug when using matrices containing more than 2**31 elements.  On many systems, one will have to recompile Octave and enable 64-bit indexing.  Fortunately, the heavy lifting for this has already been done.  See the following links:

http://calaba.tumblr.com/post/107087607479/octave-64
https://github.com/calaba/octave-3.8.2-enable-64-ubuntu-14.04

On Linux Mint 17.2, in addition to the scripts described above, I also had to run the following:
'''
sudo apt-get install libblas-dev liblapack-dev libqhull-dev libglpk-dev libqrupdate-dev libsuitesparse-dev libarpack2 libarpack2-dev
sudo apt-get install freeglut3-dev mesa-common-dev libfltk1.3-dev
export JAVA_HOME=/usr/lib/jvm/default-java
'''

