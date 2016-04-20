# Multizone Phonon Refinement

Traditionally, phonon measurements have been made at
specific wave vectors in one or two BZs. This works well when
phonons are well separated in energy, but in cases of multiple
overlapping peaks it is often impossible to distinguish a single
broad peak from two or more narrow peaks. TOF inelastic
neutron-scattering instruments capture data over many BZs.
Thus phonon peaks belonging to a multiplet of overlapping
peaks at a specific reduced wave vector will have different
intensities in different BZs, but they will have the same
energies and intrinsic linewidths. To perform the fitting, we
first start with predictions of DFT calculations, which give
approximate phonon frequencies and structure factors. Some
of these are very accurate and agree with the data without any
fitting, whereas others are significantly off (see the left panel
in Fig. 3). By simultaneously fitting the experimental phonon
spectrum at the same reduced wave vector in multiple BZs, we
can identify the experimental positions of individual phonons
quite well if we constrain their widths and energies to be BZ
independent while letting their intensities vary from zone to
zone (see the right panel in Fig. 3).

http://dx.doi.org/10.1103/PhysRevB.89.064310

"van_var.m" is from the Horace package ( http://horace.isis.rl.ac.uk/Main_Page )
"merchop.m" is an brutal hack of MCHOP, also part of Horace.
Both of these functions will go away in future releases, to be replaced with a
resolution-width calculator.

This refinement program uses the "leasqr.m" file from Octave's "optim" package.
You can install optim by typing this on the octave command line :
```
pkg install -forge optim
```
You may have to install dependencies first.
