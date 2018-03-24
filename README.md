# 5dpes

For molecules whose symmetry is described by 3 translational and 2 rotational modes, this code produces the potential energy surface as a function of those internal variables.  The 5D PES resulting from this calculation can include an external influence, such as a nano surface or mean field description of neighbors.  With the 5D PES in hand, additional calculations (such as coupled quantum translations/rotations can be applied).

For use of this 5dpes code, please cite I. Matanovic, J.L. Belof, B. Space, K. Sillar, J. Sauer, J. Eckert and B. Zlatko, "Hydrogen adsorbed in a metal organic framework-5: Coupled translation-rotation eigenstates from quantum five-dimensional calculations", J. Chem. Phys. 137:014701 (2012),  http://doi.org/10.1063/1.4730906

For use of the hydrogen potential used in an example, please cite J.L. Belof, A.C. Stern and B. Space, "An accurate and transferable intermolecular diatomic hydrogen potential for condensed phase simulation", J. Chem. Theory. Comput., 4:1332 (2008), http://doi.org/10.1021/ct800155q

For use of the nanomaterial potential (MOF-5) used in an example, please cite J.L. Belof, A.C. Abraham and B. Space, "A predictive model of hydrogen sorption for metal-organic materials", J. Phys. Chem. C., 113:9316 (2009), http://doi.org/10.1021/jp901988e


## Getting Started

After obtaining the source code, please consult the Makefile to set any specific compiler flags (defaults are gcc with std C library).

The code implements a potential surface that described long-range electrostatics, electronic repulsion and dispersion forces and many-body polarization.


## Installing

Compilation is simple and relies on only standard libraries:

$ make  
gcc -c -O3 -DDEBUG -I. main.c  
gcc -c -O3 -DDEBUG -I. cleanup.c  
gcc -c -O3 -DDEBUG -I. input.c  
gcc -c -O3 -DDEBUG -I. pairs.c  
gcc -c -O3 -DDEBUG -I. pbc.c  
gcc -c -O3 -DDEBUG -I. surface.c  
gcc -c -O3 -DDEBUG -I. energy.c  
gcc -O3 -DDEBUG *.o -o 5dpes  

## Running the examples

Run the 5dpes binary without arguments to obtain the command line input:

$ ./5dpes 
usage: ./5dpes <PDB filename> [<b1_x> <b1_y> <b1_z> <b2_x> <b2_y> <b2_z> <b3_x> <b3_y> <b3_z>] [<xi> <yi> <zi> <xf> <yf> <zf>] [<dx> <dy> <dz> <dtheta> <dphi>]
	[b1,...,b3] : the cartesian basis vectors (in A) of the unit cell, each b vector is a row of the basis matrix
	[xi,...,zf] : the initial and final xyz c.o.m. coordinates for the PES generation
	[dx,...,dphi] : the step size for c.o.m. coordinates and for the spherical polar angles (theta angle of rotation around the Y axis, phi is the angle around the Z axis)
	columnar output is [x,y,z,theta,phi,E(in Kelvin)]
MOF5 example:
	$ ./5dpes MOF5+H2.pdb 25.669 0.0 0.0 0.0 25.669 0.0 0.0 0.0 25.669 -10.0 -10.0 -10.0 10.0 10.0 10.0 0.001 0.001 0.001 0.1 0.1
would map the PES at each point within a subcube of the unit cell spanned by [-10,-10,-10]x[10,10,10] for all angles theta=0-pi,phi=0-2*pi

This main example consists of an orthorhombic unit cell containing the metal-organic framework MOF-5 with a single hydrogen molecule that is then imposed on a grid and rotated.  The resulting output is then:

\# X Y Z THETA PHI ENERGY  
-10.0000000000000000 -10.0000000000000000 -10.0000000000000000 0.0000000000000000 0.0000000000000000 -183.4211687345789414  
-10.0000000000000000 -10.0000000000000000 -10.0000000000000000 0.0000000000000000 0.1000000000000000 -184.2287123170468135  
-10.0000000000000000 -10.0000000000000000 -10.0000000000000000 0.0000000000000000 0.2000000000000000 -184.9762340141523111  

...

with the first 3 numbers the (x,y,z) coordinates, the following 2 numbers the Euler angles for the molecular rotation and finally the energy (in Kelvin).


## Authors

* **Jon Belof** [jbelof@github](https://github.com/jbelof)  

[google scholar](https://scholar.google.com/citations?user=gNrlNbwAAAAJ&hl=en)  
[research gate](https://www.researchgate.net/profile/Jon_Belof)  
[linkedin](http://www.linkedin.com/in/jbelof)  
[web profile](http://jbelof.academia.edu)  


## License

This proejct is licensed under the GNU General Public License v3, please see GPL_license.txt for details.


