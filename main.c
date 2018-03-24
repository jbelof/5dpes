/* 

@2009, Jonathan Belof

*/


#include <5dpes.h>


void usage(char *progname) {

	fprintf(stderr, "usage: %s <PDB filename> [<b1_x> <b1_y> <b1_z> <b2_x> <b2_y> <b2_z> <b3_x> <b3_y> <b3_z>] [<xi> <yi> <zi> <xf> <yf> <zf>] [<dx> <dy> <dz> <dtheta> <dphi>]\n", progname);
	fprintf(stderr, "\t[b1,...,b3] : the cartesian basis vectors (in A) of the unit cell, each b vector is a row of the basis matrix\n");
	fprintf(stderr, "\t[xi,...,zf] : the initial and final xyz c.o.m. coordinates for the PES generation\n");
	fprintf(stderr, "\t[dx,...,dphi] : the step size for c.o.m. coordinates and for the spherical polar angles (theta angle of rotation around the Y axis, phi is the angle around the Z axis)\n");
	fprintf(stderr, "\tcolumnar output is [x,y,z,theta,phi,E(in Kelvin)]\n");
	fprintf(stderr, "MOF5 example:\n");
	fprintf(stderr, "\t$ %s MOF5+H2.pdb 25.669 0.0 0.0 0.0 25.669 0.0 0.0 0.0 25.669 -10.0 -10.0 -10.0 10.0 10.0 10.0 0.001 0.001 0.001 0.1 0.1\n", progname);
	fprintf(stderr, "would map the PES at each point within a subcube of the unit cell spanned by [-10,-10,-10]x[10,10,10] for all angles theta=0-pi,phi=0-2*pi\n");
	exit(1);

}

int main(int argc, char **argv) {

	system_t *system;
	double b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z;
	double xi, yi, zi, xf, yf, zf;
	double dx, dy, dz, dtheta, dphi;

	/* check args */
	if(argc != 22) usage(argv[0]);

	if(!argv[1]) {
		fprintf(stderr, "MAIN: invalid PDB file specified");
		exit(1);
	}

	/* get the basis vectors for the unit cell */
	b1x = atof(argv[2]); b1y = atof(argv[3]); b1z = atof(argv[4]);
	b2x = atof(argv[5]); b2y = atof(argv[6]); b2z = atof(argv[7]);
	b3x = atof(argv[8]); b3y = atof(argv[9]); b3z = atof(argv[10]);

	/* read the geometry and potential parameters and setup the data structures */
	system = setup_system(argv[1], b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z);
	if(!system) {
		fprintf(stderr, "MAIN: could not initialize the data structures\n");
		exit(1);
	}

	/* get the surface boundaries and increments */
	xi = atof(argv[11]); yi = atof(argv[12]); zi = atof(argv[13]);
	xf = atof(argv[14]); yf = atof(argv[15]); zf = atof(argv[16]);
	dx = atof(argv[17]); dy = atof(argv[18]); dz = atof(argv[19]); dtheta = atof(argv[20]); dphi = atof(argv[21]);

	/* calculate the PES at the requested points and print the energy (in K) to stdout */
	surface_scan(system, xi, yi, zi, xf, yf, zf, dx, dy, dz, dtheta, dphi);


	/* cleanup */
	cleanup(system);
	exit(0);

}

