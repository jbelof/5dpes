/* 

@2009, Jonathan Belof

*/

#include <5dpes.h>

/* returns the total potential energy for the system */
double energy(system_t *system) {

	double potential_energy, rd_energy, coulombic_energy;

	/* get the pairwise terms necessary for the energy calculation */
	pairs(system);

	/* get the repulsion/dispersion potential */
	rd_energy = lj(system);

	/* get the Ewald summed electrostatic potential */
	coulombic_energy = coulombic(system);

	/* sum the total potential energy */
	potential_energy = rd_energy + coulombic_energy;

	return(potential_energy);


}

/* Lennard-Jones repulsion/dispersion */
double lj(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	double sigma_over_r, term12, term6;
	double potential;

	potential = 0;
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				/* make sure we're not excluded or beyond the cutoff */
				if(!((pair_ptr->rimg > system->pbc->cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {

					/* the LJ potential */
					sigma_over_r = pair_ptr->sigma/pair_ptr->rimg;
					term6 = pow(sigma_over_r, 6.0);
					term12 = pow(term6, 2.0);
					potential += 4.0*pair_ptr->epsilon*(term12 - term6);


				}


			} /* pair */
		} /* atom */
	} /* molecule */

	return(potential);

}


/* total ES energy term */
double coulombic(system_t *system) {

	double real, reciprocal, self;
	double potential;

	real = coulombic_real(system);
	reciprocal = coulombic_reciprocal(system);
	self = coulombic_self_point(system) + coulombic_self_intra(system);

	potential = real + reciprocal - self;

	return(potential);

}


/* fourier space sum */
double coulombic_reciprocal(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	int q, p;
	int kmax;
	double alpha;
	double l[3], k[3], k_squared, norm;
	double gaussian, position_product;
	double SF_real, SF_imaginary;			/* structure factor */
	double potential;

	alpha = system->ewald_alpha;
	kmax = system->ewald_kmax;

	potential = 0;
	/* perform the fourier sum over the reciprocal lattice for each particle */
	for(l[0] = 0; l[0] <= kmax; l[0]++) {
		for(l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
			for(l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {

				/* compare the norm */
				for(p = 0, norm = 0; p < 3; p++)
					norm += l[p]*l[p];

				/* get the reciprocal lattice vectors */
				for(p = 0, k_squared = 0; p < 3; p++) {
					for(q = 0, k[p] = 0; q < 3; q++)
						k[p] += system->pbc->reciprocal_basis[p][q]*2.0*M_PI*((double)l[q]);
					k_squared += k[p]*k[p];
				}

				if((norm <= kmax*kmax) && (k_squared > 0.0)) {

					gaussian = exp(-k_squared/(4.0*alpha*alpha))/k_squared;

					SF_real = 0; SF_imaginary = 0;
					for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
						for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

							if(!atom_ptr->frozen) {

								/* the inner product of the position vector and the k vector */
								/* also take the inner product of the dipole and k vector */
								for(p = 0, position_product = 0; p < 3; p++)
									position_product += k[p]*atom_ptr->pos[p];

								SF_real += atom_ptr->charge*cos(position_product);
								SF_imaginary += atom_ptr->charge*sin(position_product);

							} /* !frozen */

						} /* atom */
					} /* molecule */
					potential += gaussian*(SF_real*SF_real + SF_imaginary*SF_imaginary);

				} /* end if norm */

			} /* end for n */
		} /* end for m */
	} /* end for l */

	potential *= 4.0*M_PI/system->pbc->volume;

	return(potential);

}

/* real space sum */
double coulombic_real(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	double alpha, r, erfc_term, gaussian_term;
	double potential;

	alpha = system->ewald_alpha;

	potential = 0;
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(!pair_ptr->frozen) {

					r = pair_ptr->rimg;

					if(!((r > system->pbc->cutoff) || pair_ptr->es_excluded)) {	/* unit cell part */

						erfc_term = erfc(alpha*r);
						gaussian_term = exp(-alpha*alpha*r*r);

						potential += atom_ptr->charge*pair_ptr->charge*erfc_term/r;

					}
				} /* !frozen */

			} /* pair */
		} /* atom */
	} /* molecule */

	return(potential);

}

/* point self-interaction term */
double coulombic_self_point(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	double alpha, potential;

	alpha = system->ewald_alpha;

	potential = 0;
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			if(!atom_ptr->frozen) {

				potential += alpha*atom_ptr->charge*atom_ptr->charge/sqrt(M_PI);

			}

		}
	}

	return(potential);

}


double coulombic_self_intra(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	double alpha;
	double potential;

	alpha = system->ewald_alpha;

	potential = 0;
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(!pair_ptr->frozen) {
					if(pair_ptr->es_excluded) { 	/* include neighbors excluded from the first image */

						potential += atom_ptr->charge*pair_ptr->charge*erf(alpha*pair_ptr->r)/pair_ptr->r;

					}
				} /* !frozen */

			} /* pair */
		} /* atom */
	} /* molecule */

	return(potential);

}



