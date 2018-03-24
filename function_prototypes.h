
/* main.c */
void usage(char *);
void free_pairs(molecule_t *);
void free_atoms(molecule_t *);
void free_molecule(molecule_t *);
void free_molecules(molecule_t *);

/* cleanup.c */
void cleanup(system_t *);

/* input.c */
system_t *setup_system(char *, double, double, double, double, double, double, double, double, double);
molecule_t *read_molecules(system_t *);

/* pairs.c */
void setup_pairs(molecule_t *);
void pairs(system_t *);
void update_com(molecule_t *);

/* pbc.c */
void pbc(pbc_t *);
void pbc_reciprocal(pbc_t *);
double pbc_volume(pbc_t *);
double pbc_cutoff(pbc_t *);

/* surface.c */
void surface_scan(system_t *, double, double, double, double, double, double, double, double, double, double, double);
double surface_point(system_t *, double, double, double, double, double);
void translate(molecule_t *, double, double, double);
void rotate(molecule_t *, double, double);
molecule_t *copy_molecule(molecule_t *);

/* energy.c */
double energy(system_t *);
double lj(system_t *);
double coulombic(system_t *);
double coulombic_reciprocal(system_t *);
double coulombic_real(system_t *);
double coulombic_self_point(system_t *);
double coulombic_self_intra(system_t *);

#ifdef DEBUG
/* debugging */
void test_pairs(molecule_t *);
void test_list(molecule_t *);
void test_molecule(molecule_t *);
int write_molecules(system_t *, char *);
void write_states(FILE *, molecule_t *);
#endif /* DEBUG */


