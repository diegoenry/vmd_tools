/*
 * Test program for GROMACS topology plugin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* We'll compile the plugin directly into this test */
#define STATIC_PLUGIN
#include "grotopplugin.c"

void print_atom_info(molfile_atom_t *atoms, int natoms, int max_print) {
  printf("\n=== First %d atoms ===\n", max_print);
  printf("%-6s %-8s %-8s %-6s %-8s %-8s %-8s %-8s\n",
         "Index", "Name", "Type", "ResID", "ResName", "SegID", "Charge", "Mass");
  printf("-----------------------------------------------------------------------\n");

  for (int i = 0; i < natoms && i < max_print; i++) {
    printf("%-6d %-8s %-8s %-6d %-8s %-8s %8.3f %8.3f\n",
           i + 1,
           atoms[i].name,
           atoms[i].type,
           atoms[i].resid,
           atoms[i].resname,
           atoms[i].segid,
           atoms[i].charge,
           atoms[i].mass);
  }
}

void print_bond_info(int *from, int *to, int nbonds, int max_print) {
  printf("\n=== First %d bonds ===\n", max_print);
  printf("%-6s %-8s %-8s\n", "Index", "From", "To");
  printf("-------------------------\n");

  for (int i = 0; i < nbonds && i < max_print; i++) {
    printf("%-6d %-8d %-8d\n", i + 1, from[i], to[i]);
  }
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <topology_file.top>\n", argv[0]);
    return 1;
  }

  const char *filename = argv[1];

  /* Initialize the plugin */
  VMDPLUGIN_init();

  printf("=======================================================\n");
  printf("GROMACS Topology Plugin Test\n");
  printf("=======================================================\n");
  printf("Reading file: %s\n", filename);
  printf("=======================================================\n");

  /* Open the file */
  int natoms = 0;
  void *handle = open_grotop_read(filename, "grotop", &natoms);

  if (!handle) {
    fprintf(stderr, "ERROR: Failed to open topology file\n");
    return 1;
  }

  printf("\nSuccessfully opened file\n");
  printf("Total atoms: %d\n", natoms);

  /* Allocate atom array */
  molfile_atom_t *atoms = (molfile_atom_t *)calloc(natoms, sizeof(molfile_atom_t));
  if (!atoms) {
    fprintf(stderr, "ERROR: Failed to allocate memory for atoms\n");
    close_grotop_read(handle);
    return 1;
  }

  /* Read structure */
  int optflags = 0;
  int rc = read_grotop_structure(handle, &optflags, atoms);

  if (rc != MOLFILE_SUCCESS) {
    fprintf(stderr, "ERROR: Failed to read structure\n");
    free(atoms);
    close_grotop_read(handle);
    return 1;
  }

  printf("\nSuccessfully read structure\n");
  printf("Optional flags: 0x%x\n", optflags);
  if (optflags & MOLFILE_CHARGE) printf("  - Has charges\n");
  if (optflags & MOLFILE_MASS) printf("  - Has masses\n");

  /* Print atom information */
  int max_atoms_print = (natoms > 20) ? 20 : natoms;
  print_atom_info(atoms, natoms, max_atoms_print);

  if (natoms > max_atoms_print) {
    printf("... (%d more atoms)\n", natoms - max_atoms_print);
  }

  /* Read bonds */
  int nbonds = 0;
  int *from = NULL, *to = NULL;
  float *bondorder = NULL;
  int *bondtype = NULL;
  int nbondtypes = 0;
  char **bondtypename = NULL;

  rc = read_grotop_bonds(handle, &nbonds, &from, &to, &bondorder,
                         &bondtype, &nbondtypes, &bondtypename);

  if (rc != MOLFILE_SUCCESS) {
    fprintf(stderr, "ERROR: Failed to read bonds\n");
    free(atoms);
    close_grotop_read(handle);
    return 1;
  }

  printf("\nSuccessfully read bonds\n");
  printf("Total bonds: %d\n", nbonds);

  if (nbonds > 0) {
    int max_bonds_print = (nbonds > 20) ? 20 : nbonds;
    print_bond_info(from, to, nbonds, max_bonds_print);

    if (nbonds > max_bonds_print) {
      printf("... (%d more bonds)\n", nbonds - max_bonds_print);
    }
  }

  /* Read angles, dihedrals, impropers */
  int numangles = 0, numdihedrals = 0, numimpropers = 0, numcterms = 0;
  int *angles = NULL, *angletypes = NULL, *dihedrals = NULL, *dihedraltypes = NULL;
  int *impropers = NULL, *impropertypes = NULL, *cterms = NULL;
  int numangletypes = 0, numdihedraltypes = 0, numimpropertypes = 0;
  int ctermcols = 0, ctermrows = 0;
  char **angletypenames = NULL, **dihedraltypenames = NULL, **impropertypenames = NULL;

  /* Check if plugin has read_angles function pointer */
  molfile_plugin_t *p = (molfile_plugin_t *)&plugin;
  printf("\nChecking for read_angles function pointer: %s\n", p->read_angles ? "YES" : "NO");

  if (p->read_angles) {
    rc = p->read_angles(handle, &numangles, &angles, &angletypes,
                        &numangletypes, &angletypenames, &numdihedrals,
                        &dihedrals, &dihedraltypes, &numdihedraltypes,
                        &dihedraltypenames, &numimpropers, &impropers,
                        &impropertypes, &numimpropertypes, &impropertypenames,
                        &numcterms, &cterms, &ctermcols, &ctermrows);

    if (rc == MOLFILE_SUCCESS) {
      printf("\nSuccessfully read angles/dihedrals/impropers\n");
      printf("  Angles: %d\n", numangles);
      printf("  Dihedrals: %d\n", numdihedrals);
      printf("  Impropers: %d\n", numimpropers);
    } else {
      printf("\nFailed to read angles/dihedrals/impropers (rc=%d)\n", rc);
    }
  } else {
    printf("Plugin does not have read_angles function\n");
  }

  /* Summary statistics */
  printf("\n=======================================================\n");
  printf("Summary:\n");
  printf("  Total atoms: %d\n", natoms);
  printf("  Total bonds: %d\n", nbonds);
  if (p->read_angles && rc == MOLFILE_SUCCESS) {
    printf("  Total angles: %d\n", numangles);
    printf("  Total dihedrals: %d\n", numdihedrals);
    printf("  Total impropers: %d\n", numimpropers);
  }

  /* Count unique residues */
  int unique_residues = 0;
  int last_resid = -999999;
  char last_segid[8] = "";
  for (int i = 0; i < natoms; i++) {
    if (atoms[i].resid != last_resid || strcmp(atoms[i].segid, last_segid) != 0) {
      unique_residues++;
      last_resid = atoms[i].resid;
      strncpy(last_segid, atoms[i].segid, sizeof(last_segid) - 1);
    }
  }
  printf("  Unique residues: %d\n", unique_residues);

  /* Count unique segments */
  int unique_segments = 0;
  last_segid[0] = '\0';
  for (int i = 0; i < natoms; i++) {
    if (strcmp(atoms[i].segid, last_segid) != 0) {
      unique_segments++;
      strncpy(last_segid, atoms[i].segid, sizeof(last_segid) - 1);
    }
  }
  printf("  Unique segments: %d\n", unique_segments);

  printf("=======================================================\n");

  /* Clean up */
  free(atoms);
  close_grotop_read(handle);

  printf("\nTest completed successfully!\n");

  return 0;
}
