/*
 * Test program: Read GROMACS topology and write PSF
 *
 * This validates the grotopplugin by using it to read .top files
 * and writing the result to .psf using the standard psfplugin.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* We'll compile both plugins directly into this test */
#define STATIC_PLUGIN
#include "grotopplugin.c"

/* Include PSF plugin */
/* We just redefine the VMDPLUGIN macros to avoid conflicts with grotop */
#undef VMDPLUGIN_init
#undef VMDPLUGIN_register
#undef VMDPLUGIN_fini

#define VMDPLUGIN_init VMDPLUGIN_psf_init
#define VMDPLUGIN_register VMDPLUGIN_psf_register
#define VMDPLUGIN_fini VMDPLUGIN_psf_fini

/* Include the PSF plugin source */
/* Suppress warnings from third-party plugin code */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "/Users/deb0054/gitlab/vmd2prototype/plugins/molfile_plugin/src/psfplugin.c"
#pragma GCC diagnostic pop

/* Now get pointers to the PSF functions from the plugin structure */
/* We need to initialize the PSF plugin to populate the structure */
extern molfile_plugin_t plugin;  /* This is from psfplugin.c */

static void *(*psf_open_write)(const char *, const char *, int);
static int (*psf_write_structure)(void *, int, const molfile_atom_t *);
static int (*psf_write_bonds)(void *, int, int *, int *, float *, int *, int, char **);
static int (*psf_write_angles)(void *, int, const int *, const int *, int, const char **,
                                int, const int *, const int *, int, const char **,
                                int, const int *, const int *, int, const char **,
                                int, const int *, int, int);
static void (*psf_close_write)(void *);

static void init_psf_plugin(void) {
  VMDPLUGIN_psf_init();
  psf_open_write = plugin.open_file_write;
  psf_write_structure = plugin.write_structure;
  psf_write_bonds = plugin.write_bonds;
  psf_write_angles = plugin.write_angles;
  psf_close_write = plugin.close_file_write;
}


int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s <input.top> <output.psf>\n", argv[0]);
    return 1;
  }

  const char *input_file = argv[1];
  const char *output_file = argv[2];

  /* Initialize PSF plugin */
  init_psf_plugin();

  printf("=======================================================\n");
  printf("GROMACS Topology to PSF Converter\n");
  printf("=======================================================\n");
  printf("Input:  %s\n", input_file);
  printf("Output: %s\n", output_file);
  printf("=======================================================\n\n");

  /* Step 1: Read GROMACS topology */
  printf("Step 1: Reading GROMACS topology...\n");

  int natoms = 0;
  void *grotop_handle = open_grotop_read(input_file, "grotop", &natoms);

  if (!grotop_handle) {
    fprintf(stderr, "ERROR: Failed to open topology file\n");
    return 1;
  }

  printf("  - Total atoms: %d\n", natoms);

  /* Allocate atom array */
  molfile_atom_t *atoms = (molfile_atom_t *)calloc(natoms, sizeof(molfile_atom_t));
  if (!atoms) {
    fprintf(stderr, "ERROR: Failed to allocate memory for atoms\n");
    close_grotop_read(grotop_handle);
    return 1;
  }

  /* Read structure */
  int optflags = 0;
  int rc = read_grotop_structure(grotop_handle, &optflags, atoms);

  if (rc != MOLFILE_SUCCESS) {
    fprintf(stderr, "ERROR: Failed to read structure\n");
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  printf("  - Read structure successfully\n");
  printf("  - Optional flags: 0x%x\n", optflags);

  /* Read bonds */
  int nbonds = 0;
  int *from = NULL, *to = NULL;
  float *bondorder = NULL;
  int *bondtype = NULL;
  int nbondtypes = 0;
  char **bondtypename = NULL;

  rc = read_grotop_bonds(grotop_handle, &nbonds, &from, &to, &bondorder,
                         &bondtype, &nbondtypes, &bondtypename);

  if (rc != MOLFILE_SUCCESS) {
    fprintf(stderr, "ERROR: Failed to read bonds\n");
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  printf("  - Total bonds: %d\n", nbonds);

  /* Read angles, dihedrals, impropers */
  int numangles = 0, numdihedrals = 0, numimpropers = 0, numcterms = 0;
  int *angles = NULL, *angletypes = NULL, *dihedrals = NULL, *dihedraltypes = NULL;
  int *impropers = NULL, *impropertypes = NULL, *cterms = NULL;
  int numangletypes = 0, numdihedraltypes = 0, numimpropertypes = 0;
  int ctermcols = 0, ctermrows = 0;
  char **angletypenames = NULL, **dihedraltypenames = NULL, **impropertypenames = NULL;

  rc = read_grotop_angles(grotop_handle, &numangles, &angles, &angletypes,
                          &numangletypes, &angletypenames, &numdihedrals,
                          &dihedrals, &dihedraltypes, &numdihedraltypes,
                          &dihedraltypenames, &numimpropers, &impropers,
                          &impropertypes, &numimpropertypes, &impropertypenames,
                          &numcterms, &cterms, &ctermcols, &ctermrows);

  if (rc != MOLFILE_SUCCESS) {
    fprintf(stderr, "ERROR: Failed to read angles/dihedrals/impropers\n");
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  printf("  - Total angles: %d\n", numangles);
  printf("  - Total dihedrals: %d\n", numdihedrals);
  printf("  - Total impropers: %d\n\n", numimpropers);

  /* Step 2: Write PSF file */
  printf("Step 2: Writing PSF file...\n");

  void *psf_handle = psf_open_write(output_file, "psf", natoms);

  if (!psf_handle) {
    fprintf(stderr, "ERROR: Failed to open PSF file for writing\n");
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  /* IMPORTANT: Write bonds BEFORE structure! */
  /* The PSF plugin's write_bonds() saves bond data for later use by write_structure() */
  if (nbonds > 0) {
    rc = psf_write_bonds(psf_handle, nbonds, from, to, bondorder,
                         bondtype, nbondtypes, bondtypename);

    if (rc != MOLFILE_SUCCESS) {
      fprintf(stderr, "ERROR: Failed to save PSF bonds\n");
      psf_close_write(psf_handle);
      free(atoms);
      close_grotop_read(grotop_handle);
      return 1;
    }

    printf("  - Saved bonds to PSF structure\n");
  }

  /* Write structure (this actually writes the complete PSF file including bonds) */
  rc = psf_write_structure(psf_handle, optflags, atoms);

  if (rc != MOLFILE_SUCCESS) {
    fprintf(stderr, "ERROR: Failed to write PSF structure\n");
    psf_close_write(psf_handle);
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  printf("  - Wrote structure successfully\n");

  /* Write angles, dihedrals, and impropers */
  if (numangles > 0 || numdihedrals > 0 || numimpropers > 0) {
    rc = psf_write_angles(psf_handle,
                          numangles, angles, angletypes, numangletypes,
                          (const char **)angletypenames,
                          numdihedrals, dihedrals, dihedraltypes, numdihedraltypes,
                          (const char **)dihedraltypenames,
                          numimpropers, impropers, impropertypes, numimpropertypes,
                          (const char **)impropertypenames,
                          numcterms, cterms, ctermcols, ctermrows);

    if (rc != MOLFILE_SUCCESS) {
      fprintf(stderr, "ERROR: Failed to write angles/dihedrals/impropers to PSF\n");
      psf_close_write(psf_handle);
      free(atoms);
      close_grotop_read(grotop_handle);
      return 1;
    }

    printf("  - Wrote angles/dihedrals/impropers successfully\n");
  }

  printf("  - Wrote complete PSF file successfully\n");

  /* Step 3: Clean up */
  printf("\nStep 3: Cleaning up...\n");

  psf_close_write(psf_handle);
  close_grotop_read(grotop_handle);
  free(atoms);

  printf("\n=======================================================\n");
  printf("SUCCESS: PSF file written to %s\n", output_file);
  printf("=======================================================\n");
  printf("\nSummary:\n");
  printf("  Atoms:      %d\n", natoms);
  printf("  Bonds:      %d\n", nbonds);
  printf("  Angles:     %d\n", numangles);
  printf("  Dihedrals:  %d\n", numdihedrals);
  printf("  Impropers:  %d\n", numimpropers);
  printf("\n");

  return 0;
}
