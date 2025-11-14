/*
 * Test program: Read GROMACS topology and coordinates, write JS format
 *
 * This validates the grotopplugin by using it to read .top files,
 * reads coordinates from .gro file using gromacsplugin,
 * and writes the complete result to .js using the jsplugin.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* We'll compile both plugins directly into this test */
#define STATIC_PLUGIN
#include "grotopplugin.c"

/* Include gromacsplugin for reading coordinates */
/* Declare the wrapper functions from gromacs_wrapper.cpp */
extern void *open_gro_read_wrapper(const char *filename, const char *filetype, int *natoms);
extern int read_gro_structure_wrapper(void *v, int *optflags, molfile_atom_t *atoms);
extern int read_gro_timestep_wrapper(void *v, int natoms, molfile_timestep_t *ts);
extern void close_gro_read_wrapper(void *v);

/* Include JS plugin */
/* We just redefine the VMDPLUGIN macros to avoid conflicts with grotop */
#undef VMDPLUGIN_init
#undef VMDPLUGIN_register
#undef VMDPLUGIN_fini

#define VMDPLUGIN_init VMDPLUGIN_js_init
#define VMDPLUGIN_register VMDPLUGIN_js_register
#define VMDPLUGIN_fini VMDPLUGIN_js_fini

/* Include the JS plugin source */
/* Suppress warnings from third-party plugin code */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wformat"
#include "/Users/deb0054/gitlab/vmd2prototype/plugins/molfile_plugin/src/jsplugin.c"
#pragma GCC diagnostic pop

/* Now get pointers to the JS functions from the plugin structure */
/* We need to initialize the JS plugin to populate the structure */
extern molfile_plugin_t plugin;  /* This is from jsplugin.c */

static void *(*js_open_write)(const char *, const char *, int);
static int (*js_write_structure)(void *, int, const molfile_atom_t *);
static int (*js_write_bonds)(void *, int, int *, int *, float *, int *, int, char **);
static int (*js_write_angles)(void *, int, const int *, const int *, int, const char **,
                               int, const int *, const int *, int, const char **,
                               int, const int *, const int *, int, const char **,
                               int, const int *, int, int);
static int (*js_write_timestep)(void *, const molfile_timestep_t *);
static void (*js_close_write)(void *);

static void init_js_plugin(void) {
  VMDPLUGIN_js_init();
  js_open_write = plugin.open_file_write;
  js_write_structure = plugin.write_structure;
  js_write_bonds = plugin.write_bonds;
  js_write_angles = plugin.write_angles;
  js_write_timestep = plugin.write_timestep;
  js_close_write = plugin.close_file_write;
}


int main(int argc, char *argv[]) {
  if (argc < 4) {
    fprintf(stderr, "Usage: %s <input.top> <input.gro> <output.js>\n", argv[0]);
    return 1;
  }

  const char *top_file = argv[1];
  const char *gro_file = argv[2];
  const char *output_file = argv[3];

  /* Initialize JS plugin */
  init_js_plugin();

  printf("=======================================================\n");
  printf("GROMACS Topology and Coordinates to JS Converter\n");
  printf("=======================================================\n");
  printf("Topology:    %s\n", top_file);
  printf("Coordinates: %s\n", gro_file);
  printf("Output:      %s\n", output_file);
  printf("=======================================================\n\n");

  /* Step 1: Read GROMACS topology */
  printf("Step 1: Reading GROMACS topology...\n");

  int natoms = 0;
  void *grotop_handle = open_grotop_read(top_file, "grotop", &natoms);

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

  /* Step 2: Read coordinates from GRO file */
  printf("Step 2: Reading coordinates from GRO file...\n");

  int gro_natoms = 0;
  void *gro_handle = open_gro_read_wrapper(gro_file, "gro", &gro_natoms);

  if (!gro_handle) {
    fprintf(stderr, "ERROR: Failed to open GRO file\n");
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  if (gro_natoms != natoms) {
    fprintf(stderr, "ERROR: Atom count mismatch between topology (%d) and GRO (%d)\n",
            natoms, gro_natoms);
    close_gro_read_wrapper(gro_handle);
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  printf("  - GRO file contains %d atoms (matches topology)\n", gro_natoms);

  /* Allocate separate atom array for GRO file (we'll extract coords from it) */
  molfile_atom_t *gro_atoms = (molfile_atom_t *)calloc(gro_natoms, sizeof(molfile_atom_t));
  if (!gro_atoms) {
    fprintf(stderr, "ERROR: Failed to allocate memory for GRO atoms\n");
    close_gro_read_wrapper(gro_handle);
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  /* Read structure from GRO file (this includes coordinates) */
  int gro_optflags = 0;
  printf("  - Reading structure and coordinates from GRO file...\n");
  rc = read_gro_structure_wrapper(gro_handle, &gro_optflags, gro_atoms);

  if (rc != MOLFILE_SUCCESS) {
    fprintf(stderr, "ERROR: Failed to read structure from GRO file (return code %d)\n", rc);
    free(gro_atoms);
    close_gro_read_wrapper(gro_handle);
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  printf("  - Read structure successfully (optflags: 0x%x)\n", gro_optflags);
  printf("  - First atom: %s %s\n", gro_atoms[0].resname, gro_atoms[0].name);

  /* Now read timestep to get coordinates */
  molfile_timestep_t *ts = (molfile_timestep_t *)calloc(1, sizeof(molfile_timestep_t));
  if (!ts) {
    fprintf(stderr, "ERROR: Failed to allocate memory for timestep\n");
    free(gro_atoms);
    close_gro_read_wrapper(gro_handle);
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  /* Allocate coordinate array */
  ts->coords = (float *)calloc(3 * natoms, sizeof(float));
  if (!ts->coords) {
    fprintf(stderr, "ERROR: Failed to allocate memory for coordinates\n");
    free(ts);
    free(gro_atoms);
    close_gro_read_wrapper(gro_handle);
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  /* Read coordinates using timestep reader */
  printf("  - Reading coordinates...\n");
  rc = read_gro_timestep_wrapper(gro_handle, natoms, ts);

  if (rc != MOLFILE_SUCCESS) {
    fprintf(stderr, "ERROR: Failed to read coordinates (return code %d)\n", rc);
    free(ts->coords);
    free(ts);
    free(gro_atoms);
    close_gro_read_wrapper(gro_handle);
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  printf("  - Read coordinates successfully\n");
  printf("  - First atom coordinates: (%.3f, %.3f, %.3f)\n",
         ts->coords[0], ts->coords[1], ts->coords[2]);

  free(gro_atoms);
  close_gro_read_wrapper(gro_handle);
  printf("\n");

  /* Step 3: Write JS file */
  printf("Step 3: Writing JS file...\n");

  void *js_handle = js_open_write(output_file, "js", natoms);

  if (!js_handle) {
    fprintf(stderr, "ERROR: Failed to open JS file for writing\n");
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  /* IMPORTANT: Write bonds BEFORE structure! */
  /* The molfile API requires write_bonds() to be called before write_structure() */
  if (nbonds > 0) {
    rc = js_write_bonds(js_handle, nbonds, from, to, bondorder,
                        bondtype, nbondtypes, bondtypename);

    if (rc != MOLFILE_SUCCESS) {
      fprintf(stderr, "ERROR: Failed to write JS bonds\n");
      js_close_write(js_handle);
      free(ts->coords);
      free(ts);
      free(atoms);
      close_grotop_read(grotop_handle);
      return 1;
    }

    printf("  - Saved bonds to JS structure\n");
  }

  /* Write structure (this includes bond data saved above) */
  rc = js_write_structure(js_handle, optflags, atoms);

  if (rc != MOLFILE_SUCCESS) {
    fprintf(stderr, "ERROR: Failed to write JS structure\n");
    js_close_write(js_handle);
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  printf("  - Wrote structure successfully\n");

  /* Write angles, dihedrals, and impropers */
  if (numangles > 0 || numdihedrals > 0 || numimpropers > 0) {
    rc = js_write_angles(js_handle,
                         numangles, angles, angletypes, numangletypes,
                         (const char **)angletypenames,
                         numdihedrals, dihedrals, dihedraltypes, numdihedraltypes,
                         (const char **)dihedraltypenames,
                         numimpropers, impropers, impropertypes, numimpropertypes,
                         (const char **)impropertypenames,
                         numcterms, cterms, ctermcols, ctermrows);

    if (rc != MOLFILE_SUCCESS) {
      fprintf(stderr, "ERROR: Failed to write angles/dihedrals/impropers to JS\n");
      js_close_write(js_handle);
      free(ts->coords);
      free(ts);
      free(atoms);
      close_grotop_read(grotop_handle);
      return 1;
    }

    printf("  - Wrote angles/dihedrals/impropers successfully\n");
  }

  /* Write coordinates */
  rc = js_write_timestep(js_handle, ts);

  if (rc != MOLFILE_SUCCESS) {
    fprintf(stderr, "ERROR: Failed to write coordinates\n");
    js_close_write(js_handle);
    free(ts->coords);
    free(ts);
    free(atoms);
    close_grotop_read(grotop_handle);
    return 1;
  }

  printf("  - Wrote coordinates successfully\n");

  /* Step 4: Clean up */
  printf("\nStep 4: Cleaning up...\n");

  js_close_write(js_handle);
  close_grotop_read(grotop_handle);
  free(ts->coords);
  free(ts);
  free(atoms);

  printf("\n=======================================================\n");
  printf("SUCCESS: JS file written to %s\n", output_file);
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
