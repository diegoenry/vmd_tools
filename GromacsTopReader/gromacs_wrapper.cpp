/*
 * C++ wrapper for gromacsplugin
 * Exposes C-compatible interface for reading GRO files
 */

/* Prevent duplicate symbol errors by disabling VMDPLUGIN functions */
/* The main program already defines these from grotopplugin and jsplugin */
/* We must define these before any includes */
#define VMDPLUGIN vmdplugin_gromacs
#define STATIC_PLUGIN

/* Suppress warnings from third-party plugin code */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "/Users/deb0054/gitlab/vmd2prototype/plugins/molfile_plugin/src/gromacsplugin.C"
#pragma GCC diagnostic pop

extern "C" {

/* Re-export the functions we need with C linkage */
void *open_gro_read_wrapper(const char *filename, const char *filetype, int *natoms) {
  return open_gro_read(filename, filetype, natoms);
}

int read_gro_structure_wrapper(void *v, int *optflags, molfile_atom_t *atoms) {
  return read_gro_structure(v, optflags, atoms);
}

int read_gro_timestep_wrapper(void *v, int natoms, molfile_timestep_t *ts) {
  return read_gro_timestep(v, natoms, ts);
}

void close_gro_read_wrapper(void *v) {
  close_gro_read(v);
}

}  /* extern "C" */
