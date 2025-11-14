# VMD Molfile Plugin Development Guide

## Converting Python Implementation to VMD Plugin

This guide explains how to port the GROMACS to PSF Python code to a VMD molfile plugin in C.

## VMD Molfile Plugin Architecture

VMD molfile plugins are shared libraries that implement a standard API for reading/writing molecular file formats. The plugin registers callbacks for various operations.

### Required Plugin Interface

```c
// Plugin registration structure
typedef struct {
  vmdplugin_HEAD
  const char *filename_extension;  // e.g., "top,itp"
  void *(*open_file_read)(const char *filepath, const char *filetype, int *natoms);
  int (*read_structure)(void *handle, int *optflags, molfile_atom_t *atoms);
  int (*read_bonds)(void *handle, int *nbonds, int **from, int **to, float **bondorder,
                    int **bondtype, int *nbondtypes, char ***bondtypename);
  int (*read_angles)(void *handle, int *numangles, int **angles, int **angletypes,
                     int *numangletypes, char ***angletypenames, float **angleforces);
  int (*read_dihedrals)(void *handle, int *numdihedrals, int **dihedrals, int **dihedraltypes,
                        int *numdihedraltypes, char ***dihedraltypenames, float **dihedralforces);
  void (*close_file_read)(void *handle);
} molfile_plugin_t;
```

## Data Structure Mapping

### Python to C Struct Conversion

#### Atom Structure
```c
// Python: Atom dataclass
// C equivalent:
typedef struct {
    int atom_id;
    char atom_type[8];
    int residue_number;
    char residue_name[8];
    char atom_name[8];
    int charge_group;
    float charge;
    float mass;
    char segment_id[8];
    int original_id;
    char molecule_type[32];
} gromacs_atom_t;
```

#### MoleculeType Structure
```c
// Python: MoleculeType dataclass
// C equivalent:
typedef struct {
    char name[64];
    int nrexcl;
    int natoms;
    gromacs_atom_t *atoms;
    int nbonds;
    int *bonds;  // Pairs of atom indices
    int nangles;
    int *angles;  // Triples of atom indices
    int ndihedrals;
    int *dihedrals;  // Quads of atom indices
    int nimpropers;
    int *impropers;  // Quads of atom indices
} molecule_type_t;
```

#### Plugin Handle Structure
```c
// Python: Topology dataclass
// C plugin handle:
typedef struct {
    FILE *fp;
    char system_name[256];

    // Molecule type definitions
    int nmoltypes;
    molecule_type_t *moltypes;

    // System instantiation
    int nmolecules;
    char **molecule_names;
    int *molecule_counts;

    // Complete system
    int natoms;
    gromacs_atom_t *atoms;
    int nbonds;
    int *bonds;
    int nangles;
    int *angles;
    int ndihedrals;
    int *dihedrals;
    int nimpropers;
    int *impropers;
} gromacs_handle_t;
```

## Key Functions to Implement

### 1. File Opening
```c
void *open_gromacs_read(const char *filepath, const char *filetype, int *natoms) {
    gromacs_handle_t *handle = calloc(1, sizeof(gromacs_handle_t));

    // Parse .top file (analogous to parse_top_file())
    parse_top_file(filepath, handle);

    // Build system (analogous to build_system())
    build_system(handle);

    *natoms = handle->natoms;
    return handle;
}
```

### 2. Reading Structure
```c
int read_gromacs_structure(void *v, int *optflags, molfile_atom_t *atoms) {
    gromacs_handle_t *handle = (gromacs_handle_t *)v;

    *optflags = MOLFILE_MASS | MOLFILE_CHARGE | MOLFILE_ATOMICNUMBER;

    for (int i = 0; i < handle->natoms; i++) {
        gromacs_atom_t *ga = &handle->atoms[i];
        molfile_atom_t *a = &atoms[i];

        strncpy(a->name, ga->atom_name, sizeof(a->name));
        strncpy(a->type, ga->atom_type, sizeof(a->type));
        strncpy(a->resname, ga->residue_name, sizeof(a->resname));
        strncpy(a->segid, ga->segment_id, sizeof(a->segid));
        strncpy(a->chain, ga->segment_id, sizeof(a->chain));

        a->resid = ga->residue_number;
        a->charge = ga->charge;
        a->mass = ga->mass;
        a->atomicnumber = guess_atomic_number(ga->atom_name);
    }

    return MOLFILE_SUCCESS;
}
```

### 3. Reading Bonds
```c
int read_gromacs_bonds(void *v, int *nbonds, int **from, int **to,
                       float **bondorder, int **bondtype,
                       int *nbondtypes, char ***bondtypename) {
    gromacs_handle_t *handle = (gromacs_handle_t *)v;

    *nbonds = handle->nbonds;
    *from = malloc(*nbonds * sizeof(int));
    *to = malloc(*nbonds * sizeof(int));

    for (int i = 0; i < *nbonds; i++) {
        (*from)[i] = handle->bonds[2*i] - 1;      // Convert to 0-based
        (*to)[i] = handle->bonds[2*i + 1] - 1;    // Convert to 0-based
    }

    *bondorder = NULL;
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;

    return MOLFILE_SUCCESS;
}
```

### 4. Reading Angles
```c
int read_gromacs_angles(void *v, int *numangles, int **angles,
                        int **angletypes, int *numangletypes,
                        char ***angletypenames, float **angleforces) {
    gromacs_handle_t *handle = (gromacs_handle_t *)v;

    *numangles = handle->nangles;
    *angles = malloc(*numangles * 3 * sizeof(int));

    for (int i = 0; i < *numangles; i++) {
        (*angles)[3*i]     = handle->angles[3*i] - 1;      // Convert to 0-based
        (*angles)[3*i + 1] = handle->angles[3*i + 1] - 1;
        (*angles)[3*i + 2] = handle->angles[3*i + 2] - 1;
    }

    *angletypes = NULL;
    *numangletypes = 0;
    *angletypenames = NULL;
    *angleforces = NULL;

    return MOLFILE_SUCCESS;
}
```

### 5. File Closing
```c
void close_gromacs_read(void *v) {
    gromacs_handle_t *handle = (gromacs_handle_t *)v;

    // Free all allocated memory
    for (int i = 0; i < handle->nmoltypes; i++) {
        free(handle->moltypes[i].atoms);
        free(handle->moltypes[i].bonds);
        free(handle->moltypes[i].angles);
        free(handle->moltypes[i].dihedrals);
        free(handle->moltypes[i].impropers);
    }
    free(handle->moltypes);

    free(handle->atoms);
    free(handle->bonds);
    free(handle->angles);
    free(handle->dihedrals);
    free(handle->impropers);

    for (int i = 0; i < handle->nmolecules; i++) {
        free(handle->molecule_names[i]);
    }
    free(handle->molecule_names);
    free(handle->molecule_counts);

    free(handle);
}
```

## Parsing Functions

The Python parsing functions translate directly to C:

### Line Processing
```c
// Python: strip_comments()
char *strip_comments(char *line) {
    char *comment = strchr(line, ';');
    if (comment) *comment = '\0';

    // Trim whitespace
    while (isspace(*line)) line++;
    char *end = line + strlen(line) - 1;
    while (end > line && isspace(*end)) *end-- = '\0';

    return line;
}

// Python: is_section_header()
int is_section_header(const char *line, char *section_name) {
    if (line[0] == '[' && line[strlen(line)-1] == ']') {
        strncpy(section_name, line + 1, strlen(line) - 2);
        section_name[strlen(line) - 2] = '\0';
        return 1;
    }
    return 0;
}
```

### Atom Parsing
```c
// Python: parse_atom_line()
int parse_atom_line(const char *line, gromacs_atom_t *atom) {
    int nr, resnr, cgnr;
    char type[16], residue[16], atomname[16];
    float charge, mass;

    int n = sscanf(line, "%d %s %d %s %s %d %f %f",
                   &nr, type, &resnr, residue, atomname, &cgnr, &charge, &mass);

    if (n < 7) return 0;

    atom->original_id = nr;
    strncpy(atom->atom_type, type, sizeof(atom->atom_type));
    atom->residue_number = resnr;
    strncpy(atom->residue_name, residue, sizeof(atom->residue_name));
    strncpy(atom->atom_name, atomname, sizeof(atom->atom_name));
    atom->charge_group = cgnr;
    atom->charge = charge;
    atom->mass = (n >= 8) ? mass : 0.0;

    return 1;
}
```

## File Handling with #include

```c
void parse_top_file_recursive(const char *filepath, gromacs_handle_t *handle,
                              const char *base_dir) {
    FILE *fp = fopen(filepath, "r");
    if (!fp) return;

    char line[1024];
    char section[64] = "";

    while (fgets(line, sizeof(line), fp)) {
        char *clean_line = strip_comments(line);

        if (strlen(clean_line) == 0) continue;

        // Handle #include
        if (strncmp(clean_line, "#include", 8) == 0) {
            char include_file[256];
            if (sscanf(clean_line, "#include \"%[^\"]\"", include_file) == 1 ||
                sscanf(clean_line, "#include '%[^']'", include_file) == 1) {

                char full_path[512];
                snprintf(full_path, sizeof(full_path), "%s/%s", base_dir, include_file);

                if (strstr(include_file, ".itp")) {
                    parse_itp_file_recursive(full_path, handle, base_dir);
                }
            }
            continue;
        }

        // Handle section headers
        if (is_section_header(clean_line, section)) {
            continue;
        }

        // Parse section content
        if (strcmp(section, "system") == 0) {
            strncpy(handle->system_name, clean_line, sizeof(handle->system_name));
        }
        else if (strcmp(section, "molecules") == 0) {
            // Parse molecule instances
            // ... (similar to Python implementation)
        }
    }

    fclose(fp);
}
```

## Build System

```bash
# Makefile for VMD plugin
CC = gcc
CFLAGS = -fPIC -Wall -O2 -I$(VMD_PLUGIN_DIR)/include
LDFLAGS = -shared

all: gmxtopplugin.so

gmxtopplugin.so: gmxtopplugin.o
    $(CC) $(LDFLAGS) -o $@ $<

gmxtopplugin.o: gmxtopplugin.c
    $(CC) $(CFLAGS) -c $<

clean:
    rm -f *.o *.so

install: gmxtopplugin.so
    cp gmxtopplugin.so $(VMD_PLUGIN_DIR)/plugins/MACOSXX86_64/molfile/
```

## Testing the Plugin

```tcl
# VMD script to test plugin
mol new example.top type gmxtop
mol modstyle 0 0 CPK
```

## Key Differences Between Python and C

1. **Memory Management**: C requires explicit malloc/free
2. **String Handling**: Use strncpy, snprintf instead of Python string operations
3. **Dynamic Arrays**: Need to track sizes and reallocate as needed
4. **Error Handling**: Return codes instead of exceptions
5. **File I/O**: Use FILE* and fgets() instead of Python's file objects

## Development Workflow

1. Start with Python implementation to verify logic
2. Create C structs matching Python dataclasses
3. Port parsing functions one at a time
4. Test each function independently
5. Integrate into VMD plugin framework
6. Test with VMD

## Resources

- VMD Plugin Developer Guide: https://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/
- VMD Molfile Plugin Source: `$VMDDIR/plugins/molfile_plugin/src/`
- Example plugins: `psfplugin.c`, `gromacsplugin.c` (for .gro coordinates)

## Notes

- The Python version handles topology only; coordinates come from .gro/.pdb files
- VMD plugins typically separate structure (topology) and trajectory (coordinates) reading
- Consider implementing both .top (topology) and enhanced .gro (with topology hints) support
- GROMACS force field parameters are not needed for visualization, only connectivity
