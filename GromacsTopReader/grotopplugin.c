/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: grotopplugin.c,v $
 *      $Author: Generated with Claude Code $
 *      $Date: 2025/11/07 $
 *
 ***************************************************************************/

/*
 * GROMACS Topology (.top, .itp) file reader plugin
 *
 * Reads GROMACS topology files including:
 * - Molecule type definitions
 * - Atom types, charges, masses
 * - Bond connectivity
 * - Include directives
 * - Multiple molecule instantiation
 *
 * Forcefield Support:
 * - Automatically loads atom masses from forcefield files (MARTINI, CHARMM, AMBER, etc.)
 * - Supports both MARTINI v3 format (name mass charge ...)
 *   and GROMACS standard format (name bond_type atomic_num mass ...)
 * - These forcefield .itp files could be shipped with VMD for convenience
 */

#include "molfile_plugin.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#define GROTOP_RECORD_LENGTH 512
#define MAX_INCLUDES 100
#define MAX_MOLTYPES 500      /* Increased for large force fields like Martini */
#define MAX_MOLECULES 1000
#define MAX_ATOMTYPES 1000    /* Increased for large force fields */
#define MAX_DEFINES 100       /* Maximum number of #define symbols */
#define MAX_IFDEF_DEPTH 20    /* Maximum nesting depth for #ifdef */

/* Forward declarations */
typedef struct moltype_t moltype_t;
typedef struct atom_data_t atom_data_t;
typedef struct bond_data_t bond_data_t;
typedef struct atomtype_t atomtype_t;

/* Atom data within a molecule type */
struct atom_data_t {
  int id;                    /* Atom ID within molecule (1-based) */
  char atom_type[16];        /* Atom type */
  int resnr;                 /* Residue number */
  char residue[8];           /* Residue name */
  char atom_name[16];        /* Atom name */
  int cgnr;                  /* Charge group number */
  float charge;              /* Partial charge */
  float mass;                /* Atomic mass */
};

/* Bond data within a molecule type */
struct bond_data_t {
  int ai, aj;                /* Atom indices (1-based within molecule) */
};

/* Angle data within a molecule type */
typedef struct {
  int ai, aj, ak;            /* Atom indices (1-based within molecule) */
} angle_data_t;

/* Dihedral data within a molecule type */
typedef struct {
  int ai, aj, ak, al;        /* Atom indices (1-based within molecule) */
  int funct;                 /* Function type (for identifying impropers) */
} dihedral_data_t;

/* Atom type definition */
struct atomtype_t {
  char name[16];             /* Atom type name */
  float mass;                /* Atomic mass */
};

/* Molecule type definition */
struct moltype_t {
  char name[32];             /* Molecule type name */
  int nrexcl;                /* Number of exclusions */
  atom_data_t *atoms;        /* Array of atoms */
  int natoms;                /* Number of atoms */
  int atoms_allocated;       /* Allocated size */
  bond_data_t *bonds;        /* Array of bonds */
  int nbonds;                /* Number of bonds */
  int bonds_allocated;       /* Allocated size */
  angle_data_t *angles;      /* Array of angles */
  int nangles;               /* Number of angles */
  int angles_allocated;      /* Allocated size */
  dihedral_data_t *dihedrals; /* Array of dihedrals */
  int ndihedrals;            /* Number of dihedrals */
  int dihedrals_allocated;   /* Allocated size */
};

/* Main topology data structure */
typedef struct {
  FILE *fp;
  char filepath[512];

  /* Molecule type definitions */
  moltype_t *moltypes[MAX_MOLTYPES];
  int num_moltypes;

  /* Atom type definitions */
  atomtype_t atomtypes[MAX_ATOMTYPES];
  int num_atomtypes;

  /* Molecules section (instances) */
  char molnames[MAX_MOLECULES][32];
  int molcounts[MAX_MOLECULES];
  int num_molecules;

  /* Instantiated system */
  int total_atoms;
  int total_bonds;
  int total_angles;
  int total_dihedrals;
  int total_impropers;

  /* For returning to VMD */
  int *bond_from;
  int *bond_to;
  int *angles;        /* 3 ints per angle: i, j, k */
  int *dihedrals;     /* 4 ints per dihedral: i, j, k, l */
  int *impropers;     /* 4 ints per improper: i, j, k, l */

  /* Preprocessor state */
  char defines[MAX_DEFINES][64];  /* Defined symbols */
  int num_defines;

} grotop_data;


/*
 * Helper Functions
 */

/* Strip comments (everything after semicolon) and trim whitespace */
static char* strip_comments(char *line) {
  char *comment = strchr(line, ';');
  if (comment) *comment = '\0';

  /* Trim trailing whitespace */
  int len = strlen(line);
  while (len > 0 && isspace(line[len-1])) {
    line[--len] = '\0';
  }

  /* Trim leading whitespace */
  char *start = line;
  while (*start && isspace(*start)) start++;

  if (start != line) {
    memmove(line, start, strlen(start) + 1);
  }

  return line;
}

/* Check if line is a section header [section_name] */
static int is_section_header(const char *line, char *section_name) {
  const char *p = line;

  /* Skip leading whitespace */
  while (*p && isspace(*p)) p++;

  if (*p != '[') return 0;
  p++;

  /* Extract section name */
  const char *start = p;
  while (*p && *p != ']' && *p != '\0') p++;

  if (*p != ']') return 0;

  int len = p - start;
  if (len == 0 || len >= 64) return 0;

  strncpy(section_name, start, len);
  section_name[len] = '\0';

  /* Trim spaces from section name */
  char *s = section_name;
  while (*s && isspace(*s)) s++;
  if (s != section_name) {
    memmove(section_name, s, strlen(s) + 1);
  }

  len = strlen(section_name);
  while (len > 0 && isspace(section_name[len-1])) {
    section_name[--len] = '\0';
  }

  return 1;
}

/* Check if a symbol is defined */
static int is_defined(grotop_data *data, const char *symbol) {
  for (int i = 0; i < data->num_defines; i++) {
    if (strcmp(data->defines[i], symbol) == 0) {
      return 1;
    }
  }
  return 0;
}

/* Add a defined symbol */
static void add_define(grotop_data *data, const char *symbol) {
  if (data->num_defines >= MAX_DEFINES) {
    fprintf(stderr, "grotopplugin) WARNING: Maximum number of #define symbols (%d) exceeded\n", MAX_DEFINES);
    return;
  }

  /* Check if already defined */
  if (is_defined(data, symbol)) {
    return;
  }

  strncpy(data->defines[data->num_defines], symbol, 63);
  data->defines[data->num_defines][63] = '\0';
  data->num_defines++;
  printf("grotopplugin) Defined symbol: %s\n", symbol);
}

/* Parse #define directive */
static int parse_define(const char *line, grotop_data *data) {
  const char *p = line;

  /* Skip whitespace */
  while (*p && isspace(*p)) p++;

  /* Check for #define */
  if (strncmp(p, "#define", 7) != 0) return 0;
  p += 7;

  /* Skip whitespace */
  while (*p && isspace(*p)) p++;

  /* Extract symbol name */
  const char *start = p;
  while (*p && !isspace(*p)) p++;

  int len = p - start;
  if (len == 0 || len >= 64) return 0;

  char symbol[64];
  strncpy(symbol, start, len);
  symbol[len] = '\0';

  add_define(data, symbol);
  return 1;
}

/* Parse #include directive */
static int parse_include(const char *line, const char *base_path, char *include_path) {
  const char *p = line;

  /* Skip whitespace */
  while (*p && isspace(*p)) p++;

  /* Check for #include */
  if (strncmp(p, "#include", 8) != 0) return 0;
  p += 8;

  /* Skip whitespace */
  while (*p && isspace(*p)) p++;

  /* Find quoted filename */
  char quote = *p;
  if (quote != '"' && quote != '\'') return 0;
  p++;

  const char *start = p;
  while (*p && *p != quote) p++;
  if (*p != quote) return 0;

  int len = p - start;
  if (len == 0 || len >= 256) return 0;

  char filename[256];
  strncpy(filename, start, len);
  filename[len] = '\0';

  /* Build full path: base_path/filename */
  /* Extract directory from base_path */
  const char *last_slash = strrchr(base_path, '/');
  if (last_slash) {
    int dir_len = last_slash - base_path + 1;
    strncpy(include_path, base_path, dir_len);
    include_path[dir_len] = '\0';
    strcat(include_path, filename);
  } else {
    strcpy(include_path, filename);
  }

  return 1;
}

/* Parse #ifdef or #ifndef directive and return symbol name */
static int parse_ifdef(const char *line, char *symbol, int *is_ifndef) {
  const char *p = line;

  /* Skip whitespace */
  while (*p && isspace(*p)) p++;

  /* Check for #ifdef or #ifndef */
  *is_ifndef = 0;
  if (strncmp(p, "#ifdef", 6) == 0) {
    p += 6;
  } else if (strncmp(p, "#ifndef", 7) == 0) {
    p += 7;
    *is_ifndef = 1;
  } else {
    return 0;
  }

  /* Skip whitespace */
  while (*p && isspace(*p)) p++;

  /* Extract symbol name */
  const char *start = p;
  while (*p && !isspace(*p)) p++;

  int len = p - start;
  if (len == 0 || len >= 64) return 0;

  strncpy(symbol, start, len);
  symbol[len] = '\0';

  return 1;
}

/* Check if line is #else directive */
static int is_else_directive(const char *line) {
  const char *p = line;

  /* Skip whitespace */
  while (*p && isspace(*p)) p++;

  return strncmp(p, "#else", 5) == 0;
}

/* Check if line is #endif directive */
static int is_endif_directive(const char *line) {
  const char *p = line;

  /* Skip whitespace */
  while (*p && isspace(*p)) p++;

  return strncmp(p, "#endif", 6) == 0;
}

/* Check if line is any preprocessor directive */
static int is_preprocessor_directive(const char *line) {
  const char *p = line;

  /* Skip whitespace */
  while (*p && isspace(*p)) p++;

  return (*p == '#');
}

/* Find molecule type by name */
static moltype_t* find_moltype(grotop_data *data, const char *name) {
  for (int i = 0; i < data->num_moltypes; i++) {
    if (strcmp(data->moltypes[i]->name, name) == 0) {
      return data->moltypes[i];
    }
  }
  return NULL;
}

/* Find atom type mass */
static float find_atomtype_mass(grotop_data *data, const char *type) {
  for (int i = 0; i < data->num_atomtypes; i++) {
    if (strcmp(data->atomtypes[i].name, type) == 0) {
      return data->atomtypes[i].mass;
    }
  }
  return 0.0f;  /* Default if not found */
}

/* Create new molecule type */
static moltype_t* create_moltype(void) {
  moltype_t *mt = (moltype_t *)calloc(1, sizeof(moltype_t));
  if (!mt) return NULL;

  mt->atoms_allocated = 100;
  mt->atoms = (atom_data_t *)calloc(mt->atoms_allocated, sizeof(atom_data_t));

  mt->bonds_allocated = 100;
  mt->bonds = (bond_data_t *)calloc(mt->bonds_allocated, sizeof(bond_data_t));

  mt->angles_allocated = 100;
  mt->angles = (angle_data_t *)calloc(mt->angles_allocated, sizeof(angle_data_t));

  mt->dihedrals_allocated = 100;
  mt->dihedrals = (dihedral_data_t *)calloc(mt->dihedrals_allocated, sizeof(dihedral_data_t));

  if (!mt->atoms || !mt->bonds || !mt->angles || !mt->dihedrals) {
    if (mt->atoms) free(mt->atoms);
    if (mt->bonds) free(mt->bonds);
    if (mt->angles) free(mt->angles);
    if (mt->dihedrals) free(mt->dihedrals);
    free(mt);
    return NULL;
  }

  return mt;
}

/* Free molecule type */
static void free_moltype(moltype_t *mt) {
  if (!mt) return;
  if (mt->atoms) free(mt->atoms);
  if (mt->bonds) free(mt->bonds);
  if (mt->angles) free(mt->angles);
  if (mt->dihedrals) free(mt->dihedrals);
  free(mt);
}


/*
 * Section Parsing Functions
 */

/* Parse [ atomtypes ] section */
static int parse_atomtypes_section(FILE *fp, grotop_data *data) {
  char line[GROTOP_RECORD_LENGTH];
  int atomtypes_start = data->num_atomtypes;

  while (1) {
    long pos_before = ftell(fp);  /* Save position BEFORE reading this line */
    if (!fgets(line, sizeof(line), fp)) break;

    /* Check for preprocessor directives BEFORE stripping comments */
    if (is_preprocessor_directive(line)) {
      fseek(fp, pos_before, SEEK_SET);
      printf("grotopplugin)   Loaded %d atomtypes\n", data->num_atomtypes - atomtypes_start);
      return 1;
    }

    strip_comments(line);
    if (strlen(line) == 0) continue;

    /* Check for new section */
    char section[64];
    if (is_section_header(line, section)) {
      printf("grotopplugin)   Loaded %d atomtypes\n", data->num_atomtypes - atomtypes_start);
      /* Seek back and return so main parser can handle new section */
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    /* Parse atomtype line - support multiple formats:
     * MARTINI v3:          name mass charge ptype sigma epsilon
     * GROMACS (full):      name bond_type atomic_num mass charge ptype sigma epsilon
     * GROMACS (simple):    name mass
     */
    char name[16];
    float mass;
    int parsed = 0;

    /* Try MARTINI format: name mass ... (mass is 2nd field) */
    if (sscanf(line, "%15s %f", name, &mass) == 2) {
      parsed = 1;
    }
    /* Try GROMACS full format: name bond_type atomic_num mass ... (mass is 4th field) */
    else if (sscanf(line, "%15s %*s %*s %f", name, &mass) == 2) {
      parsed = 1;
    }

    if (parsed && data->num_atomtypes < MAX_ATOMTYPES) {
      strncpy(data->atomtypes[data->num_atomtypes].name, name, 15);
      data->atomtypes[data->num_atomtypes].name[15] = '\0';
      data->atomtypes[data->num_atomtypes].mass = mass;
      data->num_atomtypes++;
    }
  }

  printf("grotopplugin)   Loaded %d atomtypes\n", data->num_atomtypes - atomtypes_start);
  return 1;
}

/* Parse [ moleculetype ] section header */
static int parse_moleculetype_header(FILE *fp, moltype_t *mt) {
  char line[GROTOP_RECORD_LENGTH];

  /* Next non-empty line should be: name nrexcl */
  while (fgets(line, sizeof(line), fp)) {
    strip_comments(line);
    if (strlen(line) == 0) continue;

    /* Check if it's another section */
    char section[64];
    if (is_section_header(line, section)) {
      return 0;  /* No molecule type data found */
    }

    /* Parse: name nrexcl */
    char name[32];
    int nrexcl = 3;
    if (sscanf(line, "%31s %d", name, &nrexcl) >= 1) {
      strncpy(mt->name, name, 31);
      mt->name[31] = '\0';
      mt->nrexcl = nrexcl;
      return 1;
    }
  }

  return 0;
}

/* Forward declare for recursion */
static int process_section(FILE *fp, const char *section, grotop_data *data, moltype_t **current_mt);

/* Parse [ atoms ] section within a moleculetype */
static int parse_atoms_section(FILE *fp, moltype_t *mt, grotop_data *data, moltype_t **current_mt) {
  char line[GROTOP_RECORD_LENGTH];

  while (1) {
    long pos_before = ftell(fp);  /* Save position BEFORE reading this line */
    if (!fgets(line, sizeof(line), fp)) break;

    /* Check for preprocessor directives BEFORE stripping comments */
    /* If found, seek back and return so main parser can handle it */
    if (is_preprocessor_directive(line)) {
      /* Seek back to before this line */
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    strip_comments(line);
    if (strlen(line) == 0) continue;

    /* Check for new section */
    char section[64];
    if (is_section_header(line, section)) {
      /* Seek back and return so main parser can handle new section */
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    /* Parse atom line: id type resnr residue atom cgnr charge [mass] */
    atom_data_t atom;
    memset(&atom, 0, sizeof(atom));

    int n = sscanf(line, "%d %15s %d %7s %15s %d %f %f",
                   &atom.id, atom.atom_type, &atom.resnr,
                   atom.residue, atom.atom_name, &atom.cgnr,
                   &atom.charge, &atom.mass);

    if (n >= 7) {
      /* Expand array if needed */
      if (mt->natoms >= mt->atoms_allocated) {
        mt->atoms_allocated *= 2;
        mt->atoms = (atom_data_t *)realloc(mt->atoms,
                                           mt->atoms_allocated * sizeof(atom_data_t));
        if (!mt->atoms) return 0;
      }

      mt->atoms[mt->natoms] = atom;
      mt->natoms++;
    }
  }

  return 1;
}

/* Parse [ bonds ] section within a moleculetype */
static int parse_bonds_section(FILE *fp, moltype_t *mt, grotop_data *data, moltype_t **current_mt) {
  char line[GROTOP_RECORD_LENGTH];

  while (1) {
    long pos_before = ftell(fp);  /* Save position BEFORE reading this line */
    if (!fgets(line, sizeof(line), fp)) break;

    /* Check for preprocessor directives BEFORE stripping comments */
    if (is_preprocessor_directive(line)) {
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    strip_comments(line);
    if (strlen(line) == 0) continue;

    /* Check for new section */
    char section[64];
    if (is_section_header(line, section)) {
      /* Seek back and return so main parser can handle new section */
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    /* Parse bond line: ai aj [func params...] */
    int ai, aj;
    if (sscanf(line, "%d %d", &ai, &aj) == 2) {
      /* Expand array if needed */
      if (mt->nbonds >= mt->bonds_allocated) {
        mt->bonds_allocated *= 2;
        mt->bonds = (bond_data_t *)realloc(mt->bonds,
                                           mt->bonds_allocated * sizeof(bond_data_t));
        if (!mt->bonds) return 0;
      }

      mt->bonds[mt->nbonds].ai = ai;
      mt->bonds[mt->nbonds].aj = aj;
      mt->nbonds++;
    }
  }

  return 1;
}

/* Parse [ constraints ] section within a moleculetype - treat as bonds */
static int parse_constraints_section(FILE *fp, moltype_t *mt, grotop_data *data, moltype_t **current_mt) {
  char line[GROTOP_RECORD_LENGTH];

  while (1) {
    long pos_before = ftell(fp);  /* Save position BEFORE reading this line */
    if (!fgets(line, sizeof(line), fp)) break;

    /* Check for preprocessor directives BEFORE stripping comments */
    if (is_preprocessor_directive(line)) {
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    strip_comments(line);
    if (strlen(line) == 0) continue;

    /* Check for new section */
    char section[64];
    if (is_section_header(line, section)) {
      /* Seek back and return so main parser can handle new section */
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    /* Parse constraint line: ai aj [func params...] */
    int ai, aj;
    if (sscanf(line, "%d %d", &ai, &aj) == 2) {
      /* Expand array if needed */
      if (mt->nbonds >= mt->bonds_allocated) {
        mt->bonds_allocated *= 2;
        mt->bonds = (bond_data_t *)realloc(mt->bonds,
                                           mt->bonds_allocated * sizeof(bond_data_t));
        if (!mt->bonds) return 0;
      }

      /* Add constraint as a bond */
      mt->bonds[mt->nbonds].ai = ai;
      mt->bonds[mt->nbonds].aj = aj;
      mt->nbonds++;
    }
  }

  return 1;
}

/* Parse [ angles ] section within a moleculetype */
static int parse_angles_section(FILE *fp, moltype_t *mt, grotop_data *data, moltype_t **current_mt) {
  char line[GROTOP_RECORD_LENGTH];

  while (1) {
    long pos_before = ftell(fp);  /* Save position BEFORE reading this line */
    if (!fgets(line, sizeof(line), fp)) break;

    /* Check for preprocessor directives BEFORE stripping comments */
    if (is_preprocessor_directive(line)) {
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    strip_comments(line);
    if (strlen(line) == 0) continue;

    /* Check for new section */
    char section[64];
    if (is_section_header(line, section)) {
      /* Seek back and return so main parser can handle new section */
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    /* Parse angle line: ai aj ak [func params...] */
    int ai, aj, ak;
    if (sscanf(line, "%d %d %d", &ai, &aj, &ak) == 3) {
      /* Expand array if needed */
      if (mt->nangles >= mt->angles_allocated) {
        mt->angles_allocated *= 2;
        mt->angles = (angle_data_t *)realloc(mt->angles,
                                             mt->angles_allocated * sizeof(angle_data_t));
        if (!mt->angles) return 0;
      }

      mt->angles[mt->nangles].ai = ai;
      mt->angles[mt->nangles].aj = aj;
      mt->angles[mt->nangles].ak = ak;
      mt->nangles++;
    }
  }

  return 1;
}

/* Parse [ dihedrals ] section within a moleculetype */
static int parse_dihedrals_section(FILE *fp, moltype_t *mt, grotop_data *data, moltype_t **current_mt) {
  char line[GROTOP_RECORD_LENGTH];

  while (1) {
    long pos_before = ftell(fp);  /* Save position BEFORE reading this line */
    if (!fgets(line, sizeof(line), fp)) break;

    /* Check for preprocessor directives BEFORE stripping comments */
    if (is_preprocessor_directive(line)) {
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    strip_comments(line);
    if (strlen(line) == 0) continue;

    /* Check for new section */
    char section[64];
    if (is_section_header(line, section)) {
      /* Seek back and return so main parser can handle new section */
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    /* Parse dihedral line: ai aj ak al [func params...] */
    int ai, aj, ak, al, funct;
    /* Try to parse with function type */
    int n = sscanf(line, "%d %d %d %d %d", &ai, &aj, &ak, &al, &funct);
    if (n >= 4) {
      /* Expand array if needed */
      if (mt->ndihedrals >= mt->dihedrals_allocated) {
        mt->dihedrals_allocated *= 2;
        mt->dihedrals = (dihedral_data_t *)realloc(mt->dihedrals,
                                                    mt->dihedrals_allocated * sizeof(dihedral_data_t));
        if (!mt->dihedrals) return 0;
      }

      mt->dihedrals[mt->ndihedrals].ai = ai;
      mt->dihedrals[mt->ndihedrals].aj = aj;
      mt->dihedrals[mt->ndihedrals].ak = ak;
      mt->dihedrals[mt->ndihedrals].al = al;
      mt->dihedrals[mt->ndihedrals].funct = (n == 5) ? funct : 0;
      mt->ndihedrals++;
    }
  }

  return 1;
}

/* Parse [ molecules ] section */
static int parse_molecules_section(FILE *fp, grotop_data *data) {
  char line[GROTOP_RECORD_LENGTH];

  printf("grotopplugin) Parsing [molecules] section\n");

  while (1) {
    long pos_before = ftell(fp);  /* Save position BEFORE reading this line */
    if (!fgets(line, sizeof(line), fp)) break;

    /* Check for preprocessor directives BEFORE stripping comments */
    if (is_preprocessor_directive(line)) {
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    strip_comments(line);
    if (strlen(line) == 0) continue;

    /* Check for new section */
    char section[64];
    if (is_section_header(line, section)) {
      fseek(fp, pos_before, SEEK_SET);
      return 1;
    }

    /* Parse: molname count */
    char molname[32];
    int count;
    if (sscanf(line, "%31s %d", molname, &count) == 2) {
      printf("grotopplugin)   Found molecule: %s x %d\n", molname, count);
      if (data->num_molecules < MAX_MOLECULES) {
        strncpy(data->molnames[data->num_molecules], molname, 31);
        data->molnames[data->num_molecules][31] = '\0';
        data->molcounts[data->num_molecules] = count;
        data->num_molecules++;
      }
    }
  }

  return 1;
}

/* Parse a topology file (recursively handles includes) */
static int parse_topology_file(const char *filepath, grotop_data *data, int depth);

/* Process current section based on name */
static int process_section(FILE *fp, const char *section, grotop_data *data, moltype_t **current_mt) {
  printf("grotopplugin) Processing section: [%s]\n", section);

  if (strcmp(section, "atomtypes") == 0) {
    return parse_atomtypes_section(fp, data);
  }
  else if (strcmp(section, "moleculetype") == 0) {
    moltype_t *mt = create_moltype();
    if (!mt) return 0;

    if (!parse_moleculetype_header(fp, mt)) {
      free_moltype(mt);
      return 0;
    }

    if (data->num_moltypes < MAX_MOLTYPES) {
      data->moltypes[data->num_moltypes++] = mt;
      *current_mt = mt;
      return 1;
    }

    fprintf(stderr, "grotopplugin) WARNING: Maximum number of molecule types (%d) exceeded\n", MAX_MOLTYPES);
    free_moltype(mt);
    return 0;
  }
  else if (strcmp(section, "atoms") == 0 && *current_mt) {
    return parse_atoms_section(fp, *current_mt, data, current_mt);
  }
  else if (strcmp(section, "bonds") == 0 && *current_mt) {
    return parse_bonds_section(fp, *current_mt, data, current_mt);
  }
  else if (strcmp(section, "constraints") == 0 && *current_mt) {
    return parse_constraints_section(fp, *current_mt, data, current_mt);
  }
  else if (strcmp(section, "angles") == 0 && *current_mt) {
    return parse_angles_section(fp, *current_mt, data, current_mt);
  }
  else if (strcmp(section, "dihedrals") == 0 && *current_mt) {
    return parse_dihedrals_section(fp, *current_mt, data, current_mt);
  }
  else if (strcmp(section, "molecules") == 0) {
    return parse_molecules_section(fp, data);
  }
  else if (strcmp(section, "system") == 0 ||
           strcmp(section, "defaults") == 0 ||
           strcmp(section, "pairs") == 0 ||
           strcmp(section, "exclusions") == 0 ||
           strcmp(section, "settles") == 0 ||
           strcmp(section, "position_restraints") == 0) {
    /* Skip these sections for now */
    char line[GROTOP_RECORD_LENGTH];
    while (1) {
      long pos_before = ftell(fp);  /* Save position BEFORE reading this line */
      if (!fgets(line, sizeof(line), fp)) break;

      /* Check for preprocessor directives BEFORE stripping comments */
      if (is_preprocessor_directive(line)) {
        fseek(fp, pos_before, SEEK_SET);
        return 1;
      }

      strip_comments(line);
      if (strlen(line) == 0) continue;

      char next_section[64];
      if (is_section_header(line, next_section)) {
        /* Found next section - seek back and return so main parser can handle it */
        fseek(fp, pos_before, SEEK_SET);
        return 1;
      }
    }
  }

  return 1;
}

/* Parse a topology file (recursively handles includes) */
static int parse_topology_file(const char *filepath, grotop_data *data, int depth) {
  if (depth > MAX_INCLUDES) {
    fprintf(stderr, "grotopplugin) Too many nested includes (depth %d)\n", depth);
    return 0;
  }

  FILE *fp = fopen(filepath, "r");
  if (!fp) {
    fprintf(stderr, "grotopplugin) Cannot open file '%s': %s\n", filepath, strerror(errno));
    return 0;
  }

  printf("grotopplugin) Parsing file: %s (depth %d)\n", filepath, depth);

  char line[GROTOP_RECORD_LENGTH];
  moltype_t *current_mt = NULL;
  char current_section[64] = "";  /* Track current section we're parsing */

  /* Conditional compilation state stack */
  int ifdef_stack[MAX_IFDEF_DEPTH];  /* 1 if condition is true, 0 if false */
  int ifdef_depth = 0;

  while (fgets(line, sizeof(line), fp)) {
    /* Check for #define directive */
    if (parse_define(line, data)) {
      /* If we're in a section, re-enter it after processing the directive */
      if (strlen(current_section) > 0) {
        if (!process_section(fp, current_section, data, &current_mt)) {
          fclose(fp);
          return 0;
        }
      }
      continue;
    }

    /* Check for #ifdef / #ifndef directive */
    char symbol[64];
    int is_ifndef;
    if (parse_ifdef(line, symbol, &is_ifndef)) {
      if (ifdef_depth >= MAX_IFDEF_DEPTH) {
        fprintf(stderr, "grotopplugin) ERROR: Too many nested #ifdef directives (max %d)\n", MAX_IFDEF_DEPTH);
        fclose(fp);
        return 0;
      }

      int condition = is_defined(data, symbol);
      if (is_ifndef) {
        condition = !condition;  /* #ifndef is the inverse of #ifdef */
      }

      ifdef_stack[ifdef_depth] = condition;
      ifdef_depth++;

      printf("grotopplugin) %s %s -> %s\n",
             is_ifndef ? "#ifndef" : "#ifdef",
             symbol,
             condition ? "true (processing)" : "false (skipping)");

      /* If we're in a section, re-enter it after processing the directive */
      if (strlen(current_section) > 0) {
        if (!process_section(fp, current_section, data, &current_mt)) {
          fclose(fp);
          return 0;
        }
      }
      continue;
    }

    /* Check for #else directive */
    if (is_else_directive(line)) {
      if (ifdef_depth == 0) {
        fprintf(stderr, "grotopplugin) ERROR: #else without matching #ifdef\n");
        fclose(fp);
        return 0;
      }

      /* Flip the condition */
      ifdef_stack[ifdef_depth - 1] = !ifdef_stack[ifdef_depth - 1];
      printf("grotopplugin) #else -> %s\n",
             ifdef_stack[ifdef_depth - 1] ? "true (processing)" : "false (skipping)");

      /* If we're in a section, re-enter it after processing the directive */
      if (strlen(current_section) > 0) {
        if (!process_section(fp, current_section, data, &current_mt)) {
          fclose(fp);
          return 0;
        }
      }
      continue;
    }

    /* Check for #endif directive */
    if (is_endif_directive(line)) {
      if (ifdef_depth == 0) {
        fprintf(stderr, "grotopplugin) ERROR: #endif without matching #ifdef\n");
        fclose(fp);
        return 0;
      }

      ifdef_depth--;
      printf("grotopplugin) #endif (depth now %d)\n", ifdef_depth);

      /* If we're in a section, re-enter it after processing the directive */
      if (strlen(current_section) > 0) {
        if (!process_section(fp, current_section, data, &current_mt)) {
          fclose(fp);
          return 0;
        }
      }
      continue;
    }

    /* Check if we should process this line based on conditional compilation state */
    int should_process = 1;
    for (int i = 0; i < ifdef_depth; i++) {
      if (!ifdef_stack[i]) {
        should_process = 0;
        break;
      }
    }

    if (!should_process) {
      /* Skip this line - we're in a false conditional block */
      continue;
    }

    /* Handle includes */
    char include_path[512];
    if (parse_include(line, filepath, include_path)) {
      /* Parse included file without closing this one */
      if (!parse_topology_file(include_path, data, depth + 1)) {
        fclose(fp);
        return 0;
      }
      /* Continue reading current file */
      continue;
    }

    strip_comments(line);
    if (strlen(line) == 0) continue;

    /* Check for section header */
    char section[64];
    if (is_section_header(line, section)) {
      /* New section encountered - clear previous section and enter new one */
      strncpy(current_section, section, 63);
      current_section[63] = '\0';
      if (!process_section(fp, section, data, &current_mt)) {
        fclose(fp);
        return 0;
      }
      /* Keep current_section set - it will be cleared when we encounter the next section */
    }
    else if (strlen(current_section) > 0) {
      /* We have a current section but this line is not a section header
       * This shouldn't happen in normal flow - this line was already processed
       * Skip it */
    }
  }

  /* Check for unmatched #ifdef */
  if (ifdef_depth != 0) {
    fprintf(stderr, "grotopplugin) WARNING: %d unmatched #ifdef directive(s) in file %s\n",
            ifdef_depth, filepath);
  }

  fclose(fp);
  return 1;
}

/* Calculate total atoms from molecules section */
static int calculate_total_atoms(grotop_data *data) {
  int total = 0;

  for (int i = 0; i < data->num_molecules; i++) {
    moltype_t *mt = find_moltype(data, data->molnames[i]);
    if (!mt) {
      fprintf(stderr, "grotopplugin) Unknown molecule type '%s' in [molecules] section\n",
              data->molnames[i]);
      return -1;
    }
    total += mt->natoms * data->molcounts[i];
  }

  return total;
}

/* Calculate total bonds from molecules section */
static int calculate_total_bonds(grotop_data *data) {
  int total = 0;

  for (int i = 0; i < data->num_molecules; i++) {
    moltype_t *mt = find_moltype(data, data->molnames[i]);
    if (mt) {
      total += mt->nbonds * data->molcounts[i];
    }
  }

  return total;
}

/* Calculate total angles from molecules section */
static int calculate_total_angles(grotop_data *data) {
  int total = 0;

  for (int i = 0; i < data->num_molecules; i++) {
    moltype_t *mt = find_moltype(data, data->molnames[i]);
    if (mt) {
      total += mt->nangles * data->molcounts[i];
    }
  }

  return total;
}

/* Calculate total dihedrals from molecules section */
static int calculate_total_dihedrals(grotop_data *data) {
  int total = 0;

  for (int i = 0; i < data->num_molecules; i++) {
    moltype_t *mt = find_moltype(data, data->molnames[i]);
    if (mt) {
      total += mt->ndihedrals * data->molcounts[i];
    }
  }

  return total;
}

/* Calculate total impropers from molecules section (counted from dihedrals with funct type 2 or 4) */
static int calculate_total_impropers(grotop_data *data) {
  int total = 0;

  for (int i = 0; i < data->num_molecules; i++) {
    moltype_t *mt = find_moltype(data, data->molnames[i]);
    if (mt) {
      /* Count improper dihedrals (function types 2 and 4) */
      for (int j = 0; j < mt->ndihedrals; j++) {
        if (mt->dihedrals[j].funct == 2 || mt->dihedrals[j].funct == 4) {
          total += data->molcounts[i];
        }
      }
    }
  }

  return total;
}


/*
 * Main Plugin API Functions
 */

static void *open_grotop_read(const char *filepath, const char *filetype, int *natoms) {
  grotop_data *data = (grotop_data *)calloc(1, sizeof(grotop_data));
  if (!data) return NULL;

  strncpy(data->filepath, filepath, sizeof(data->filepath) - 1);

  /* Parse the topology file */
  if (!parse_topology_file(filepath, data, 0)) {
    fprintf(stderr, "grotopplugin) Failed to parse topology file\n");
    free(data);
    return NULL;
  }

  /* Calculate total atoms */
  data->total_atoms = calculate_total_atoms(data);
  if (data->total_atoms < 0) {
    fprintf(stderr, "grotopplugin) Failed to calculate total atoms\n");
    free(data);
    return NULL;
  }

  data->total_bonds = calculate_total_bonds(data);
  data->total_angles = calculate_total_angles(data);
  data->total_dihedrals = calculate_total_dihedrals(data);
  data->total_impropers = calculate_total_impropers(data);

  *natoms = data->total_atoms;

  printf("grotopplugin) Parsed %d molecule types, %d atoms total, %d bonds, %d angles, %d dihedrals, %d impropers\n",
         data->num_moltypes, data->total_atoms, data->total_bonds,
         data->total_angles, data->total_dihedrals, data->total_impropers);

  return data;
}

static int read_grotop_structure(void *mydata, int *optflags, molfile_atom_t *atoms) {
  grotop_data *data = (grotop_data *)mydata;

  *optflags = MOLFILE_CHARGE | MOLFILE_MASS;

  int global_atom_idx = 0;
  int residue_offset = 0;  /* Track residue numbering across all molecules */

  /* Instantiate molecules */
  for (int mol_idx = 0; mol_idx < data->num_molecules; mol_idx++) {
    moltype_t *mt = find_moltype(data, data->molnames[mol_idx]);
    if (!mt) continue;

    int count = data->molcounts[mol_idx];

    /* Generate segment ID - same for all copies of this molecule type */
    /* Truncate to 4 characters and convert to uppercase (match Python behavior) */
    char segid[8];
    snprintf(segid, sizeof(segid), "%.4s", mt->name);
    for (int j = 0; segid[j]; j++) {
      segid[j] = toupper(segid[j]);
    }

    /* Find min and max residue numbers in this molecule type */
    int min_resid = mt->natoms > 0 ? mt->atoms[0].resnr : 1;
    int max_resid = min_resid;
    for (int i = 0; i < mt->natoms; i++) {
      if (mt->atoms[i].resnr < min_resid) min_resid = mt->atoms[i].resnr;
      if (mt->atoms[i].resnr > max_resid) max_resid = mt->atoms[i].resnr;
    }
    int num_residues = max_resid - min_resid + 1;

    /* Create multiple copies */
    for (int copy = 0; copy < count; copy++) {
      /* Calculate residue offset for this molecule instance */
      int resid_offset = residue_offset - min_resid + 1;

      /* Copy atoms */
      for (int i = 0; i < mt->natoms; i++) {
        atom_data_t *src = &mt->atoms[i];
        molfile_atom_t *dst = &atoms[global_atom_idx];

        strncpy(dst->name, src->atom_name, sizeof(dst->name) - 1);
        dst->name[sizeof(dst->name) - 1] = '\0';

        strncpy(dst->type, src->atom_type, sizeof(dst->type) - 1);
        dst->type[sizeof(dst->type) - 1] = '\0';

        strncpy(dst->resname, src->residue, sizeof(dst->resname) - 1);
        dst->resname[sizeof(dst->resname) - 1] = '\0';

        /* Adjust residue ID to be continuous across all molecule copies */
        dst->resid = src->resnr + resid_offset;

        strncpy(dst->segid, segid, sizeof(dst->segid) - 1);
        dst->segid[sizeof(dst->segid) - 1] = '\0';

        dst->chain[0] = '\0';

        dst->charge = src->charge;

        /* Use mass from atom or lookup from atomtypes */
        if (src->mass > 0.0f) {
          dst->mass = src->mass;
        } else {
          dst->mass = find_atomtype_mass(data, src->atom_type);
        }

        global_atom_idx++;
      }

      /* Update residue offset for next molecule copy */
      residue_offset += num_residues;
    }
  }

  return MOLFILE_SUCCESS;
}

static int read_grotop_bonds(void *mydata, int *nbonds, int **fromptr, int **toptr,
                             float **bondorderptr, int **bondtype,
                             int *nbondtypes, char ***bondtypename) {
  grotop_data *data = (grotop_data *)mydata;

  if (data->total_bonds == 0) {
    *nbonds = 0;
    *fromptr = NULL;
    *toptr = NULL;
    *bondorderptr = NULL;
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;
    return MOLFILE_SUCCESS;
  }

  /* Allocate bond arrays */
  data->bond_from = (int *)malloc(data->total_bonds * sizeof(int));
  data->bond_to = (int *)malloc(data->total_bonds * sizeof(int));

  if (!data->bond_from || !data->bond_to) {
    if (data->bond_from) free(data->bond_from);
    if (data->bond_to) free(data->bond_to);
    return MOLFILE_ERROR;
  }

  int bond_idx = 0;
  int global_atom_offset = 0;

  /* Instantiate bonds */
  for (int mol_idx = 0; mol_idx < data->num_molecules; mol_idx++) {
    moltype_t *mt = find_moltype(data, data->molnames[mol_idx]);
    if (!mt) continue;

    int count = data->molcounts[mol_idx];

    for (int copy = 0; copy < count; copy++) {
      for (int i = 0; i < mt->nbonds; i++) {
        /* Convert to 1-based global indices */
        data->bond_from[bond_idx] = global_atom_offset + mt->bonds[i].ai;
        data->bond_to[bond_idx] = global_atom_offset + mt->bonds[i].aj;
        bond_idx++;
      }

      global_atom_offset += mt->natoms;
    }
  }

  *nbonds = data->total_bonds;
  *fromptr = data->bond_from;
  *toptr = data->bond_to;
  *bondorderptr = NULL;
  *bondtype = NULL;
  *nbondtypes = 0;
  *bondtypename = NULL;

  return MOLFILE_SUCCESS;
}

static int read_grotop_angles(void *handle, int *numangles, int **angles, int **angletypes,
                               int *numangletypes, char ***angletypenames, int *numdihedrals,
                               int **dihedrals, int **dihedraltypes, int *numdihedraltypes,
                               char ***dihedraltypenames, int *numimpropers, int **impropers,
                               int **impropertypes, int *numimpropertypes, char ***impropertypenames,
                               int *numcterms, int **cterms, int *ctermcols, int *ctermrows) {
  grotop_data *data = (grotop_data *)handle;

  /* Initialize all pointers to NULL/0 */
  *numangles = 0;
  *angles = NULL;
  *angletypes = NULL;
  *numangletypes = 0;
  *angletypenames = NULL;

  *numdihedrals = 0;
  *dihedrals = NULL;
  *dihedraltypes = NULL;
  *numdihedraltypes = 0;
  *dihedraltypenames = NULL;

  *numimpropers = 0;
  *impropers = NULL;
  *impropertypes = NULL;
  *numimpropertypes = 0;
  *impropertypenames = NULL;

  *numcterms = 0;
  *cterms = NULL;
  *ctermcols = 0;
  *ctermrows = 0;

  /* Allocate angle arrays if needed */
  if (data->total_angles > 0) {
    data->angles = (int *)malloc(data->total_angles * 3 * sizeof(int));
    if (!data->angles) return MOLFILE_ERROR;

    int angle_idx = 0;
    int global_atom_offset = 0;

    /* Instantiate angles */
    for (int mol_idx = 0; mol_idx < data->num_molecules; mol_idx++) {
      moltype_t *mt = find_moltype(data, data->molnames[mol_idx]);
      if (!mt) continue;

      int count = data->molcounts[mol_idx];

      for (int copy = 0; copy < count; copy++) {
        for (int i = 0; i < mt->nangles; i++) {
          data->angles[angle_idx * 3 + 0] = global_atom_offset + mt->angles[i].ai;
          data->angles[angle_idx * 3 + 1] = global_atom_offset + mt->angles[i].aj;
          data->angles[angle_idx * 3 + 2] = global_atom_offset + mt->angles[i].ak;
          angle_idx++;
        }

        global_atom_offset += mt->natoms;
      }
    }

    *numangles = data->total_angles;
    *angles = data->angles;
  }

  /* Allocate dihedral arrays if needed */
  if (data->total_dihedrals > 0) {
    data->dihedrals = (int *)malloc(data->total_dihedrals * 4 * sizeof(int));
    if (!data->dihedrals) {
      if (data->angles) free(data->angles);
      return MOLFILE_ERROR;
    }

    int dihedral_idx = 0;
    int global_atom_offset = 0;

    /* Instantiate dihedrals (only proper dihedrals, not impropers) */
    for (int mol_idx = 0; mol_idx < data->num_molecules; mol_idx++) {
      moltype_t *mt = find_moltype(data, data->molnames[mol_idx]);
      if (!mt) continue;

      int count = data->molcounts[mol_idx];

      for (int copy = 0; copy < count; copy++) {
        for (int i = 0; i < mt->ndihedrals; i++) {
          /* Skip improper dihedrals (function types 2 and 4) */
          if (mt->dihedrals[i].funct == 2 || mt->dihedrals[i].funct == 4) continue;

          data->dihedrals[dihedral_idx * 4 + 0] = global_atom_offset + mt->dihedrals[i].ai;
          data->dihedrals[dihedral_idx * 4 + 1] = global_atom_offset + mt->dihedrals[i].aj;
          data->dihedrals[dihedral_idx * 4 + 2] = global_atom_offset + mt->dihedrals[i].ak;
          data->dihedrals[dihedral_idx * 4 + 3] = global_atom_offset + mt->dihedrals[i].al;
          dihedral_idx++;
        }

        global_atom_offset += mt->natoms;
      }
    }

    *numdihedrals = dihedral_idx;
    *dihedrals = data->dihedrals;
  }

  /* Allocate improper arrays if needed */
  if (data->total_impropers > 0) {
    data->impropers = (int *)malloc(data->total_impropers * 4 * sizeof(int));
    if (!data->impropers) {
      if (data->angles) free(data->angles);
      if (data->dihedrals) free(data->dihedrals);
      return MOLFILE_ERROR;
    }

    int improper_idx = 0;
    int global_atom_offset = 0;

    /* Instantiate impropers (function types 2 and 4) */
    for (int mol_idx = 0; mol_idx < data->num_molecules; mol_idx++) {
      moltype_t *mt = find_moltype(data, data->molnames[mol_idx]);
      if (!mt) continue;

      int count = data->molcounts[mol_idx];

      for (int copy = 0; copy < count; copy++) {
        for (int i = 0; i < mt->ndihedrals; i++) {
          /* Only process improper dihedrals (function types 2 and 4) */
          if (mt->dihedrals[i].funct != 2 && mt->dihedrals[i].funct != 4) continue;

          data->impropers[improper_idx * 4 + 0] = global_atom_offset + mt->dihedrals[i].ai;
          data->impropers[improper_idx * 4 + 1] = global_atom_offset + mt->dihedrals[i].aj;
          data->impropers[improper_idx * 4 + 2] = global_atom_offset + mt->dihedrals[i].ak;
          data->impropers[improper_idx * 4 + 3] = global_atom_offset + mt->dihedrals[i].al;
          improper_idx++;
        }

        global_atom_offset += mt->natoms;
      }
    }

    *numimpropers = improper_idx;
    *impropers = data->impropers;
  }

  return MOLFILE_SUCCESS;
}

static void close_grotop_read(void *mydata) {
  grotop_data *data = (grotop_data *)mydata;
  if (!data) return;

  /* Free molecule types */
  for (int i = 0; i < data->num_moltypes; i++) {
    free_moltype(data->moltypes[i]);
  }

  /* Free bond arrays */
  if (data->bond_from) free(data->bond_from);
  if (data->bond_to) free(data->bond_to);

  /* Free angle/dihedral/improper arrays */
  if (data->angles) free(data->angles);
  if (data->dihedrals) free(data->dihedrals);
  if (data->impropers) free(data->impropers);

  free(data);
}


/*
 * Plugin Registration
 */

static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "grotop";
  plugin.prettyname = "GROMACS Topology";
  plugin.author = "Generated with Claude Code";
  plugin.majorv = 0;
  plugin.minorv = 1;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "top,itp";
  plugin.open_file_read = open_grotop_read;
  plugin.read_structure = read_grotop_structure;
  plugin.read_bonds = read_grotop_bonds;
  plugin.read_angles = read_grotop_angles;
  plugin.close_file_read = close_grotop_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}
