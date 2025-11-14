# GROMACS to PSF Converter

A Python tool to convert GROMACS topology files (.top, .itp) to PSF (Protein Structure File) format. This tool is designed with a modular architecture that facilitates future development of a VMD molfile plugin in C.

## Features

- Reads GROMACS .top (topology) files
- Automatically processes #include directives for .itp files
- Supports nested includes
- Handles multiple molecule types and instances
- Generates PSF files compatible with VMD and NAMD
- Clean, modular code design for easy translation to C

## Architecture

The code is organized into distinct modules that mirror what would be needed for a VMD molfile plugin:

1. **Data Structures**: Clear classes (Atom, MoleculeType, Topology) that map to C structs
2. **Parsers**: Separate functions for parsing different file types and sections
3. **System Builder**: Instantiates molecules from templates
4. **PSF Writer**: Outputs standard PSF format

## GROMACS Topology Structure

GROMACS topologies typically consist of:
- **`.top` file**: Main topology file that defines the system
  - Contains `[ system ]` section with system name
  - Contains `[ molecules ]` section listing molecule instances
  - Includes force field and molecule definition files via `#include`

- **`.itp` files**: Include topology files with molecule definitions
  - `[ moleculetype ]`: Defines a molecule type
  - `[ atoms ]`: Atom definitions within the molecule
  - `[ bonds ]`, `[ angles ]`, `[ dihedrals ]`: Bonded interactions
  - Can be nested (one .itp can include others)

## Usage

```bash
# Basic usage
python gromacs_to_psf.py -s topol.top -o system.psf

# With additional .itp files
python gromacs_to_psf.py -s topol.top -p protein.itp lipid.itp -o system.psf

# Verbose output
python gromacs_to_psf.py -s topol.top -o system.psf -v
```

## Command-line Arguments

- `-s, --top`: Input GROMACS .top file (required)
- `-p, --itp`: Additional .itp files (optional, can specify multiple)
- `-o, --output`: Output PSF file (required)
- `-v, --verbose`: Enable verbose output

## Example GROMACS Files

### Minimal .top file (topol.top)
```
#include "forcefield.itp"
#include "protein.itp"

[ system ]
; Name
Protein in water

[ molecules ]
; Compound        #mols
Protein             1
```

### Minimal .itp file (protein.itp)
```
[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass
     1         CA      1    ALA     CA      1      0.000     12.010
     2          C      1    ALA      C      1      0.500     12.010
     3          O      1    ALA      O      1     -0.500     15.999

[ bonds ]
;  ai    aj funct
    1     2     1
    2     3     1
```

## PSF Output Format

The generated PSF file contains:
- **!NATOM**: Atom definitions with segment IDs, residue info, charges, masses
- **!NBOND**: Bond connectivity
- **!NTHETA**: Angle definitions
- **!NPHI**: Dihedral definitions
- **!NIMPHI**: Improper dihedral definitions

## Design Considerations for VMD Plugin

The code is structured to facilitate porting to a VMD molfile plugin:

1. **Modular functions**: Each parsing function has a single responsibility
2. **Clear data flow**: Parse → Build → Write
3. **No complex Python features**: Uses basic data structures and control flow
4. **Comment annotations**: Indicates C struct equivalents

### Future VMD Plugin Development

Key components for a VMD plugin:
- `open_file_read()`: Initialize and parse topology
- `read_structure()`: Return atom data (coordinates would come from .gro)
- `read_bonds()`: Return bond topology
- `close_file_read()`: Cleanup

The current Python implementation can serve as a reference for the C implementation.

## Limitations

- Does not read coordinate information (use .gro file separately in VMD)
- Does not process force field parameters (only topology)
- Assumes sequential atom numbering within molecules
- Does not handle all GROMACS topology features (e.g., virtual sites, constraints)

## Testing

To test the converter, you need:
1. A GROMACS .top file
2. Associated .itp files (if not already included)
3. Run the conversion and load the PSF in VMD

```bash
# Run conversion
python gromacs_to_psf.py -s your_system.top -o output.psf

# View in VMD (you'll also need a .gro or .pdb for coordinates)
vmd your_system.gro output.psf
```

## Dependencies

- Python 3.7+
- No external packages required (uses only standard library)

## Author

Generated with Claude Code

## License

Free to use and modify for research purposes.
