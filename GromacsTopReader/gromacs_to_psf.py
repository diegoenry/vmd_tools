#!/usr/bin/env python3
"""
GROMACS to PSF Converter

This module converts GROMACS topology files (.top, .itp) to PSF format.
Designed with modularity in mind for future VMD molfile plugin implementation in C.

Author: Diego E.B. Gomes
"""

import re
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass, field


# =============================================================================
# Data Structures (analogous to C structs for VMD plugin)
# =============================================================================

@dataclass
class Atom:
    """Atom data structure - maps to PSF atom record"""
    atom_id: int = 0
    atom_type: str = ""
    residue_number: int = 0
    residue_name: str = ""
    atom_name: str = ""
    charge_group: int = 0
    charge: float = 0.0
    mass: float = 0.0
    segment_id: str = ""

    # GROMACS specific (not in PSF)
    original_id: int = 0  # ID within molecule type
    molecule_type: str = ""


@dataclass
class MoleculeType:
    """Molecule type definition from .itp file"""
    name: str = ""
    nrexcl: int = 3
    atoms: List[Atom] = field(default_factory=list)
    bonds: List[Tuple[int, int]] = field(default_factory=list)
    angles: List[Tuple[int, int, int]] = field(default_factory=list)
    dihedrals: List[Tuple[int, int, int, int]] = field(default_factory=list)
    impropers: List[Tuple[int, int, int, int]] = field(default_factory=list)


@dataclass
class Topology:
    """Complete topology data structure"""
    system_name: str = "SYSTEM"
    molecule_types: Dict[str, MoleculeType] = field(default_factory=dict)
    molecules: List[Tuple[str, int]] = field(default_factory=list)  # (name, count)
    atom_type_masses: Dict[str, float] = field(default_factory=dict)  # atom_type -> mass

    # Preprocessor state
    defines: set = field(default_factory=set)  # Defined symbols

    # Instantiated system
    atoms: List[Atom] = field(default_factory=list)
    bonds: List[Tuple[int, int]] = field(default_factory=list)
    angles: List[Tuple[int, int, int]] = field(default_factory=list)
    dihedrals: List[Tuple[int, int, int, int]] = field(default_factory=list)
    impropers: List[Tuple[int, int, int, int]] = field(default_factory=list)


# =============================================================================
# File Parsing Functions
# =============================================================================

def strip_comments(line: str) -> str:
    """Remove comments from a line (everything after semicolon)"""
    if ';' in line:
        line = line[:line.index(';')]
    return line.strip()


def is_section_header(line: str) -> Optional[str]:
    """Check if line is a section header and return section name"""
    line = line.strip()
    if line.startswith('[') and line.endswith(']'):
        return line[1:-1].strip()
    return None


def parse_include_directive(line: str, base_path: Path) -> Optional[Path]:
    """Parse #include directive and return the file path"""
    if not line.startswith('#include'):
        return None

    match = re.search(r'["\'](.+?)["\']', line)
    if match:
        include_file = match.group(1)
        return base_path / include_file
    return None


def parse_define_directive(line: str) -> Optional[str]:
    """
    Parse #define directive and return the symbol name
    Returns: symbol name or None
    """
    if not line.startswith('#define'):
        return None

    parts = line.split()
    if len(parts) >= 2:
        return parts[1]  # Return symbol name
    return None


def parse_ifdef_directive(line: str) -> Optional[Tuple[str, bool]]:
    """
    Parse #ifdef or #ifndef directive
    Returns: (symbol, is_ifndef) or None
    """
    line_stripped = line.strip()

    if line_stripped.startswith('#ifdef'):
        parts = line_stripped.split()
        if len(parts) >= 2:
            return (parts[1], False)  # (symbol, is_ifndef=False)
    elif line_stripped.startswith('#ifndef'):
        parts = line_stripped.split()
        if len(parts) >= 2:
            return (parts[1], True)  # (symbol, is_ifndef=True)

    return None


def is_else_directive(line: str) -> bool:
    """Check if line is #else directive"""
    return line.strip().startswith('#else')


def is_endif_directive(line: str) -> bool:
    """Check if line is #endif directive"""
    return line.strip().startswith('#endif')


def parse_atomtype_line(line: str) -> Optional[Tuple[str, float]]:
    """
    Parse an atomtype line from [ atomtypes ] section
    Format: type mass charge ptype sigma epsilon
    Returns: (atom_type, mass) or None
    """
    parts = line.split()
    if len(parts) < 2:
        return None

    try:
        atom_type = parts[0]
        mass = float(parts[1])
        return (atom_type, mass)
    except (ValueError, IndexError):
        return None


def parse_atom_line(line: str, atom_type_masses: Dict[str, float]) -> Optional[Atom]:
    """
    Parse an atom line from [ atoms ] section
    Format: nr type resnr residue atom cgnr charge [mass]
    If mass is not provided, look it up in atom_type_masses dictionary
    """
    parts = line.split()
    if len(parts) < 7:
        return None

    atom_type = parts[1]
    mass = 0.0

    # Check if mass is explicitly provided in atom definition
    if len(parts) > 7:
        mass = float(parts[7])
    # Otherwise, look it up from atomtypes
    elif atom_type in atom_type_masses:
        mass = atom_type_masses[atom_type]
    else:
        # If not found, default to 0.0 and warn
        print(f"Warning: Mass not found for atom type '{atom_type}', using 0.0")

    atom = Atom(
        original_id=int(parts[0]),
        atom_type=atom_type,
        residue_number=int(parts[2]),
        residue_name=parts[3],
        atom_name=parts[4],
        charge_group=int(parts[5]),
        charge=float(parts[6]),
        mass=mass
    )
    return atom


def parse_bond_line(line: str) -> Optional[Tuple[int, int]]:
    """Parse a bond line from [ bonds ] section"""
    parts = line.split()
    if len(parts) < 2:
        return None
    return (int(parts[0]), int(parts[1]))


def parse_constraint_line(line: str) -> Optional[Tuple[int, int]]:
    """
    Parse a constraint line from [ constraints ] section
    Constraints are treated as bonds (fixed-length connections)
    """
    parts = line.split()
    if len(parts) < 2:
        return None
    return (int(parts[0]), int(parts[1]))


def parse_angle_line(line: str) -> Optional[Tuple[int, int, int]]:
    """Parse an angle line from [ angles ] section"""
    parts = line.split()
    if len(parts) < 3:
        return None
    return (int(parts[0]), int(parts[1]), int(parts[2]))


def parse_dihedral_line(line: str) -> Tuple[Optional[Tuple[int, int, int, int]], bool]:
    """
    Parse a dihedral line from [ dihedrals ] section
    Returns: (dihedral_tuple, is_improper)
    Function types 2 and 4 are impropers in GROMACS
    """
    parts = line.split()
    if len(parts) < 4:
        return (None, False)

    dihedral = (int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]))
    is_improper = len(parts) > 4 and parts[4] in ['2', '4']

    return (dihedral, is_improper)


# =============================================================================
# ITP File Parser
# =============================================================================

def parse_itp_file(filepath: Path, topology: Topology) -> None:
    """
    Parse a GROMACS .itp (include topology) file
    This function can be called recursively for nested includes
    """
    if not filepath.exists():
        print(f"Warning: File not found: {filepath}")
        return

    print(f"  Parsing: {filepath.name}")

    with open(filepath, 'r') as f:
        lines = f.readlines()

    current_section = None
    current_moltype = None

    # Conditional compilation state stack
    ifdef_stack = []  # List of booleans: True if condition is true, False if false

    for line in lines:
        original_line = line
        line_stripped = original_line.strip()

        # Handle #define directive
        symbol = parse_define_directive(line_stripped)
        if symbol:
            topology.defines.add(symbol)
            print(f"    Defined symbol: {symbol}")
            continue

        # Handle #ifdef / #ifndef directive
        ifdef_result = parse_ifdef_directive(line_stripped)
        if ifdef_result:
            symbol, is_ifndef = ifdef_result
            condition = symbol in topology.defines
            if is_ifndef:
                condition = not condition

            ifdef_stack.append(condition)
            print(f"    {'#ifndef' if is_ifndef else '#ifdef'} {symbol} -> {'true (processing)' if condition else 'false (skipping)'}")
            continue

        # Handle #else directive
        if is_else_directive(line_stripped):
            if not ifdef_stack:
                print(f"ERROR: #else without matching #ifdef in {filepath}")
                continue
            # Flip the condition
            ifdef_stack[-1] = not ifdef_stack[-1]
            print(f"    #else -> {'true (processing)' if ifdef_stack[-1] else 'false (skipping)'}")
            continue

        # Handle #endif directive
        if is_endif_directive(line_stripped):
            if not ifdef_stack:
                print(f"ERROR: #endif without matching #ifdef in {filepath}")
                continue
            ifdef_stack.pop()
            print(f"    #endif (depth now {len(ifdef_stack)})")
            continue

        # Check if we should process this line based on conditional compilation state
        should_process = all(ifdef_stack)  # Process only if all conditions are True

        if not should_process:
            # Skip this line - we're in a false conditional block
            continue

        # Check for include directives (only process if not in false conditional)
        include_path = parse_include_directive(line_stripped, filepath.parent)
        if include_path:
            parse_itp_file(include_path, topology)
            continue

        # Now process the line with comments stripped
        line = strip_comments(line)

        if not line:
            continue

        # Check for section headers
        section = is_section_header(line)
        if section:
            current_section = section
            continue

        # Parse content based on current section
        if current_section == 'atomtypes':
            atomtype_data = parse_atomtype_line(line)
            if atomtype_data:
                atom_type, mass = atomtype_data
                topology.atom_type_masses[atom_type] = mass

        elif current_section == 'moleculetype':
            parts = line.split()
            if len(parts) >= 1:
                current_moltype = MoleculeType(
                    name=parts[0],
                    nrexcl=int(parts[1]) if len(parts) > 1 else 3
                )
                topology.molecule_types[current_moltype.name] = current_moltype

        elif current_section == 'atoms' and current_moltype:
            atom = parse_atom_line(line, topology.atom_type_masses)
            if atom:
                atom.molecule_type = current_moltype.name
                current_moltype.atoms.append(atom)

        elif current_section == 'bonds' and current_moltype:
            bond = parse_bond_line(line)
            if bond:
                current_moltype.bonds.append(bond)

        elif current_section == 'constraints' and current_moltype:
            # Constraints are treated as bonds (fixed-length connections)
            constraint = parse_constraint_line(line)
            if constraint:
                current_moltype.bonds.append(constraint)

        elif current_section == 'angles' and current_moltype:
            angle = parse_angle_line(line)
            if angle:
                current_moltype.angles.append(angle)

        elif current_section == 'dihedrals' and current_moltype:
            dihedral, is_improper = parse_dihedral_line(line)
            if dihedral:
                if is_improper:
                    current_moltype.impropers.append(dihedral)
                else:
                    current_moltype.dihedrals.append(dihedral)

    # Check for unmatched #ifdef
    if ifdef_stack:
        print(f"WARNING: {len(ifdef_stack)} unmatched #ifdef directive(s) in {filepath}")


# =============================================================================
# TOP File Parser
# =============================================================================

def parse_top_file(filepath: Path, topology: Optional[Topology] = None) -> Topology:
    """
    Parse main GROMACS .top (topology) file
    This is the entry point for topology parsing
    If topology is provided, it will be used (e.g., with pre-loaded atom types)
    """
    print(f"Reading GROMACS topology: {filepath}")

    if topology is None:
        topology = Topology()

    with open(filepath, 'r') as f:
        lines = f.readlines()

    current_section = None
    current_moltype = None  # Track current molecule type being defined

    # Conditional compilation state stack
    ifdef_stack = []  # List of booleans: True if condition is true, False if false

    for line in lines:
        original_line = line
        line_stripped = original_line.strip()

        # Handle #define directive
        symbol = parse_define_directive(line_stripped)
        if symbol:
            topology.defines.add(symbol)
            print(f"  Defined symbol: {symbol}")
            continue

        # Handle #ifdef / #ifndef directive
        ifdef_result = parse_ifdef_directive(line_stripped)
        if ifdef_result:
            symbol, is_ifndef = ifdef_result
            condition = symbol in topology.defines
            if is_ifndef:
                condition = not condition

            ifdef_stack.append(condition)
            print(f"  {'#ifndef' if is_ifndef else '#ifdef'} {symbol} -> {'true (processing)' if condition else 'false (skipping)'}")
            continue

        # Handle #else directive
        if is_else_directive(line_stripped):
            if not ifdef_stack:
                print(f"ERROR: #else without matching #ifdef in {filepath}")
                continue
            # Flip the condition
            ifdef_stack[-1] = not ifdef_stack[-1]
            print(f"  #else -> {'true (processing)' if ifdef_stack[-1] else 'false (skipping)'}")
            continue

        # Handle #endif directive
        if is_endif_directive(line_stripped):
            if not ifdef_stack:
                print(f"ERROR: #endif without matching #ifdef in {filepath}")
                continue
            ifdef_stack.pop()
            print(f"  #endif (depth now {len(ifdef_stack)})")
            continue

        # Check if we should process this line based on conditional compilation state
        should_process = all(ifdef_stack)  # Process only if all conditions are True

        if not should_process:
            # Skip this line - we're in a false conditional block
            continue

        # Check for include directives (only process if not in false conditional)
        include_path = parse_include_directive(line_stripped, filepath.parent)
        if include_path:
            parse_itp_file(include_path, topology)
            continue

        # Now process the line with comments stripped
        line = strip_comments(line)

        if not line:
            continue

        # Check for section headers
        section = is_section_header(line)
        if section:
            current_section = section
            continue

        # Parse content based on current section
        if current_section == 'atomtypes':
            atomtype_data = parse_atomtype_line(line)
            if atomtype_data:
                atom_type, mass = atomtype_data
                topology.atom_type_masses[atom_type] = mass

        elif current_section == 'moleculetype':
            parts = line.split()
            if len(parts) >= 1:
                current_moltype = MoleculeType(
                    name=parts[0],
                    nrexcl=int(parts[1]) if len(parts) > 1 else 3
                )
                topology.molecule_types[current_moltype.name] = current_moltype

        elif current_section == 'atoms' and current_moltype:
            atom = parse_atom_line(line, topology.atom_type_masses)
            if atom:
                atom.molecule_type = current_moltype.name
                current_moltype.atoms.append(atom)

        elif current_section == 'bonds' and current_moltype:
            bond = parse_bond_line(line)
            if bond:
                current_moltype.bonds.append(bond)

        elif current_section == 'constraints' and current_moltype:
            # Constraints are treated as bonds (fixed-length connections)
            constraint = parse_constraint_line(line)
            if constraint:
                current_moltype.bonds.append(constraint)

        elif current_section == 'angles' and current_moltype:
            angle = parse_angle_line(line)
            if angle:
                current_moltype.angles.append(angle)

        elif current_section == 'dihedrals' and current_moltype:
            dihedral, is_improper = parse_dihedral_line(line)
            if dihedral:
                if is_improper:
                    current_moltype.impropers.append(dihedral)
                else:
                    current_moltype.dihedrals.append(dihedral)

        elif current_section == 'system':
            topology.system_name = line

        elif current_section == 'molecules':
            parts = line.split()
            if len(parts) >= 2:
                mol_name = parts[0]
                mol_count = int(parts[1])
                topology.molecules.append((mol_name, mol_count))

    # Check for unmatched #ifdef
    if ifdef_stack:
        print(f"WARNING: {len(ifdef_stack)} unmatched #ifdef directive(s) in {filepath}")

    return topology


# =============================================================================
# System Builder
# =============================================================================

def make_segment_id(mol_name: str) -> str:
    """
    Generate a clean segment ID from molecule type name.
    PSF segment IDs are typically 4 characters max.
    """
    # Truncate to 4 characters and convert to uppercase
    return mol_name[:4].upper()


def build_system(topology: Topology) -> None:
    """
    Build the complete system by instantiating molecules
    This expands molecule types into actual atoms with proper numbering
    """
    print("\nBuilding system from molecule definitions...")

    atom_offset = 0
    residue_offset = 0  # Track global residue numbering

    for mol_name, mol_count in topology.molecules:
        if mol_name not in topology.molecule_types:
            print(f"Warning: Molecule type '{mol_name}' not found in topology")
            continue

        moltype = topology.molecule_types[mol_name]
        print(f"  Adding {mol_count} x {mol_name} ({len(moltype.atoms)} atoms each)")

        # Use same segment ID for all instances of this molecule type
        segment_id = make_segment_id(mol_name)

        for mol_idx in range(mol_count):
            # Calculate residue offset for this molecule instance
            # Find min residue number in template (usually 1, but could be different)
            if moltype.atoms:
                min_resid = min(a.residue_number for a in moltype.atoms)
                max_resid = max(a.residue_number for a in moltype.atoms)
                resid_offset = residue_offset - min_resid + 1
            else:
                resid_offset = residue_offset

            # Add atoms for this molecule instance
            for atom_template in moltype.atoms:
                new_atom = Atom(
                    atom_id=len(topology.atoms) + 1,
                    atom_type=atom_template.atom_type,
                    residue_number=atom_template.residue_number + resid_offset,
                    residue_name=atom_template.residue_name,
                    atom_name=atom_template.atom_name,
                    charge_group=atom_template.charge_group,
                    charge=atom_template.charge,
                    mass=atom_template.mass,
                    segment_id=segment_id,
                    original_id=atom_template.original_id,
                    molecule_type=mol_name
                )
                topology.atoms.append(new_atom)

            # Add bonds with proper atom offset
            for bond in moltype.bonds:
                new_bond = (bond[0] + atom_offset, bond[1] + atom_offset)
                topology.bonds.append(new_bond)

            # Add angles with proper atom offset
            for angle in moltype.angles:
                new_angle = (
                    angle[0] + atom_offset,
                    angle[1] + atom_offset,
                    angle[2] + atom_offset
                )
                topology.angles.append(new_angle)

            # Add dihedrals with proper atom offset
            for dihedral in moltype.dihedrals:
                new_dihedral = (
                    dihedral[0] + atom_offset,
                    dihedral[1] + atom_offset,
                    dihedral[2] + atom_offset,
                    dihedral[3] + atom_offset
                )
                topology.dihedrals.append(new_dihedral)

            # Add impropers with proper atom offset
            for improper in moltype.impropers:
                new_improper = (
                    improper[0] + atom_offset,
                    improper[1] + atom_offset,
                    improper[2] + atom_offset,
                    improper[3] + atom_offset
                )
                topology.impropers.append(new_improper)

            # Update offsets for next molecule instance
            atom_offset = len(topology.atoms)
            if moltype.atoms:
                # Update residue offset based on the highest residue number used
                num_residues = max_resid - min_resid + 1
                residue_offset += num_residues


# =============================================================================
# PSF Writer
# =============================================================================

def write_psf(topology: Topology, output_path: Path) -> None:
    """Write topology to PSF format file"""
    print(f"\nWriting PSF file: {output_path}")

    with open(output_path, 'w') as f:
        # PSF header
        f.write("PSF CMAP\n\n")
        f.write("       1 !NTITLE\n")
        f.write(f" REMARKS {topology.system_name}\n\n")

        # Atoms section
        natoms = len(topology.atoms)
        f.write(f"{natoms:8d} !NATOM\n")

        for atom in topology.atoms:
            f.write(f"{atom.atom_id:8d} "
                   f"{atom.segment_id:4s} "
                   f"{atom.residue_number:<4d} "
                   f"{atom.residue_name:4s} "
                   f"{atom.atom_name:4s} "
                   f"{atom.atom_type:6s} "
                   f"{atom.charge:10.6f} "
                   f"{atom.mass:13.4f} "
                   f"           0\n")

        f.write("\n")

        # Bonds section
        nbonds = len(topology.bonds)
        f.write(f"{nbonds:8d} !NBOND: bonds\n")

        for i, bond in enumerate(topology.bonds):
            f.write(f"{bond[0]:8d}{bond[1]:8d}")
            if (i + 1) % 4 == 0:
                f.write("\n")
        if nbonds % 4 != 0:
            f.write("\n")

        f.write("\n")

        # Angles section
        nangles = len(topology.angles)
        f.write(f"{nangles:8d} !NTHETA: angles\n")

        for i, angle in enumerate(topology.angles):
            f.write(f"{angle[0]:8d}{angle[1]:8d}{angle[2]:8d}")
            if (i + 1) % 3 == 0:
                f.write("\n")
        if nangles % 3 != 0:
            f.write("\n")

        f.write("\n")

        # Dihedrals section
        ndihedrals = len(topology.dihedrals)
        f.write(f"{ndihedrals:8d} !NPHI: dihedrals\n")

        for i, dihedral in enumerate(topology.dihedrals):
            f.write(f"{dihedral[0]:8d}{dihedral[1]:8d}"
                   f"{dihedral[2]:8d}{dihedral[3]:8d}")
            if (i + 1) % 2 == 0:
                f.write("\n")
        if ndihedrals % 2 != 0:
            f.write("\n")

        f.write("\n")

        # Impropers section
        nimpropers = len(topology.impropers)
        f.write(f"{nimpropers:8d} !NIMPHI: impropers\n")

        for i, improper in enumerate(topology.impropers):
            f.write(f"{improper[0]:8d}{improper[1]:8d}"
                   f"{improper[2]:8d}{improper[3]:8d}")
            if (i + 1) % 2 == 0:
                f.write("\n")
        if nimpropers % 2 != 0:
            f.write("\n")

        f.write("\n")

        # Empty sections (required by PSF format)
        f.write("       0 !NDON: donors\n\n")
        f.write("       0 !NACC: acceptors\n\n")
        f.write("       0 !NNB: non-bonded exclusions\n\n")
        f.write("       0       0 !NGRP: groups\n\n")


def print_topology_summary(topology: Topology) -> None:
    """Print summary of topology contents"""
    print("\n" + "="*60)
    print("TOPOLOGY SUMMARY")
    print("="*60)
    print(f"System: {topology.system_name}")

    if topology.atom_type_masses:
        print(f"\nAtom types loaded: {len(topology.atom_type_masses)}")
        # Print unique mass values
        unique_masses = sorted(set(topology.atom_type_masses.values()))
        print(f"  Unique masses: {', '.join(f'{m:.1f}' for m in unique_masses)} amu")

    print(f"\nMolecule types defined: {len(topology.molecule_types)}")
    for name, moltype in topology.molecule_types.items():
        print(f"  {name}: {len(moltype.atoms)} atoms, "
              f"{len(moltype.bonds)} bonds, "
              f"{len(moltype.angles)} angles, "
              f"{len(moltype.dihedrals)} dihedrals, "
              f"{len(moltype.impropers)} impropers")

    print(f"\nMolecules in system:")
    for mol_name, count in topology.molecules:
        print(f"  {count} x {mol_name}")

    print(f"\nInstantiated system:")
    print(f"  {len(topology.atoms):8d} atoms")
    print(f"  {len(topology.bonds):8d} bonds")
    print(f"  {len(topology.angles):8d} angles")
    print(f"  {len(topology.dihedrals):8d} dihedrals")
    print(f"  {len(topology.impropers):8d} impropers")
    print("="*60)


# =============================================================================
# Main Program
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Convert GROMACS topology files to PSF format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -s topol.top -o system.psf
  %(prog)s -s topol.top -p additional.itp -o system.psf

Notes:
  The .top file typically includes .itp files automatically via #include
  directives. Additional .itp files can be specified with -p if needed.

  This tool is designed with modularity for future VMD molfile plugin
  development in C.
        """
    )

    parser.add_argument('-s', '--top', required=True,
                       help='Input GROMACS .top file')
    parser.add_argument('-p', '--itp', nargs='*',
                       help='Additional .itp files (optional)')
    parser.add_argument('-f', '--forcefield',
                       help='Force field .itp file for atom type masses (e.g., martini_v3.0.0.itp)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output PSF file')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    # Create topology object
    topology = Topology()

    # Parse force field file for atom type masses FIRST if provided
    if args.forcefield:
        print(f"Loading atom type masses from force field: {args.forcefield}")
        parse_itp_file(Path(args.forcefield), topology)
        print(f"  Loaded {len(topology.atom_type_masses)} atom types\n")

    # Parse main topology file (pass existing topology to preserve atom types)
    topology = parse_top_file(Path(args.top), topology)

    # Parse additional ITP files if provided
    if args.itp:
        for itp_file in args.itp:
            parse_itp_file(Path(itp_file), topology)

    # Build the complete system
    build_system(topology)

    # Print summary
    print_topology_summary(topology)

    # Write PSF output
    write_psf(topology, Path(args.output))

    print("\nConversion complete!")
    return 0


if __name__ == '__main__':
    exit(main())
