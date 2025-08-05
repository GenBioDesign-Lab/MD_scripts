# TCL Script Templates

This directory contains template files for TCL scripts used in the membrane preparation workflow. These templates use placeholder substitution to generate customized scripts for each run.

## Template Files

### `build_psf.tcl`
Template for PSF building using CHARMM topology files.

**Placeholders:**
- `{log_file}` - PSF generation log file path
- `{topology_commands}` - Multiple topology file commands
- `{cnt_pdb}` - CNT PDB file path
- `{membrane_pdb}` - Membrane PDB file path  
- `{output_psf}` - Output PSF file path
- `{output_pdb}` - Output PDB file path

### `solvate_ionize.tcl`
Template for system solvation and ionization with integration to modular analysis scripts.

**Placeholders:**
- `{xsc_file}` - XSC file for box dimensions (optional, can be empty string)
- `{system_psf}` - Input PSF file path
- `{system_pdb}` - Input PDB file path
- `{salt_concentration}` - Salt concentration for ionization

**Features:**
- Integrates with `analyze_atoms.tcl` and `calculate_padding.tcl` modules
- Handles XSC file reading for custom box dimensions
- Creates temporary and final output directories
- Performs solvation with calculated padding
- Ionizes system with specified salt concentration

### `create_constraints.tcl`
Template for creating constraint files for CNT atoms.

**Placeholders:**
- `{system_psf}` - System PSF file path
- `{system_pdb}` - System PDB file path
- `{constraint_force}` - Constraint force constant (kcal/mol/Å²)
- `{output_pdb}` - Output constraint PDB file path

### `create_langevin.tcl`
Template for creating Langevin damping files for water molecules.

**Placeholders:**
- `{system_psf}` - System PSF file path
- `{system_pdb}` - System PDB file path
- `{output_pdb}` - Output Langevin PDB file path

## Usage

Templates are automatically loaded and processed by the `cnt_membrane_preparation.py` script using the `load_template()` method. Placeholders are replaced using Python's string formatting.

## Editing Templates

You can modify these templates to:
- Change VMD/TCL commands
- Add new functionality
- Adjust selection criteria
- Modify output formats
- Change solvation parameters
- Adjust constraint selections

When editing, preserve the `{placeholder}` syntax for variables that need to be substituted by the Python script.

## Integration with Modular Scripts

The `solvate_ionize.tcl` template works with the modular solvation system:
- Copies `analyze_atoms.tcl` and `calculate_padding.tcl` to the working directory
- Sources these modules for system analysis and padding calculation
- Uses calculated padding values for solvation box sizing

## Placeholder Format

- Use curly braces: `{variable_name}`
- Variable names should match the parameters passed from the Python script
- No spaces inside braces: `{good}` not `{ bad }`
- Case-sensitive: `{File}` ≠ `{file}`
- Empty values: Use empty string `""` for optional placeholders like `{xsc_file}` 