#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import shutil
import json
from datetime import datetime
from pathlib import Path

class MembranePreparationWorkflow:
    """
    Workflow for membrane preparation from raw CNT and membrane PDB files
    to final ionized system ready for MD simulation.
    """
    
    def __init__(self, config):
        self.config = config
        self.verbose = config.get('verbose', False)
        self.setup_directories()
        self.setup_logging()
        
    def setup_directories(self):
        """Create organized output directory structure"""
        self.base_dir = Path(self.config['output_dir'])
        
        # Directory structure
        self.dirs = {
            'base': self.base_dir,
            'logs': self.base_dir / 'logs',
            'step1_centering': self.base_dir / '01_centering',
            'step2_lipid_removal': self.base_dir / '02_lipid_removal', 
            'step3_psf_building': self.base_dir / '03_psf_building',
            'step4_solvation': self.base_dir / '04_solvation',
            'step5_constraints': self.base_dir / '05_constraints',
            'final_system': self.base_dir / 'final_system',
            'intermediate': self.base_dir / 'intermediate_files'
        }
        
        # Create all directories
        for dir_path in self.dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)
            
        if self.verbose:
            print(f"Created output directory structure in: {self.base_dir}")
        
    def setup_logging(self):
        """Setup comprehensive logging"""
        self.log_file = self.dirs['logs'] / f"membrane_prep_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        print(f"Detailed logs will be saved to: {self.log_file}")
    
    def load_template(self, template_name, **kwargs):
        """Load a TCL template file and substitute placeholders"""
        template_dir = Path(__file__).parent / 'one_click_templates'
        template_file = template_dir / template_name
        
        if not template_file.exists():
            self.log_and_print(f"Error: Template file not found: {template_file}", important=True)
            self.log_and_print(f"Template directory: {template_dir}", console=False)
            self.log_and_print(f"Available templates: {list(template_dir.glob('*.tcl')) if template_dir.exists() else 'Directory does not exist'}", console=False)
            raise FileNotFoundError(f"Template file not found: {template_file}")
        
        try:
            with open(template_file, 'r') as f:
                template_content = f.read()
            
            # Replace placeholders with actual values
            return template_content.format(**kwargs)
        except KeyError as e:
            self.log_and_print(f"Error: Missing placeholder in template {template_name}: {e}", important=True)
            self.log_and_print(f"Available placeholders: {list(kwargs.keys())}", console=False)
            raise
        
    def log_and_print(self, message, console=True, important=False):
        """Log message and optionally print to console"""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_message = f"[{timestamp}] {message}"
        
        # Always log to file
        with open(self.log_file, 'a') as f:
            f.write(log_message + '\n')
        
        # Print to console based on verbosity settings
        if console and (self.verbose or important):
            print(log_message)
        elif important:
            print(message)  # Print without timestamp for important messages
    
    def run_command(self, cmd, step_name, cwd=None):
        """Run a command and handle errors"""
        self.log_and_print(f"Running {step_name}...", important=True)
        self.log_and_print(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}", console=False)
        
        try:
            if cwd:
                result = subprocess.run(cmd, shell=isinstance(cmd, str), cwd=cwd, 
                                      capture_output=True, text=True, check=True)
            else:
                result = subprocess.run(cmd, shell=isinstance(cmd, str), 
                                      capture_output=True, text=True, check=True)
            
            if result.stdout:
                self.log_and_print(f"STDOUT:\n{result.stdout}", console=False)
            if result.stderr:
                self.log_and_print(f"STDERR:\n{result.stderr}", console=False)
                
            self.log_and_print(f"{step_name} completed successfully", important=True)
            return result
            
        except subprocess.CalledProcessError as e:
            self.log_and_print(f"✗ {step_name} failed!", important=True)
            self.log_and_print(f"Error code: {e.returncode}", console=False)
            self.log_and_print(f"STDOUT: {e.stdout}", console=False)
            self.log_and_print(f"STDERR: {e.stderr}", console=False)
            print(f"✗ {step_name} failed! Check log file for details: {self.log_file}")
            raise
    
    def step1_center_cnt_in_membrane(self):
        """Step 1: Center CNT in membrane"""
        self.log_and_print("="*60, console=False)
        self.log_and_print("STEP 1: CENTERING CNT IN MEMBRANE", important=True)
        self.log_and_print("="*60, console=False)
        
        # Input files
        cnt_file = self.config['cnt_pdb']
        membrane_file = self.config['membrane_pdb']
        
        # Output file
        centered_cnt = self.dirs['step1_centering'] / 'cnt_centered.pdb'
        
        # Prepare command
        cmd = [
            'python3', str(Path(__file__).parent / 'center_membrane.py'),
            membrane_file, cnt_file,
            '-o', str(centered_cnt),
            '-r', str(self.config.get('rotation_x', 90)), 
                  str(self.config.get('rotation_y', 0)), 
                  str(self.config.get('rotation_z', 0))
        ]
        
        self.run_command(cmd, "CNT centering")
        
        self.files = {'centered_cnt': centered_cnt}
        return centered_cnt
    
    def step2_remove_overlapping_lipids(self):
        """Step 2: Remove overlapping lipids and create complex"""
        self.log_and_print("="*60, console=False)
        self.log_and_print("STEP 2: REMOVING OVERLAPPING LIPIDS", important=True)
        self.log_and_print("="*60, console=False)
        
        # Input files
        centered_cnt = self.files['centered_cnt']
        membrane_file = self.config['membrane_pdb']
        
        # Output files
        cleaned_membrane = self.dirs['step2_lipid_removal'] / 'cleaned_membrane.pdb'
        complex_structure = self.dirs['step2_lipid_removal'] / 'cnt_membrane_complex.pdb'
        removal_log = self.dirs['step2_lipid_removal'] / 'lipid_removal.log'
        
        # Prepare command
        cmd = [
            'python3', str(Path(__file__).parent / 'remove_lipid.py'),
            str(centered_cnt), membrane_file,
            '-o', str(cleaned_membrane),
            '-c', str(complex_structure),
            '-r', str(self.config.get('radius_factor', 1.0)),
            '-b', str(self.config.get('buffer_distance', 0.8)),
            '-l', str(removal_log)
        ]
        
        if self.config.get('symmetric_removal', True):
            cmd.append('-s')
            
        self.run_command(cmd, "Lipid removal")
        
        self.files.update({
            'cleaned_membrane': cleaned_membrane,
            'complex_structure': complex_structure
        })
        
        return cleaned_membrane, complex_structure
    
    def step3_build_psf(self):
        """Step 3: Build PSF file using CHARMM topology"""
        self.log_and_print("="*60, console=False)
        self.log_and_print("STEP 3: BUILDING PSF FILE", important=True)
        self.log_and_print("="*60, console=False)
        
        # Create custom TCL script for PSF building
        psf_script = self.dirs['step3_psf_building'] / 'build_psf.tcl'
        
        # Prepare file paths
        centered_cnt = self.files['centered_cnt']
        cleaned_membrane = self.files['cleaned_membrane']
        output_psf = self.dirs['step3_psf_building'] / 'cnt_membrane.psf'
        output_pdb = self.dirs['step3_psf_building'] / 'cnt_membrane.pdb'
        
        # Prepare topology commands
        topology_commands = ""
        if 'topology_files' in self.config and self.config['topology_files']:
            for topo_file in self.config['topology_files']:
                topology_commands += f"topology {topo_file}\n"
        else:
            # Default topology files (update paths as needed)
            default_topos = [
                "/data01/genbiolab/mdanh/data/CHARMMFF/cnt_ff/final/multi_wall/cnt_5_45_none.rtf",
                "/data01/genbiolab/mdanh/data/CHARMMFF/cnt_ff/final/multi_wall/cnt_10_45_none.rtf", 
                "/data01/genbiolab/mdanh/data/CHARMMFF/cnt_ff/final/multi_wall/cnt_15_45_oh.rtf",
                "/data01/genbiolab/mdanh/data/CHARMMFF/toppar/top_all36_lipid.rtf"
            ]
            for topo in default_topos:
                if os.path.exists(topo):
                    topology_commands += f"topology {topo}\n"
        
        # Generate TCL script content from template
        tcl_content = self.load_template('build_psf.tcl',
            log_file=self.dirs['step3_psf_building'] / 'psfgen.log',
            topology_commands=topology_commands.rstrip(),
            cnt_pdb=centered_cnt,
            membrane_pdb=cleaned_membrane,
            output_psf=output_psf,
            output_pdb=output_pdb
        )
        
        # Write TCL script
        with open(psf_script, 'w') as f:
            f.write(tcl_content)
            
        # Run VMD with TCL script
        cmd = ['vmd', '-dispdev', 'text', '-e', str(psf_script)]
        self.run_command(cmd, "PSF building", cwd=self.dirs['step3_psf_building'])
        
        self.files.update({
            'system_psf': output_psf,
            'system_pdb': output_pdb
        })
        
        return output_psf, output_pdb
    
    def step4_solvate_and_ionize(self):
        """Step 4: Solvate and ionize the system"""
        self.log_and_print("="*60, console=False)
        self.log_and_print("STEP 4: SOLVATION AND IONIZATION", important=True)
        self.log_and_print("="*60, console=False)
        
        # Copy solvation scripts to step directory
        script_dir = Path(__file__).parent / 'solvate_ionize_membrane_system'
        step_dir = self.dirs['step4_solvation']
        
        for script_file in ['solvate_ion_modular.tcl', 'analyze_atoms.tcl', 'calculate_padding.tcl']:
            shutil.copy2(script_dir / script_file, step_dir / script_file)
        
        # Create modified solvation script based on the modular template
        solvate_script = step_dir / 'solvate_system.tcl'
        
        # Prepare paths
        system_psf = self.files['system_psf']
        system_pdb = self.files['system_pdb']
        xsc_file = self.config.get('xsc_file', None)
        
        # Generate solvation script from template
        tcl_content = self.load_template('solvate_ionize.tcl',
            xsc_file=str(xsc_file) if xsc_file else "",
            system_psf=system_psf,
            system_pdb=system_pdb,
            salt_concentration=self.config.get('salt_concentration', 0.15)
        )
        
        with open(solvate_script, 'w') as f:
            f.write(tcl_content)
            
        # Run solvation
        cmd = ['vmd', '-dispdev', 'text', '-e', str(solvate_script)]
        self.run_command(cmd, "Solvation and ionization", cwd=step_dir)
        
        # Final system files
        final_psf = step_dir / 'final' / 'system_ionized.psf'
        final_pdb = step_dir / 'final' / 'system_ionized.pdb'
        
        # Check if solvation produced the expected output files
        if not final_psf.exists():
            self.log_and_print(f"Error: Solvation did not produce PSF file: {final_psf}", important=True)
            self.log_and_print(f"Check solvation log and VMD output for errors", console=False)
            raise FileNotFoundError(f"Solvation failed to create PSF file: {final_psf}")
        
        if not final_pdb.exists():
            self.log_and_print(f"Error: Solvation did not produce PDB file: {final_pdb}", important=True)
            self.log_and_print(f"Check solvation log and VMD output for errors", console=False)
            raise FileNotFoundError(f"Solvation failed to create PDB file: {final_pdb}")
        
        self.log_and_print(f"Solvation output verified: PSF={final_psf.name}, PDB={final_pdb.name}", console=False)
        
        self.files.update({
            'final_psf': final_psf,
            'final_pdb': final_pdb
        })
        
        return final_psf, final_pdb
    
    def step5_create_constraint_files(self):
        """Step 5: Create constraint and Langevin files"""
        self.log_and_print("="*60, console=False)
        self.log_and_print("STEP 5: CREATING CONSTRAINT FILES", important=True)
        self.log_and_print("="*60, console=False)
        
        final_psf = self.files['final_psf']
        final_pdb = self.files['final_pdb']
        
        # Create constraint file
        constraint_script = self.dirs['step5_constraints'] / 'create_constraints.tcl'
        constraint_pdb = self.dirs['step5_constraints'] / 'cnt_constraints.pdb'
        
        tcl_content = self.load_template('create_constraints.tcl',
            system_psf=final_psf,
            system_pdb=final_pdb,
            constraint_force=self.config.get('constraint_force', 2.5),
            output_pdb=constraint_pdb
        )
        
        with open(constraint_script, 'w') as f:
            f.write(tcl_content)
            
        self.run_command(['vmd', '-dispdev', 'text', '-e', str(constraint_script)], 
                        "Constraint file creation")
        
        # Create Langevin file
        langevin_script = self.dirs['step5_constraints'] / 'create_langevin.tcl'
        langevin_pdb = self.dirs['step5_constraints'] / 'langevin_water.pdb'
        
        tcl_content = self.load_template('create_langevin.tcl',
            system_psf=final_psf,
            system_pdb=final_pdb,
            output_pdb=langevin_pdb
        )
        
        with open(langevin_script, 'w') as f:
            f.write(tcl_content)
            
        self.run_command(['vmd', '-dispdev', 'text', '-e', str(langevin_script)], 
                        "Langevin file creation")
        
        self.files.update({
            'constraint_pdb': constraint_pdb,
            'langevin_pdb': langevin_pdb
        })
    
    def finalize_outputs(self):
        """Copy final outputs to organized final directory"""
        self.log_and_print("="*60, console=False)
        self.log_and_print("FINALIZING OUTPUTS", important=True)
        self.log_and_print("="*60, console=False)
        
        final_dir = self.dirs['final_system']
        
        # Check if final system files exist before copying
        if not self.files['final_psf'].exists():
            self.log_and_print(f"Error: Final PSF file not found: {self.files['final_psf']}", important=True)
            raise FileNotFoundError(f"Final PSF file not found: {self.files['final_psf']}")
        
        if not self.files['final_pdb'].exists():
            self.log_and_print(f"Error: Final PDB file not found: {self.files['final_pdb']}", important=True)
            raise FileNotFoundError(f"Final PDB file not found: {self.files['final_pdb']}")
        
        # Copy final system files
        shutil.copy2(self.files['final_psf'], final_dir / 'system_final.psf')
        shutil.copy2(self.files['final_pdb'], final_dir / 'system_final.pdb')
        
        # Copy constraint files
        if 'constraint_pdb' in self.files:
            shutil.copy2(self.files['constraint_pdb'], final_dir / 'constraints.pdb')
        if 'langevin_pdb' in self.files:
            shutil.copy2(self.files['langevin_pdb'], final_dir / 'langevin.pdb')
        
        # Create summary file
        summary = {
            'workflow_completed': datetime.now().isoformat(),
            'input_files': {
                'cnt_pdb': str(self.config['cnt_pdb']),
                'membrane_pdb': str(self.config['membrane_pdb']),
                'xsc_file': str(self.config.get('xsc_file', 'None'))
            },
            'final_outputs': {
                'system_psf': str(final_dir / 'system_final.psf'),
                'system_pdb': str(final_dir / 'system_final.pdb'),
                'constraints_pdb': str(final_dir / 'constraints.pdb'),
                'langevin_pdb': str(final_dir / 'langevin.pdb')
            },
            'parameters': self.config
        }
        
        with open(final_dir / 'workflow_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
            
        self.log_and_print(f"Final system files saved to: {final_dir}", important=True)
        self.log_and_print(f"Workflow summary saved to: {final_dir / 'workflow_summary.json'}", console=False)
    
    def run_complete_workflow(self):
        """Run the complete membrane preparation workflow"""
        self.log_and_print("="*80, console=False)
        print("MEMBRANE PREPARATION WORKFLOW STARTED")
        self.log_and_print("MEMBRANE PREPARATION WORKFLOW STARTED", console=False)
        self.log_and_print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", console=False)
        self.log_and_print("="*80, console=False)
        
        try:
            # Step 1: Center CNT in membrane
            self.step1_center_cnt_in_membrane()
            
            # Step 2: Remove overlapping lipids
            self.step2_remove_overlapping_lipids()
            
            # Step 3: Build PSF
            self.step3_build_psf()
            
            # Step 4: Solvate and ionize
            self.step4_solvate_and_ionize()
            
            # Step 5: Create constraint files
            if self.config.get('create_constraints', True):
                self.step5_create_constraint_files()
            
            # Finalize outputs
            self.finalize_outputs()
            
            self.log_and_print("="*80, console=False)
            print("WORKFLOW COMPLETED SUCCESSFULLY!")
            self.log_and_print("WORKFLOW COMPLETED SUCCESSFULLY!", console=False)
            print(f"Final system ready for MD simulation")
            self.log_and_print(f"Final system ready for MD simulation", console=False)
            print(f"Output directory: {self.dirs['final_system']}")
            self.log_and_print(f"Output directory: {self.dirs['final_system']}", console=False)
            self.log_and_print("="*80, console=False)
            
        except Exception as e:
            self.log_and_print("="*80, console=False)
            print("✗ WORKFLOW FAILED!")
            self.log_and_print("✗ WORKFLOW FAILED!", console=False)
            print(f"Error: {str(e)}")
            self.log_and_print(f"Error: {str(e)}", console=False)
            print(f"Check detailed log file: {self.log_file}")
            self.log_and_print("="*80, console=False)
            raise

def main():
    parser = argparse.ArgumentParser(description='Complete membrane preparation workflow')
    
    # Required arguments
    parser.add_argument('cnt_pdb', help='Input CNT PDB file')
    parser.add_argument('membrane_pdb', help='Input membrane PDB file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for all results')
    
    # Optional topology files
    parser.add_argument('-t', '--topology_files', nargs='+', 
                       help='CHARMM topology files (RTF format)')
    
    # Optional XSC file for box dimensions
    parser.add_argument('--xsc_file', help='XSC file for target box dimensions')
    
    # CNT rotation parameters
    parser.add_argument('--rotation_x', type=float, default=90, help='X-axis rotation (degrees)')
    parser.add_argument('--rotation_y', type=float, default=0, help='Y-axis rotation (degrees)')  
    parser.add_argument('--rotation_z', type=float, default=0, help='Z-axis rotation (degrees)')
    
    # Lipid removal parameters
    parser.add_argument('--radius_factor', type=float, default=1.0, help='CNT radius scaling factor')
    parser.add_argument('--buffer_distance', type=float, default=0.8, help='Buffer distance (Å)')
    parser.add_argument('--symmetric_removal', action='store_true', default=True, 
                       help='Enable symmetric lipid removal')
    
    # Ionization parameters
    parser.add_argument('--salt_concentration', type=float, default=0.15, 
                       help='Salt concentration for ionization')
    
    # Constraint parameters
    parser.add_argument('--constraint_force', type=float, default=2.5, 
                       help='Constraint force constant (kcal/mol/Å²)')
    parser.add_argument('--no_constraints', action='store_true', 
                       help='Skip constraint file creation')
    
    # Verbosity control
    parser.add_argument('-v', '--verbose', action='store_true', 
                       help='Enable verbose console output (default: minimal output)')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.cnt_pdb):
        print(f"Error: CNT PDB file not found: {args.cnt_pdb}")
        sys.exit(1)
        
    if not os.path.exists(args.membrane_pdb):
        print(f"Error: Membrane PDB file not found: {args.membrane_pdb}")
        sys.exit(1)
    
    if args.xsc_file and not os.path.exists(args.xsc_file):
        print(f"Error: XSC file not found: {args.xsc_file}")
        sys.exit(1)
    
    # Prepare configuration
    config = {
        'cnt_pdb': os.path.abspath(args.cnt_pdb),
        'membrane_pdb': os.path.abspath(args.membrane_pdb),
        'output_dir': os.path.abspath(args.output_dir),
        'xsc_file': os.path.abspath(args.xsc_file) if args.xsc_file else None,
        'topology_files': [os.path.abspath(f) for f in args.topology_files] if args.topology_files else None,
        'rotation_x': args.rotation_x,
        'rotation_y': args.rotation_y,
        'rotation_z': args.rotation_z,
        'radius_factor': args.radius_factor,
        'buffer_distance': args.buffer_distance,
        'symmetric_removal': args.symmetric_removal,
        'salt_concentration': args.salt_concentration,
        'constraint_force': args.constraint_force,
        'create_constraints': not args.no_constraints,
        'verbose': args.verbose
    }
    
    # Run workflow
    workflow = MembranePreparationWorkflow(config)
    workflow.run_complete_workflow()

if __name__ == "__main__":
    main() 