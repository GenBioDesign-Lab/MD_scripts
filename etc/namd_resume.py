#!/usr/bin/env python3
"""
Usage:
    python namd_resume.py -i <previous_config_file> -n <new_numstep_value>
"""

import os
import sys
import re
import argparse
from pathlib import Path

def parse_namd_config(config_file):
    """Parse NAMD configuration file and extract key parameters."""
    
    config = {}
    
    with open(config_file, 'r') as f:
        content = f.read()
    
    # Define patterns for key parameters
    patterns = {
        'structure': r'structure\s+(.+)',
        'coordinates': r'coordinates\s+(.+)',
        'extendedSystem': r'extendedSystem\s+(.+)',
        'velocities': r'velocities\s+(.+)',
        'outputname': r'set\s+outputname\s+(.+)',
        'firsttimestep': r'firsttimestep\s+(\d+)',
        'numsteps': r'numsteps\s+(\d+)',
        'temperature': r'set\s+temperature\s+(\d+)',
        'parameters': r'parameters\s+(.+)',
    }
    
    # Extract single-line parameters
    for key, pattern in patterns.items():
        match = re.search(pattern, content, re.IGNORECASE)
        if match:
            config[key] = match.group(1).strip()
    
    # Extract all parameter files
    param_matches = re.findall(r'parameters\s+(.+)', content, re.IGNORECASE)
    if param_matches:
        config['parameters_files'] = [p.strip() for p in param_matches]
    
    # Get the directory of the config file for relative path resolution
    config['config_dir'] = os.path.dirname(os.path.abspath(config_file))
    
    return config

def determine_output_paths(config):
    """Determine the output file paths from the previous simulation."""
    
    if 'outputname' not in config:
        raise ValueError("Could not find outputname in the configuration file")
    
    outputname = config['outputname']
    config_dir = config['config_dir']
    
    # Construct output file paths
    output_files = {
        'coor': f"{outputname}.coor",
        'xsc': f"{outputname}.xsc", 
        'vel': f"{outputname}.vel",
        'restart_coor': f"{outputname}.restart.coor",
        'restart_xsc': f"{outputname}.restart.xsc",
        'restart_vel': f"{outputname}.restart.vel"
    }
    
    # Check which files exist (prefer .restart files if available)
    existing_files = {}
    for file_type, filename in output_files.items():
        filepath = os.path.join(config_dir, filename)
        if os.path.exists(filepath):
            existing_files[file_type] = filename
    
    # Determine best files to use
    coordinates_file = existing_files.get('restart_coor', existing_files.get('coor'))
    xsc_file = existing_files.get('restart_xsc', existing_files.get('xsc'))
    vel_file = existing_files.get('restart_vel', existing_files.get('vel'))
    
    if not coordinates_file:
        raise FileNotFoundError(f"Could not find output coordinate file for {outputname}")
    if not xsc_file:
        raise FileNotFoundError(f"Could not find output extended system file for {outputname}")
    
    return {
        'coordinates': coordinates_file,
        'extendedSystem': xsc_file,
        'velocities': vel_file
    }

def get_final_timestep(config, output_files):
    """Determine the final timestep from the previous simulation."""
    
    # Try to read from .xsc file to get the timestep
    xsc_file = os.path.join(config['config_dir'], output_files['extendedSystem'])
    
    try:
        with open(xsc_file, 'r') as f:
            lines = f.readlines()
            # The timestep is usually in the second line of .xsc file
            if len(lines) >= 2:
                timestep_match = re.search(r'(\d+)', lines[1])
                if timestep_match:
                    return int(timestep_match.group(1))
    except:
        pass
    
    # Fallback: calculate from original firsttimestep + numsteps
    first_timestep = int(config.get('firsttimestep', 0))
    num_steps = int(config.get('numsteps', 0))
    
    return first_timestep + num_steps

def generate_resume_config(config_file, config, output_files, new_outputname, new_numsteps, new_firsttimestep):
    """Generate the resume configuration content."""
    
    # Read original config file
    with open(config_file, 'r') as f:
        content = f.read()
    
    # Update coordinates
    if 'coordinates' in config:
        content = re.sub(
            r'coordinates\s+.+',
            f"coordinates        {output_files['coordinates']}",
            content,
            flags=re.IGNORECASE
        )
    
    # Update extendedSystem
    if 'extendedSystem' in config:
        content = re.sub(
            r'extendedSystem\s+.+',
            f"extendedSystem     {output_files['extendedSystem']}",
            content,
            flags=re.IGNORECASE
        )
    elif output_files['extendedSystem']:
        # Add extendedSystem if it wasn't in the original (for minimization -> equilibration)
        coordinates_line = re.search(r'coordinates\s+.+', content, re.IGNORECASE)
        if coordinates_line:
            content = content.replace(
                coordinates_line.group(0),
                coordinates_line.group(0) + f"\nextendedSystem     {output_files['extendedSystem']}"
            )
    
    # Update or add velocities
    if output_files['velocities']:
        if 'velocities' in config:
            content = re.sub(
                r'velocities\s+.+',
                f"velocities         {output_files['velocities']}",
                content,
                flags=re.IGNORECASE
            )
        else:
            # Add velocities line after extendedSystem
            extended_line = re.search(r'extendedSystem\s+.+', content, re.IGNORECASE)
            if extended_line:
                content = content.replace(
                    extended_line.group(0),
                    extended_line.group(0) + f"\nvelocities         {output_files['velocities']}"
                )
    
    # Update outputname
    content = re.sub(
        r'set\s+outputname\s+.+',
        f"set outputname     {new_outputname}",
        content,
        flags=re.IGNORECASE
    )
    
    # Update firsttimestep
    content = re.sub(
        r'firsttimestep\s+\d+',
        f"firsttimestep      {new_firsttimestep}",
        content,
        flags=re.IGNORECASE
    )
    
    # Update numsteps
    content = re.sub(
        r'numsteps\s+\d+',
        f"numsteps           {new_numsteps}",
        content,
        flags=re.IGNORECASE
    )
    
    return content

def main():
    parser = argparse.ArgumentParser(description='Generate NAMD resume configuration')
    parser.add_argument('-i', '--input', required=True, help='Previous NAMD config file')
    parser.add_argument('-n', '--numsteps', required=True, type=int, help='Number of steps for resume simulation')
    parser.add_argument('-o', '--output', help='Output config name (default: {inputname}_resume)')
    
    args = parser.parse_args()
    
    config_file = args.input
    new_numsteps = args.numsteps
    
    if not os.path.exists(config_file):
        print(f"Error: Configuration file '{config_file}' not found")
        sys.exit(1)
    
    try:
        # Parse the previous configuration
        config = parse_namd_config(config_file)
        
        # Determine output file paths
        output_files = determine_output_paths(config)
        
        # Calculate the starting timestep for resume
        final_timestep = get_final_timestep(config, output_files)
        
        # Determine output name
        if args.output:
            new_outputname = args.output
        else:
            base_name = os.path.splitext(os.path.basename(config_file))[0]
            new_outputname = f"{base_name}_resume"
        
        # Generate the resume configuration
        resume_content = generate_resume_config(
            config_file, config, output_files, new_outputname, new_numsteps, final_timestep
        )
        
        # Write the new configuration file
        output_config = f"{new_outputname}.namd"
        with open(output_config, 'w') as f:
            f.write(resume_content)
        
        print(f"Resume config: {output_config}")
        print(f"Start timestep: {final_timestep}")
        print(f"Run steps: {new_numsteps}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 