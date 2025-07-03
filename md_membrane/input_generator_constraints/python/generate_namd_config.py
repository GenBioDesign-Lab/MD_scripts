#!/usr/bin/env python3

import os
import sys
import yaml
import jinja2
import argparse
from pathlib import Path
import glob


def load_config(config_file):
    """Load the YAML configuration file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def bool_to_yesno(value):
    """Convert boolean values to yes/no strings."""
    if isinstance(value, bool):
        return "yes" if value else "no"
    return value


def render_template(template_file, config, stage, temperature=None, stage_dirs=None, use_restraints=False, eq_stage_num=None):
    """Render the Jinja2 template with the configuration for the specified stage and temperature."""
    template_dir = os.path.dirname(os.path.abspath(template_file))
    template_name = os.path.basename(template_file)
    
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(template_dir),
        trim_blocks=True,
        lstrip_blocks=True
    )
    # Add custom filter for boolean to yes/no conversion
    env.filters['yesno'] = bool_to_yesno
    template = env.get_template(template_name)
    
    # Handle multi-stage equilibration first
    if stage.startswith('eq') and eq_stage_num is not None:
        # For equilibration stages, start with common config
        merged_config = {**config['common']}
        
        eq_config = config.get('equilibration_stages', {})
        base_config = eq_config.get('base_config', {})
        stage_specific = eq_config.get('stage_specific', {}).get(stage, {})
        
        # Merge base config with stage-specific overrides
        merged_config.update(base_config)
        merged_config.update(stage_specific)
         
        # Set stage name for template
        merged_config['stage'] = stage
        
        # Map Langevin parameters from common config
        if 'langevin_file' in config['common']:
            merged_config['langevin_file'] = config['common']['langevin_file']
        if 'langevin_col' in config['common']:
            merged_config['langevin_col'] = config['common']['langevin_col']
        if 'langevin_damping' in config['common']:
            merged_config['langevin_damping'] = config['common']['langevin_damping']
    else:
        # For standard stages (minimization, equilibration, production)
        merged_config = {**config['common'], **config[stage]}
    
    if temperature is not None:
        merged_config['temperature'] = temperature
    
    # Add restraints flag to the template context
    merged_config['use_restraints'] = use_restraints
    
    # Add constraints if they are enabled
    constraints_config = config.get('constraints', {})
    if constraints_config.get('enabled', False):
        merged_config['constraints'] = constraints_config
    
    # Update file paths for cross-directory references
    if stage_dirs:
        if stage == 'equilibration':
            rel_path = os.path.relpath(stage_dirs['minimization'], stage_dirs['equilibration'])
            merged_config['coordinates_file'] = os.path.join(rel_path, "final_min.coor")
            merged_config['extended_system_file'] = os.path.join(rel_path, "final_min.xsc")
        elif stage == 'production':
            # For multi-stage equilibration, use last equilibration stage
            if 'equilibration_stages' in config:
                num_stages = config['equilibration_stages'].get('num_stages', 1)
                last_eq_stage = f'eq{num_stages}'
                if last_eq_stage in stage_dirs:
                    rel_path = os.path.relpath(stage_dirs[last_eq_stage], stage_dirs['production'])
                    merged_config['coordinates_file'] = os.path.join(rel_path, f"final_eq{num_stages}.coor")
                    merged_config['extended_system_file'] = os.path.join(rel_path, f"final_eq{num_stages}.xsc")
                    merged_config['velocities_file'] = os.path.join(rel_path, f"final_eq{num_stages}.vel")
                else:
                    # Fallback to single equilibration
                    rel_path = os.path.relpath(stage_dirs['equilibration'], stage_dirs['production'])
                    merged_config['coordinates_file'] = os.path.join(rel_path, "final_eq.coor")
                    merged_config['extended_system_file'] = os.path.join(rel_path, "final_eq.xsc")
                    merged_config['velocities_file'] = os.path.join(rel_path, "final_eq.vel")
            else:
                rel_path = os.path.relpath(stage_dirs['equilibration'], stage_dirs['production'])
                merged_config['coordinates_file'] = os.path.join(rel_path, "final_eq.coor")
                merged_config['extended_system_file'] = os.path.join(rel_path, "final_eq.xsc")
                merged_config['velocities_file'] = os.path.join(rel_path, "final_eq.vel")
        elif stage.startswith('eq') and eq_stage_num is not None:
            if eq_stage_num == 1:
                # First equilibration stage uses minimization output
                rel_path = os.path.relpath(stage_dirs['minimization'], stage_dirs[stage])
                merged_config['coordinates_file'] = os.path.join(rel_path, "final_min.coor")
                merged_config['extended_system_file'] = os.path.join(rel_path, "final_min.xsc")
            else:
                # Subsequent stages use previous equilibration output
                prev_stage = f'eq{eq_stage_num - 1}'
                if prev_stage in stage_dirs:
                    rel_path = os.path.relpath(stage_dirs[prev_stage], stage_dirs[stage])
                    merged_config['coordinates_file'] = os.path.join(rel_path, f"final_eq{eq_stage_num - 1}.coor")
                    merged_config['extended_system_file'] = os.path.join(rel_path, f"final_eq{eq_stage_num - 1}.xsc")
                    merged_config['velocities_file'] = os.path.join(rel_path, f"final_eq{eq_stage_num - 1}.vel")
            
            # Set output name for this equilibration stage
            merged_config['output_name'] = f"final_eq{eq_stage_num}"
    
    return template.render(**merged_config)


def save_file(content, output_path):
    """Save rendered content to file."""
    with open(output_path, 'w') as f:
        f.write(content)
    print(f"Created NAMD configuration file: {output_path}")


def discover_files():
    """Discover .yaml and .j2 files in the script's directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    yaml_files = glob.glob(os.path.join(script_dir, "*.yaml")) + glob.glob(os.path.join(script_dir, "*.yml"))
    j2_files = glob.glob(os.path.join(script_dir, "*.j2"))
    return yaml_files[0] if yaml_files else None, j2_files[0] if j2_files else None


def main():
    default_config, default_template = discover_files()
    
    parser = argparse.ArgumentParser(description='Generate NAMD configuration files from templates')
    parser.add_argument('--config', '-c', default=default_config,
                        help=f'Path to YAML configuration file (default: auto-discover)')
    parser.add_argument('--template', '-t', default=default_template,
                        help=f'Path to Jinja2 template file (default: auto-discover)')
    parser.add_argument('--output-dir', '-o', default='namd',
                        help='Directory to save generated NAMD files (default: namd)')
    parser.add_argument('--stages', '-s', nargs='+', 
                        choices=['minimization', 'equilibration', 'multi-equilibration', 'production', 'all'],
                        default=['all'], help='Stages to generate (default: all)')
    parser.add_argument('--temperatures', '-T', nargs='+', type=float,
                        help='List of temperatures to generate files for')
    parser.add_argument('--temp-prefix', '-p', action='store_true',
                        help='Add temperature as prefix to output directories (T{temp}_)')
    parser.add_argument('--restraints', '-r', action='store_true',
                        help='Include restraints in minimization and equilibration stages')
    parser.add_argument('--eq-stages', type=int,
                        help='Number of equilibration stages (overrides config file)')
    
    args = parser.parse_args()
    
    if args.config is None or args.template is None:
        print("Error: Required .yaml/.yml config or .j2 template file not found", file=sys.stderr)
        sys.exit(1)
    
    try:
        config = load_config(args.config)
    except Exception as e:
        print(f"Error loading configuration file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Override equilibration stages if specified
    if args.eq_stages:
        if 'equilibration_stages' not in config:
            config['equilibration_stages'] = {'base_config': {}, 'force_constants': {}}
        config['equilibration_stages']['num_stages'] = args.eq_stages
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Determine stages to generate
    if 'all' in args.stages:
        if 'equilibration_stages' in config and args.restraints:
            stages_to_generate = ['minimization', 'multi-equilibration', 'production']
        else:
            stages_to_generate = ['minimization', 'equilibration', 'production']
    else:
        stages_to_generate = args.stages
    
    temperatures = args.temperatures if args.temperatures else [None]
    
    stage_filenames = {
        'minimization': '01_Minimization.namd',
        'equilibration': '02_Equilibration.namd',
        'production': '03_Production_npt.namd'
    }
    
    for temp in temperatures:
        current_output_dir = output_dir
        if temp is not None and (len(temperatures) > 1 or args.temp_prefix):
            current_output_dir = output_dir / f"T{int(temp)}"
            current_output_dir.mkdir(exist_ok=True)
        
        # Create stage directories
        stage_dirs = {}
        
        # Handle multi-stage equilibration
        if 'multi-equilibration' in stages_to_generate:
            eq_config = config.get('equilibration_stages', {})
            num_eq_stages = eq_config.get('num_stages', 1)
            
            # Create directories for all stages
            for stage in ['minimization', 'production']:
                if stage in ['minimization'] or stage in stages_to_generate:
                    stage_dir = current_output_dir / stage
                    stage_dir.mkdir(exist_ok=True, parents=True)
                    stage_dirs[stage] = stage_dir
            
            # Create equilibration stage directories
            for i in range(1, num_eq_stages + 1):
                eq_stage = f'eq{i}'
                stage_dir = current_output_dir / eq_stage
                stage_dir.mkdir(exist_ok=True, parents=True)
                stage_dirs[eq_stage] = stage_dir
                stage_filenames[eq_stage] = f'02_{i:02d}_Equilibration_eq{i}.namd'
            
            # Generate equilibration files
            for i in range(1, num_eq_stages + 1):
                eq_stage = f'eq{i}'
                try:
                    content = render_template(args.template, config, eq_stage, temp, stage_dirs, args.restraints, i)
                    output_path = stage_dirs[eq_stage] / stage_filenames[eq_stage]
                    save_file(content, output_path)
                except Exception as e:
                    print(f"Error generating {eq_stage} configuration for temperature {temp}: {e}", file=sys.stderr)
        else:
            # Standard single-stage approach
            for stage in stages_to_generate:
                if stage not in ['multi-equilibration']:
                    stage_dir = current_output_dir / stage
                    stage_dir.mkdir(exist_ok=True, parents=True)
                    stage_dirs[stage] = stage_dir
        
        # Generate configuration files for standard stages
        for stage in stages_to_generate:
            if stage not in ['multi-equilibration'] and stage in stage_filenames:
                try:
                    content = render_template(args.template, config, stage, temp, stage_dirs, args.restraints)
                    output_path = stage_dirs[stage] / stage_filenames[stage]
                    save_file(content, output_path)
                except Exception as e:
                    print(f"Error generating {stage} configuration for temperature {temp}: {e}", file=sys.stderr)
                    
if __name__ == "__main__":
    main() 