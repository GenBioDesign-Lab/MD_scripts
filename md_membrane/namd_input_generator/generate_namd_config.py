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

def render_template(template_file, config, stage, temperature=None, stage_dirs=None):
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
    
    merged_config = {**config['common'], **config[stage]}
    
    if temperature is not None:
        merged_config['temperature'] = temperature
    
    # Update file paths for cross-directory references
    if stage_dirs:
        if stage == 'equilibration':
            rel_path = os.path.relpath(stage_dirs['minimization'], stage_dirs['equilibration'])
            merged_config['coordinates_file'] = os.path.join(rel_path, "final_min.coor")
            merged_config['extended_system_file'] = os.path.join(rel_path, "final_min.xsc")
        elif stage == 'production':
            rel_path = os.path.relpath(stage_dirs['equilibration'], stage_dirs['production'])
            merged_config['coordinates_file'] = os.path.join(rel_path, "final_eq.coor")
            merged_config['extended_system_file'] = os.path.join(rel_path, "final_eq.xsc")
            merged_config['velocities_file'] = os.path.join(rel_path, "final_eq.vel")
    
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
                        choices=['minimization', 'equilibration', 'production', 'all'],
                        default=['all'], help='Stages to generate (default: all)')
    parser.add_argument('--temperatures', '-T', nargs='+', type=float,
                        help='List of temperatures to generate files for')
    parser.add_argument('--temp-prefix', '-p', action='store_true',
                        help='Add temperature as prefix to output directories (T{temp}_)')
    
    args = parser.parse_args()
    
    if args.config is None or args.template is None:
        print("Error: Required .yaml/.yml config or .j2 template file not found", file=sys.stderr)
        sys.exit(1)
    
    try:
        config = load_config(args.config)
    except Exception as e:
        print(f"Error loading configuration file: {e}", file=sys.stderr)
        sys.exit(1)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    stages_to_generate = ['minimization', 'equilibration', 'production'] if 'all' in args.stages else args.stages
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
        for stage in stages_to_generate:
            stage_dir = current_output_dir / stage
            stage_dir.mkdir(exist_ok=True, parents=True)
            stage_dirs[stage] = stage_dir
        
        # Generate configuration files
        for stage in stages_to_generate:
            try:
                content = render_template(args.template, config, stage, temp, stage_dirs)
                output_path = stage_dirs[stage] / stage_filenames[stage]
                save_file(content, output_path)
            except Exception as e:
                print(f"Error generating {stage} configuration for temperature {temp}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main() 