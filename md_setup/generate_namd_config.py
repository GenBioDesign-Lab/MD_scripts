import os
import sys
import yaml
import jinja2
import argparse
from pathlib import Path


def load_config(config_file):
    """Load the YAML configuration file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def render_template(template_file, config, stage, temperature=None):
    """Render the Jinja2 template with the configuration for the specified stage and temperature."""
    # Load the template from file
    template_dir = os.path.dirname(os.path.abspath(template_file))
    template_name = os.path.basename(template_file)
    
    # Setup Jinja2 environment
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(template_dir),
        trim_blocks=True,
        lstrip_blocks=True
    )
    template = env.get_template(template_name)
    
    # Merge common config with stage-specific config
    merged_config = {**config['common'], **config[stage]}
    
    # Override temperature if specified
    if temperature is not None:
        merged_config['temperature'] = temperature
    
    # Render template
    return template.render(**merged_config)


def save_file(content, output_path):
    """Save rendered content to file."""
    with open(output_path, 'w') as f:
        f.write(content)
    print(f"Created NAMD configuration file: {output_path}")


def get_psf_basename(psf_path):
    """Extract the base name of the PSF file without extension."""
    return os.path.splitext(os.path.basename(psf_path))[0]


def main():
    parser = argparse.ArgumentParser(description='Generate NAMD configuration files from templates')
    parser.add_argument('--config', '-c', default='namd_config.yaml',
                        help='Path to YAML configuration file (default: namd_config.yaml)')
    parser.add_argument('--template', '-t', default='namd_template.j2',
                        help='Path to Jinja2 template file (default: namd_template.j2)')
    parser.add_argument('--output-dir', '-o', default='.',
                        help='Directory to save generated NAMD files (default: current directory)')
    parser.add_argument('--stages', '-s', nargs='+', 
                        choices=['minimization', 'equilibration', 'production', 'all'],
                        default=['all'],
                        help='Stages to generate (default: all)')
    parser.add_argument('--temperatures', '-T', nargs='+', type=float,
                        help='List of temperatures to generate files for (default: use temperature from config)')
    parser.add_argument('--temp-prefix', '-p', action='store_true',
                        help='Add temperature as prefix to output directories (T{temp}_)')
    parser.add_argument('--no-organize', '-n', action='store_true',
                        help='Disable automatic organization into subfolder based on PSF name')
    
    args = parser.parse_args()
    
    # Load configuration
    try:
        config = load_config(args.config)
    except Exception as e:
        print(f"Error loading configuration file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Get PSF basename for subfolder (default behavior)
    psf_basename = None
    if 'structure_file' in config['common'] and not args.no_organize:
        psf_basename = get_psf_basename(config['common']['structure_file'])
    
    # Determine which stages to generate
    stages_to_generate = []
    if 'all' in args.stages:
        stages_to_generate = ['minimization', 'equilibration', 'production']
    else:
        stages_to_generate = args.stages
    
    # Get temperatures to use
    temperatures = [None]  # Default to using config temperature
    if args.temperatures:
        temperatures = args.temperatures
    
    # Define stage filenames
    stage_filenames = {
        'minimization': '01_Minimization.namd',
        'equilibration': '02_Equilibration.namd',
        'production': '03_Production_npt.namd'
    }
    
    for temp in temperatures:
        # Create temperature-specific directory if multiple temperatures
        if temp is not None and (len(temperatures) > 1 or args.temp_prefix):
            temp_dir = output_dir / f"T{int(temp)}"
            temp_dir.mkdir(exist_ok=True)
            current_output_dir = temp_dir
        else:
            current_output_dir = output_dir
        
        # Create PSF-named directory (default behavior)
        if psf_basename:
            psf_dir = current_output_dir / psf_basename
            psf_dir.mkdir(exist_ok=True)
            current_output_dir = psf_dir
        
        for stage in stages_to_generate:
            try:
                # Render template for this stage and temperature
                content = render_template(args.template, config, stage, temp)
                
                # Determine output path
                output_path = current_output_dir / stage_filenames[stage]
                
                # Save to output file
                save_file(content, output_path)
            except Exception as e:
                print(f"Error generating {stage} configuration for temperature {temp}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main() 