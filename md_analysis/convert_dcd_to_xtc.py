import mdtraj as md
import argparse
import sys
import os

def convert_dcd_to_xtc(dcd_file, psf_file, output_file):
    # Check file extensions
    if not dcd_file.endswith('.dcd'):
        raise ValueError(f"DCD file should have .dcd extension: {dcd_file}")
    if not psf_file.endswith('.psf'):
        raise ValueError(f"PSF file should have .psf extension: {psf_file}")
    
    # Load the dcd file with the topology
    traj = md.load(dcd_file, top=psf_file)
    
    # Write the xtc file
    traj.save(output_file)
    print(f"Converted {dcd_file} to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert DCD trajectory to XTC format')
    parser.add_argument('input_file1', help='Input file (either DCD trajectory or PSF topology)')
    parser.add_argument('input_file2', help='Input file (either DCD trajectory or PSF topology)')
    parser.add_argument('-o', '--output', dest='output_file', help='Output XTC file name')
    
    args = parser.parse_args()
    
    # Determine which file is which based on extension
    if args.input_file1.endswith('.dcd') and args.input_file2.endswith('.psf'):
        dcd_file = args.input_file1
        psf_file = args.input_file2
    elif args.input_file1.endswith('.psf') and args.input_file2.endswith('.dcd'):
        psf_file = args.input_file1
        dcd_file = args.input_file2
    else:
        print("Error: Please provide one .dcd file and one .psf file")
        print(f"Got: {args.input_file1} and {args.input_file2}")
        sys.exit(1)
    
    # If output file is not specified, use the dcd filename but with .xtc extension
    if args.output_file is None:
        args.output_file = os.path.splitext(dcd_file)[0] + '.xtc'
    
    convert_dcd_to_xtc(dcd_file, psf_file, args.output_file)