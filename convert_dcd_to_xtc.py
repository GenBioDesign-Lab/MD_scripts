import mdtraj as md
import argparse
import sys

def convert_dcd_to_xtc(dcd_file, psf_file, output_file):
    # Load the dcd file with the topology
    traj = md.load(dcd_file, top=psf_file)
    
    # Write the xtc file
    traj.save(output_file)
    print(f"Converted {dcd_file} to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert DCD trajectory to XTC format')
    parser.add_argument('dcd_file', help='Input DCD trajectory file')
    parser.add_argument('psf_file', help='PSF topology file')
    parser.add_argument('-o', '--output', dest='output_file', help='Output XTC file name')
    
    args = parser.parse_args()
    
    # If output file is not specified, use the dcd filename but with .xtc extension
    if args.output_file is None:
        args.output_file = args.dcd_file.rsplit('.', 1)[0] + '.xtc'
    
    convert_dcd_to_xtc(args.dcd_file, args.psf_file, args.output_file)