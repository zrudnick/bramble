#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict

def parse_sam_file(filename):
    """Parse SAM file and group lines by read name."""
    reads = defaultdict(list)
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip header lines (start with @)
            if line.startswith('@'):
                continue
            
            # Skip empty lines
            if not line:
                continue
            
            # Extract read name (first column)
            fields = line.split('\t')
            if len(fields) >= 1:
                read_name = fields[0]
                reads[read_name].append(line)
    
    return reads

def main():
    parser = argparse.ArgumentParser(description='Separate SAM file lines by read name')
    parser.add_argument('sam1', help='First SAM file')
    parser.add_argument('sam2', help='Second SAM file')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    
    args = parser.parse_args()
    
    # Parse both SAM files
    print(f"Reading {args.sam1}...", file=sys.stderr)
    sam1_reads = parse_sam_file(args.sam1)
    for k,v in sam1_reads.items():
        sam1_reads[k] = list(sorted(v))
    
    print(f"Reading {args.sam2}...", file=sys.stderr)
    sam2_reads = parse_sam_file(args.sam2)
    for k,v in sam2_reads.items():
        sam2_reads[k] = list(sorted(v))
    
    # Get all unique read names from both files
    all_read_names = set(sam1_reads.keys()) | set(sam2_reads.keys())
    
    # Open output file or use stdout
    output_file = open(args.output, 'w') if args.output else sys.stdout
    
    try:
        # Process each read name
        for read_name in sorted(all_read_names):
            if (read_name in sam1_reads) and (read_name in sam2_reads):
                if sam1_reads[read_name] == sam2_reads[read_name]:
                    continue
            # Print sam1 lines for this read
            if read_name in sam1_reads:
                output_file.write(f"sam1:\n")
                for line in sam1_reads[read_name]:
                    output_file.write(f"{line}\n")
            
            # Print sam2 lines for this read
            if read_name in sam2_reads:
                output_file.write(f"sam2:\n")
                for line in sam2_reads[read_name]:
                    output_file.write(f"{line}\n")
            
            # Add blank line between different read names
            output_file.write("\n")
    
    finally:
        if args.output:
            output_file.close()
    
    print(f"Processed {len(all_read_names)} unique read names", file=sys.stderr)

if __name__ == "__main__":
    main()
