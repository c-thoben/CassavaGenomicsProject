### Cassava Genomics Project at PuckerLab ###
### PuckerLab website
### Cassava Genomics Projekt repository: https://github.com/c-thoben/CassavaGenomicsProject ###

__version__ = "v0.2"

__usage__ = """
                python genetic_map_to_fasta.py \
                    --map <FULL_PATH_TO_GENETIC_MAP_FILE> \
                    --contigs <FULL_PATH_TO_CONTIGS_FILE>
                    --output <BASE_PATH_TO_OUTPUT_FILE> \
                    [--sim <MINIMUM_SIMILARITY_BEST_HIT]
                    [--score <MINIMUM_SCORE_BEST_HIT]
                    
                Input:
                --map: genetic markers
                --fasta: FASTA file containing assembly contigs
                --output: Base path to output files (without extension)

                Output:
                <output>.fasta -> FASTA file containing marker sequences
                <output>_blastn_marker_contig_mapping.txt -> BLASTN output
                <output>_mapped_contigs.csv -> Mapped markers, compatible with ALLMAPS merge command

                Mapping of genetic markers (--map) to assembly contigs (--contigs). 
                Genetic markers are expected in the format of the composite genetic map from Manihot
                esculenta Crantz by the ICGMC (File S2, https://doi.org/10.1534/g3.114.015008).


                Please cite the COL40 assembly publication when using the script: https://doi.org/10.3390/genes13071131
                For questions, contact c.thoben@tu-braunschweig.de for help.
            """

import argparse
import sys, os, subprocess, re, csv
from operator import itemgetter


def merge_and_adjust(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            columns = line.strip().split('\t')
            identifier = f"{columns[0]}_{columns[1]}_{columns[2]}_{columns[3]}"
            sequence = re.sub(r'{([A-Z]+),.+}', r'\1', columns[5])
            outfile.write(f">{identifier}\n{sequence}\n")


def from_roman(num):
    roman_numerals = {'I':1, 'V':5, 'X':10, 'L':50, 'C':100, 'D':500, 'M':1000}
    result = 0
    for i,c in enumerate(num):
        if (i+1) == len(num) or roman_numerals[c] >= roman_numerals[num[i+1]]:
            result += roman_numerals[c]
        else:
            result -= roman_numerals[c]
    return result


def get_chromosome(marker_name):
    match = re.search(r'_chromosome([XVI]+)_', marker_name)
    if match:
        chromosome_number = from_roman(match.group(1))
        return chromosome_number
    else:
        sys.stderr.write('Chromosome not found: ' + marker_name)
        return 0
    

def parse_blastn_output(blastn_output, csv_output, min_sim, min_score):
    data = []

    with open(blastn_output, 'r') as infile:
        
        hits = {}
        for line in infile:
            fields = line.strip().split('\t')

            # Scaffold ID, scaffold position, LG, genetic position
            query_id = fields[0]
            params =  { 
                    'score': float(fields[-1]),
                    'similarity': float(fields[2]),
                    'Scaffold ID': fields[1],
                    'scaffold position': fields[8], # start: position of marker on scaffold
                    'LG': get_chromosome(query_id), 
                    'genetic position': query_id.split("_")[3]
                    }   
            try:
                hits[ query_id ].append( params )
            except:
                hits.update( {  query_id: [ params ] } )
        
    data = [] 
    sims = []
    scores = []       
    for hit in list( hits.values() ):
        sorted_hits = list( sorted( hit, key=itemgetter('score') ) )[::-1]
        
        if sorted_hits[0]['similarity'] > min_sim and sorted_hits[0]['score'] > min_score:
            data.append({
                'Scaffold ID': sorted_hits[0]['Scaffold ID'],
                'scaffold position': sorted_hits[0]['scaffold position'],
                'LG': sorted_hits[0]['LG'],
                'genetic position': sorted_hits[0]['genetic position']
            })
        sims.append(sorted_hits[0]['similarity'])
        scores.append(sorted_hits[0]['score'])
     
    print("mean similarity of best hits: ", str( sum(sims) / len(sims) )) 
    print("mean score of best hits: ", str( sum(scores) / len(scores) )) 
             
    with open(csv_output, 'w', newline='') as csvfile:
        fieldnames = ['Scaffold ID', 'scaffold position', 'LG', 'genetic position']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in data:
            writer.writerow(row)


def main():

    # get arguments
    parser = argparse.ArgumentParser(description="Create input file for ALLMAPS merge command from genetic markers.")
    parser.add_argument("--map", required=True, help="Full path to genetic markers")
    parser.add_argument("--contigs", required=True, help="Full path to FASTA file containing assembly contigs")
    parser.add_argument("--output", required=True, help="Base path to output files (without extension)")
    parser.add_argument("--sim", help="Minimum siliarity for marker best hits", type=float)
    parser.add_argument("--score", help="Minimum score for marker best hits", type=float)
    args = parser.parse_args()

    fasta_path = args.output + ".fasta"
    out = fasta_path.split(".")[0] + "_blastn_marker_contig_mapping.txt"
    out_map = fasta_path.split(".")[0] + "_mapped_contigs.csv"

    # create FASTA of genetic map
    merge_and_adjust(args.map, fasta_path)

    # BLAST with markers against contigs
    if args.contigs:
       
        cmd = f"blastn -query {fasta_path} -subject {args.contigs} -outfmt 6 -out {out}"
        if not os.path.isfile( out ):
            subprocess.run([cmd], shell=True)
        elif os.path.getsize( out ) < 1000:
            subprocess.run([cmd], shell=True)

    # create contig genetic map from blast result
    min_sim = args.sim if args.sim else 0
    min_score = args.score if args.score else 0
    parse_blastn_output(out, out_map, min_sim, min_score)


if __name__ == '__main__':
	main()	