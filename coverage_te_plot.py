### Cassava Genomics Project at PuckerLab ###
### https://www.tu-braunschweig.de/en/ifp/pbb
### At7 genome sequencing and analysis repository: https://github.com/bpucker/At7
### Cassava Genomics Projekt repository: https://github.com/c-thoben/CassavaGenomicsProject ###

### WARNING: optimized for Manihot esculenta ###

__version__ = "v0.5"

__usage__ = """
python coverage_te_plot.py \
--coverage_file <FULL_PATH_TO_COVERAGE_FILE> \
--te_file <FULL_PATH_TO_TE_FILE> \
--out <FULL_PATH_TO_OUTPUT_FILE> \
--cov <AVERAGE_COVERAGE>
[--res <RESOLUTION, WINDOW_SIZE_FOR_COVERAGE_CALCULATION> 1000]
[--sat <SATURATION, CUTOFF_FOR_MAX_COVERAGE_VALUE> 100.0]
[--num_contigs <NUMBER_OF_CONTIGS_TO_PLOT> 18]
[--max_cov <MAXIMUM_COVERAGE 600]
[--max_chromosome <MAXIMUM_CHROMOSOME_SIZE_BP 55000000]

Input:
--coverage_file: Coverage file
--te_file: Repeats TSV created with Circos genomicDensity function
--cov: Average coverage
--out: Base path to output files (without extension)

Output:
<out>.png -> Coverage plot
<out>_<contig>.png -> Histogram of coverage for each contig/chromosome
<out>_coverage_resolution<res>.tsv -> Average coverage per block

Creates a coverage plot showing the average coverage in blocks of resolution size (--res). The maximum displayed coverage is defined
by the saturation (--sat). Each chromosome is plotted separatly and the average coverage is marked by a red line. For each chromosome,
a histogram is created showing the coverage distribution.

Script adapted from cov_plot.py script, please cite the At7 publication when using the script: https://doi.org/10.1371/journal.pone.0164321
For questions, contact b.pucker@tu-braunschweig.de for help
					"""

import sys, os, re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Colour used to plot coverage
COV_COL = "#089392"
# Colour used to plot TE density
REP_COL = "#AAD9A7"


def natural_sort_key(s):
    """! @brief sort dictionary keys in numerical order """
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]


def load_cov( cov_file ):
	"""! @brief load all information from coverage file """
	
	df = pd.read_csv(cov_file, sep="\t", header=None, names = ["contig", "pos", "cov"])
 
	# only chromosomes
	df = df[df['contig'].str.contains('chr')]

	cov = df.groupby('contig')['cov'].apply(list).to_dict()
	return cov


def generate_plot( cov, out_file, saturation, coverage, resolution,
                  collected_values, ymax, max_value, chromosome_max, xticks, rep_df ):
	"""! @brief generate figure """
	
	fig, ax = plt.subplots( figsize=( 20, 8 ) )
	
	# --- plot values --- #
	max_value = float( min( [ saturation, max_value ] ) )
	factor = 1
 
	for idx, key in enumerate( sorted( cov.keys(), key=natural_sort_key )[:ymax] ):
		y = ymax - ( idx*1.3 )
		x = []

		# plot coverage + red line
		for each in collected_values[ key ]:
			x.append( (y + min( [ 1, ( each / max_value ) ] ) ) * factor )
		ax.plot( x, linestyle=":", color=COV_COL, markersize=0.2  ) #mcolors.CSS4_COLORS["lime"]
		ax.plot( [ 0, len( x ) ], [ ( y+( coverage / max_value ) ) * factor, (y+( coverage / max_value ))*factor ], color="#CF597E" , linewidth=1.5)
		ax.text(0 , (y+1) * factor, str( int( max_value ) ), ha="right", fontsize=5 )
		ax.plot( [ 0, 0 ], [ y*factor, (y+1) *factor ], color="black", linewidth=1, markersize=1 )
		
		# plot density
		repeats = rep_df[rep_df["chr"] == key ]
		x_rep = [0]
		y_rep = [y + repeats.iloc[0,3] * factor]
		min_rep = [y]
		last_val = 0
		for index, row in repeats.iterrows():
			start = row["start"]
			end = row["end"]
			x_rep.append( (start + ( (end-start)/ 2 ) ) / resolution)
			y_rep.append( y + row["value"] * factor)
			min_rep.append(y)
			last_val =  y + row["value"] * factor

		x_rep.append(len(x))
		y_rep.append(last_val)
		min_rep.append(y)
		ax.fill_between(x_rep, min_rep, y_rep, alpha=0.3, color=REP_COL) #"grey"

		# between lines
		for l in range(0,int(max_value+50), 50):
			ax.plot( [ 0, len( x ) ], [  (y + min( [ 1, ( l / max_value ) ] )) * factor , (y + min( [ 1, ( l / max_value ) ] )) * factor ], color="grey", linewidth=0.5 )
  
		# y axis label
		ax.text( 0, (y+0.5) * factor, str( int( max_value / 2 ) ), ha="right", fontsize=5 )
		ax.text( 0, y*factor, "0", ha="right", fontsize=5 )
  
		# chromosome label
		chromosome_pos = chromosome_max / 25.7
		ax.text( chromosome_pos, (y+0.75) * factor, str( key ), ha="right", fontsize=7 )
		
	# style
	ax.set_xlabel( "position on chromosome [ Mbp ]" )
	ax.set_ylabel( "Coverage" )
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.get_yaxis().set_ticks([])
	ax.yaxis.labelpad = 15
	
	# x axis label
	ax.set_xlim( 0, chromosome_max)
	ax.set_xticks([x * 1000000 / 1000 for x in xticks])
	ax.set_xticklabels( [ str(x) for x in xticks ] )

	plt.subplots_adjust( left=0.025, right=0.989, top=0.99, bottom=0.1 )
	fig.savefig( out_file, dpi=300 )


def generate_hist( cov_values, outputfile, saturation ):
	"""! @brief generate coverage histogram """
	
	values = []
	for each in cov_values:
			if each > saturation:
					values.append( saturation )
			else:
					values.append( each )
		
	fig, ax = plt.subplots()
	
	ax.hist( values, edgecolor='white', bins = 80 )
	ax.set_xlim( 0, saturation )
	ax.set_xlabel( "coverage" )
	ax.set_ylabel( "count" )
	
	fig.savefig( outputfile, dpi=300 )


def main( arguments ):
	 
	# --- script input ---
	cov_file = arguments[ arguments.index( '--coverage_file' ) + 1 ]
	te_file = arguments[ arguments.index( '--te_file' ) + 1 ]
	out_file = arguments[ arguments.index( '--out' ) + 1 ] + ".png"
	coverage = float(arguments[ arguments.index( '--cov' ) + 1 ]) 

	# --- plot options --- #
	resolution = int(arguments[ arguments.index( '--res' ) + 1 ]) if '--res' in arguments else 1000
	saturation = float(arguments[ arguments.index( '--sat' ) + 1 ]) if '--sat' in arguments else 100.0
	ymax = int(arguments[ arguments.index( '--num_contigs' ) + 1 ]) if '--num_contigs' in arguments else 18
	chromosome_max = int(arguments[ arguments.index( '--max_chromosome' ) + 1 ]) if '--max_chromosome' in arguments else 55000000
	max_value = 600 if not '--max_cov' in arguments else int( arguments[ arguments.index( '--max_cov' ) + 1 ] )
	xticks = [0, 10, 20, 30, 40, 50, 55]

	# --- adjust chromosome length to resolution --- #
	chromosome_max = chromosome_max / resolution

	# --- load repeat information --- #
	print("--- load repeat information ---")
	rep = pd.read_csv(te_file, sep="\t")
	rep.columns = ["chr", "start", "end", "value"]
	rep["chr"] = rep["chr"].str.replace('"', '')
	print(rep.head)

   	# --- get average coverage per block --- #
	print("--- get average coverage per block ---")
	out_file_list = "%s_coverage_resolution%s.tsv" % (out_file.split(".")[0], str(resolution))
 
	if not os.path.isfile(out_file_list):
		# --- load contigs coverage file --- #
		cov = load_cov( cov_file )

		# --- generate coverage histograms per chromosome --- #
		for key in cov.keys():
			outputfile = out_file.split(".")[0] + "_" + key + ".png"
			generate_hist( cov[ key ], outputfile, saturation )
     
		# --- calculate average coverage in blocks --- #
		collected_values = {}
		for idx, key in enumerate( sorted( cov.keys() ) ):
			y = ymax-idx-1
			x = []
			blocks = [ cov[ key ] [ i : i + resolution ] for i in range( 0, len( cov[ key ] ), resolution ) ]
			for block in blocks:
				x.append( np.mean( block ) )
			max_value = max( [ max_value, max( x ) ] )
			collected_values.update( { key: x } )
		df = pd.DataFrame.from_dict(collected_values, orient="index")
		df.to_csv(out_file_list, sep="\t")
	else:
		df = pd.read_csv(out_file_list, sep="\t", index_col=0)
		d = df.to_dict("split")
		max_value = max([max(v) for v in d["data"] ])
		collected_values = dict(zip(d["index"], d["data"]))
		for c in collected_values:
			collected_values[c] = [i for i in collected_values[c] if not np.isnan(i)]
		cov = collected_values
  
	# --- generate coverage plot --- #
	print("--- generate coverage plot ---")
	generate_plot( cov, out_file, saturation, coverage, resolution, collected_values, ymax, max_value, chromosome_max, xticks, rep )

	
if '--coverage_file' in sys.argv and '--te_file' in sys.argv and '--out' in sys.argv and '--cov' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
