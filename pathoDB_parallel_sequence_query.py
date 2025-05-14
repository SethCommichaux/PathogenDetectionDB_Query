import sys
import os
from Bio import SeqIO
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="Run the analysis with user-defined parameters.")

parser.add_argument("-b","--bootstrap", type=int, default=None, help="Bootstrap iteration number")
parser.add_argument("-f","--query_fasta", default=None, help="Path to the query FASTA file")
parser.add_argument("-t","--type_query", choices=["phage", "plasmid", "genes"], help="Type of query sequence")
parser.add_argument("-p", "--percent_identity", type=float, default=None, help="Minimum percent identity")
parser.add_argument("-c", "--query_coverage", type=float, default=None, help="Minimum percent query coverage")
parser.add_argument("-g", "--query_genes_found", type=float, default=None, help="Minimum proportion of query genes found")

args = parser.parse_args()

bootstrap = args.bootstrap
query_fasta = args.query_fasta
type_query = args.type_query
percent_identity = args.percent_identity
query_coverage = args.query_coverage
min_perc_genes_found = args.query_genes_found

number_genes = sum([1 for i in SeqIO.parse(query_fasta,'fasta')])

path_to_SRA_assemblies = '/isilon/NCBI/SRAassemblies/skesa_contigs/'
PathoDB_assemblies = [i for i in os.listdir(path_to_SRA_assemblies)]
slice = int(len(PathoDB_assemblies)/1000)
PathoDB_assemblies = PathoDB_assemblies[bootstrap*slice-slice:bootstrap*slice]
path_to_SRA_metadata = '/lustre/projects/SethCommichaux/SRA_metadata/Esherichia_Shigella_Listeria_Salmonella_Klebsiella_Campylobacter_SRA_metadata.txt'

SRA_metadata = {}

for h,i in enumerate(open(path_to_SRA_metadata)):
	if h == 0:
		SRA_header = i.strip()
	else:
		sra = i.strip().split('\t')[8]
		SRA_metadata[sra] = i.strip()

# blast and mash query parameters are tailored to provide similar results
for f in PathoDB_assemblies:
	os.system(f'echo {f} >> tmp{bootstrap}.out')
	if type_query == "plasmid":
		# use for plasmids when you want highly similar sequences in terms of sequence similarity, length, and synteny
		os.system(f'mash dist -d 0.01 -k 25 -s 100000000000 -i {query_fasta} {path_to_SRA_assemblies}/{f} >> tmp{bootstrap}.out')
	elif type_query == "phage":
		# use for phage queries where sequence is likely integrated into chromosome
		awk_command = "awk -F'\\t' '$5 <= 20 {print $0}'"
		command = f"blastn -query {query_fasta} -subject {path_to_SRA_assemblies}/{f} -perc_identity {percent_identity} -qcov_hsp_perc {query_coverage} -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' | {awk_command} >> tmp{bootstrap}.out"
		os.system(command)
	elif type_query == 'genes':
		command = f"blastn -query {query_fasta} -subject {path_to_SRA_assemblies}/{f} -perc_identity {percent_identity} -qcov_hsp_perc {query_coverage} -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' >> tmp{bootstrap}.out"
		os.system(command)
	else:
		sys.exit("Must select phage, plasmid, or genes as type of query sequence!")

with open(f'tmp{bootstrap}.sra_meta','w') as out:
	out.write(f'{SRA_header}\n')
	sras = set()
	results = defaultdict(set)
	for i in open(f'tmp{bootstrap}.out'):
		tmp = i.strip().split('\t')
		if len(tmp) == 1:
			sra = tmp[0].split('_')[0]
		elif len(tmp) > 1:
			if type_query == 'genes':
				results[sra].add(tmp[0])
			else:
				sras.add(sra)
	if type_query == 'genes':
		for k,v in results.items():
			if len(v)/float(number_genes) >= min_perc_genes_found:
				sras.add(k)
	for x in sras:
		sra_meta = SRA_metadata.get(x,f'No SRA data found for {x}')
		out.write(f'{sra_meta}\n')
