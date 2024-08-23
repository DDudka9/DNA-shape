### FIRST: RUN THIS ###

import subprocess
import regex
import csv
import statistics
from Bio import SeqIO
from Bio import SearchIO

# change path to folder where you downloaded fasta reads
path = "/Users/damian/Downloads/sratoolkit.3.1.1-mac-arm64/bin/SRR11606870/fasta/"

# get all records with the consensus ATTCGTTGGAAACGGGA motif (CENP-B box) specific to minor satellite
input_iterator = SeqIO.parse(path + "SRR11606870.fasta", "fasta")
Min_iterator = (record for record in input_iterator if "ATTCGTTGGAAACGGGA" in record.seq
                    and "GAAAACTGAAAA" not in record.seq)  # avoid mixed zone of Maj and Min sat
SeqIO.write(Min_iterator, path + "SRR11606870_Min.fasta", "fasta")

# make blast database locally (modify path to match the path with downloaded fasta reads)
cmd = "makeblastdb -in /Users/damian/Downloads/sratoolkit.3.1.1-mac-arm64/bin/SRR11606870/fasta/SRR11606870_Min.fasta -dbtype nucl"
subprocess.run(cmd, shell=True)

### ------------------------------------------------------------------------------------------###

### SECOND: Copy the SRR11606870_111923.fasta file from github cloned folder into the path directory
# if you havent done it yet ###

### ------------------------------------------------------------------------------------------###

### THIRD: RUN THE REST OF THE CODE ###

# blast representative Type 1 Continuous Min sat array (111923) from Packiaraj and Thakur 2024 Genome Biol (PMID: 38378611)
path2 = path + "SRR11606870_Min.fasta"
homo_cont_array_path = path + "SRR11606870_111923.fasta"  # make sure this read fasta file is in the same folder
outpath = path + "SRR11606870_Min.xml"
cmd3 = "blastn -db " + path2 + " -query " + homo_cont_array_path + " -outfmt 5 -out " + outpath
subprocess.run(cmd3, shell=True)

# gets blast results
blast_qresult = SearchIO.read(outpath, "blast-xml")
hit_keys = blast_qresult.hit_keys # gets IDs of hits
hit_iterator = SeqIO.parse(path + "SRR11606870_Min.fasta", "fasta")  # generator to iterate fast through
hit_iterator_blasted = (record for record in hit_iterator if record.id in hit_keys)  # get sequences of hits
SeqIO.write(hit_iterator_blasted, path + "SRR11606870_Min_blasted.fasta", "fasta")  # save to file

# convert reads to oneliners
with open(path + "SRR11606870_Min_blasted.fasta", "r") as f:
    reads = f.readlines()

d = {}
header = ">"
seq = ""
for line in reads:
    if line.startswith(">"):
        head = line.lstrip(">").rstrip("\n")
        seq = ""
    else:
        seq = seq + line.rstrip("\n").upper()
        d[head] = seq

with open(path + "SRR11606870_Min_blasted_oneliner.fasta", 'w') as convert_file:
    for k, v in d.items():
        convert_file.write(">" + k.rstrip("\n"))
        convert_file.write("\n")
        convert_file.write(v + "\n")

# find >=4 AT stretches
Min_sat_ATstretches = {}
density_list = []
with open(path + "SRR11606870_Min_blasted_oneliner.fasta", "r") as f:
    reads = f.readlines()
    for line in reads:
        if line.startswith(">"):
            id = line
            read = ""
        else:
            read = line.rstrip("\n")
            ATstretches = regex.findall(r"[AT]{4,}", read, overlapped=False)
            nr_stretches = len(ATstretches)
            density = len(ATstretches)/(len(read)/234)
            density_list.append(density)
            Min_sat_ATstretches[id] = read, ATstretches, nr_stretches, density

# write into file
with open(path + "SRR11606870_Min_ATstretches_average.fasta", 'w') as convert_file:
    for k, v in Min_sat_ATstretches.items():
        convert_file.write(k.rstrip("\n"))
        convert_file.write("\n")
        convert_file.write(v[0])
        convert_file.write("\n")
        convert_file.write(str(v[1]))
        convert_file.write("\n")
        convert_file.write("Number of stretches in array = %s" % str(v[2]))
        convert_file.write("\n")
        convert_file.write("Density of stretches per 234bp = %s" % str(v[3]))
        convert_file.write("\n")
    convert_file.write("\nAverage density of stretches per 234bp = %s" % statistics.mean(density_list))
    convert_file.write("\nAll densities of stretches per 234bp = %s" % density_list)

# Save densities into file
with open(path + "SRR11606870_Min_ATstretches_average.csv", "w") as f:
    write = csv.writer(f)
    write.writerow(density_list)

# Compute most narrow tetranucleotides from Rohs et al., 2009 Nature (PMID: 19865164)
ATs = ["AAAT", "AATA", "AATC", "AATT", "AAAA", "AAGT", "GAAT", "GAAA", "TAAT", "AAAC", "ATAA",
       "AGAT", "AAGA", "AGTT", "AGAA", "AAAG", "ATAG", "GAAC", "CGTT", "TATA"]
Min_sat_ATstretches = {}
density_list = []
stretches = {}
stretches_total = {}
for s in ATs:
    stretches_total[s] = 0

density_list = []
with open(path + "SRR11606870_Min_blasted_oneliner.fasta", "r") as f:
    reads = f.readlines()
    total_length = 0
    for line in reads:
        if line.startswith(">"):
            id = line
            read = ""
        else:
            read = line.rstrip("\n")
            nr_stretches = 0
            for stretch in ATs:
                ATstretches = regex.findall(stretch, read, overlapped=True)
                stretches[stretch] = len(ATstretches)
                nr_stretches += len(ATstretches)
                stretches_total[stretch] += len(ATstretches)
            length = len(read)
            total_length += length
            density = nr_stretches/(length/234)
            density_list.append(density)
            Min_sat_ATstretches[id] = read, stretches, nr_stretches, density, length

with open(path + "SRR11606870_Min_teranucleotides_average.fasta", 'w') as convert_file:
    for k, v in Min_sat_ATstretches.items():
        convert_file.write(k.rstrip("\n"))
        convert_file.write("\n")
        convert_file.write(v[0])
        convert_file.write("\n")
        convert_file.write(str(v[1]))
        convert_file.write("\n")
        convert_file.write("Number of stretches in array = %s" % str(v[2]))
        convert_file.write("\n")
        convert_file.write("Density of stretches per 234bp = %s" % str(v[3]))
        convert_file.write("\n")
    convert_file.write("\nAverage density of stretches per 234bp = %s" % statistics.mean(density_list))
    for stretch, number in stretches_total.items():
        convert_file.write("\nAverage density of %s stretches per 234bp = %s" % (stretch, str(number/(total_length/234))))

narrow_ATs = ["AAAT", "AATA", "AATC", "AATT", "AAAA", "AAGT", "GAAT", "GAAA", "TAAT", "AAAC"]

with open(path + "SRR11606870_Min_blasted_oneliner.fasta", "r") as f:
    lines = f.readlines()
    total_length = 0
    wanted_array = False
    for line in lines:
        # select the array you want to analyze (change ID number)
        if "111923" in line:
            wanted_array = True
        if not line.startswith(">") and wanted_array:
            seq = line.rstrip("\n")
            break

bin_size = 1000
num_bins = (len(seq) + bin_size - 1) // bin_size  # number of bins needed

### Use sliding window to get tetras overlapping and their positions

# Get dict where each tetra will get its frequency
Min_freq_dict = {bin_index: {tetra: 0 for tetra in narrow_ATs} for bin_index in range(num_bins)}
for i in range(len(seq) - 3):
    # Get current tetra
    current_tetra = seq[i:i + 4]
    # Put it in the current bin
    bin_index = i // bin_size
    # See if that tetra is narrow based on the list of ATs and if so, count it
    if current_tetra in narrow_ATs:
        Min_freq_dict[bin_index][current_tetra] += 1

#  write into file (with help from https://stackoverflow.com/questions/29400631/python-writing-nested-dictionary-to-csv)
with open(path + "SRR11606870_Min_tetranucleotides.csv", "w") as csvfile:
    writer = csv.DictWriter(csvfile, narrow_ATs)
    for key, val in sorted(Min_freq_dict.items()):
        row = {"AAAT": key}
        row.update(val)
        writer.writerow(row)