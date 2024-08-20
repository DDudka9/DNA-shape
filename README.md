# DNA-shape

The scripts allow computing the number of contigous A/T stretches (default set to length of 4) from SRA assemblies. 

To analyse stretches in Major and Minor satellite arrays from SRR11606870 long read _Mus musculus_ assembly 
(Ton et al., Scientific Data 2020; PMID: 33203859) follow this:
1. Download sra-tools (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
    Useful wiki: https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
    Useful blog: https://edwards.flinders.edu.au/fastq-dump/
2. Using sra-tools download SRR11606870 dataset from https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR11606870&display=metadata
3. Prefetch dataset (prefetch SRR11606870)
4. Download all reads unsorted in fasta format (fasterq-dump pathtoSRAfile --outdir pathtoSRAfile/fasta --fasta-unsorted)
5. Create a new virtual environment using miniconda (https://docs.anaconda.com/miniconda/)
6. Install packages from requirements.txt file
7. Follow instructions inside majsat.py and minsat.py scripts



