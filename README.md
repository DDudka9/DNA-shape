# DNA-shape
Based on long read SRR11606870 assembly of _Mus musculus_ (Ton et al., Scientific Data 2020; PMID: 33203859), 
and using Packiaraj and Thakur 2024 Genome Biol (PMID: 38378611) major and minor satellite reads as templates, the scripts allow computing:
1. Number of contigous A/T stretches (default set to length of 4)
2. Number of tetranucleotides associated with narrow minor groove (Rohs et al., 2009 Nature; PMID: 19865164)

Requirements:
1. Operating system: MacOS Ventura
2. Language: Python 3.11
3. Modules: biopython==1.83; mkl-service==2.4.0; regex==2024.5.15

Instructions:
1. Download sra-tools (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
    Useful wiki: https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
    Useful blog: https://edwards.flinders.edu.au/fastq-dump/
2. Using sra-tools download SRR11606870 dataset from https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR11606870&display=metadata
3. Prefetch dataset (prefetch SRR11606870)
4. Download all reads unsorted in fasta format (fasterq-dump pathtoSRAfile --outdir pathtoSRAfile/fasta --fasta-unsorted)
5. Create a new virtual environment using miniconda (https://docs.anaconda.com/miniconda/)
6. Install packages from requirements.txt file
7. Follow instructions inside majsat.py and minsat.py scripts



