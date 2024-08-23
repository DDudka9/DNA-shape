# DNA-shape
Based on long read SRR11606870 assembly of _Mus musculus_ genome (Ton et al., Scientific Data 2020; PMID: 33203859), 
and using Packiaraj and Thakur 2024 Genome Biol (PMID: 38378611) major and minor satellite reads as templates, the scripts allow computing:
1. Number of contiguous A/T stretches (default set to length of 4)
2. Number of tetranucleotides associated with narrow minor groove (Rohs et al., 2009 Nature; PMID: 19865164)

Requirements:
1. Operating system tested: MacOS (M1) Ventura
2. Language: Python 3.11
3. Installer: miniconda
4. Installation time: minutes (10min)
5. Code run time: minutes (10-30min)
6. Modules: blast - version 2.6.0, biopython - version 1.83; mkl-service - version 2.4.0; regex - version 2024.5.15

Tools to download (we recommend using terminal window/bash):
1. Download miniconda (https://docs.anaconda.com/miniconda/)
2. Download blast via miniconda (https://anaconda.org/bioconda/blast)
3. Download sra-tools (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
    Useful wiki: https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
    Useful blog: https://edwards.flinders.edu.au/fastq-dump/
4. We recommend using a free software Pycharm 2023.2.5 Community Edition (https://www.jetbrains.com/pycharm/) to run majsat.py and minsat.py scripts

Prepare dataset (we recommend using terminal window/bash):
1. Use sra-tools to download SRR11606870 dataset from https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR11606870&display=metadata
2. Prefetch dataset (prefetch SRR11606870)
3. Download all reads unsorted in fasta format (fasterq-dump pathtoSRAfile --outdir pathtoSRAfile/fasta --fasta-unsorted)
4. Rename to: SRR11606870.fasta
5. Copy the two fasta files representing major (SRR11606870_2342980.fasta) and minor satellite (SRR11606870_111923.fasta) reads from https://github.com/DDudka9/DNA-shape.git into the folder with SRR11606870.fasta file

Run the scripts (in Pycharm 2023.2.5 Community Edition)
1. Clone the Github repository: (Git / Clone / url https://github.com/DDudka9/DNA-shape.git)
2. Create new interpreter: Python interpreter (bottom right corner) / Add New Interpreter / Add Local Interpreter / Conda Environment / Create New Environment (provide path to miniconda; select Python 3.11)
3. Select the "requirement.txt" file and click "Install requirements"
4. Follow instructions inside majsat.py and minsat.py scripts (use the PyCharm in-built Python Console to run subsequent parts of the scripts by copy-pasting the code into the console -> press return)
5. The output should appear in the folder with SRR11606870_2342980.fasta and SRR11606870_111923.fasta files 

Expected output files:
1. SRR11606870_Maj_2342980_tetranucleotides.csv - Spreadsheet where each column represents a number of tetranucleotides with narrow major groove (order: AAAT; AATA; AATC; AATT; AAAA; AAGT; GAAT; GAAA; TAAT; AAAC) per 1kb along a representative major satellite array (SRR11606870_2342980)
2. SRR11606870_Min_111923_tetranucleotides.csv - Spreadsheet where each column represents a number of tetranucleotides with narrow minor groove (order: AAAT; AATA; AATC; AATT; AAAA; AAGT; GAAT; GAAA; TAAT; AAAC) per 1kb along a representative minor satellite array (SRR11606870_111923)
3. SRR11606870_Maj_tetranucleotides_average.fasta - Number of tetranucleotides with narrow major groove per 234bp of 500 major satellite arrays (find averages at the end of the file)
4. SRR11606870_Min_tetranucleotides_average.fasta - Number of tetranucleotides with narrow minor groove per 234bp of 500 minor satellite arrays (find averages at the end of the file)
5. SRR11606870_Maj_ATstretches_average.fasta - Number of AT stretches (default: minimum 4) per 234bp of 500 major satellite arrays (find averages at the end of the file)
6. SRR11606870_Min_ATstretches_average.fasta - Number of AT stretches (default: minimum 4) per 234bp of 500 minor satellite arrays (find averages at the end of the file)

You can modify the scripts (array ID) to run any other array or use a different dataset.