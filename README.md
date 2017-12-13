# CLAM_m6A_Snakemake
Snakemake pipeline of CLAM processing m6A RIP-seq (MeRIP)


## Instructions

1.a Running from download of SRA:

Configure the "config.yaml" file by running "configure.py"
```
python configure.py GSE12345
```

1.b Running from local reads/fastq files

Alternatively, just make a "{project}/{reads}" directory, and put the corresponding files in
```
python make_ln.py
```

2. Run the pipeline
```
source activate py36
PROJECT={project} snakemake -n
```
