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
## Update 08/30/2018
To accomadate the changed GEO API (the code actually updated back on March 2018), use the following command to prepare a GEO run:
```
python2 scripts/geo_downloader/getSRApath.py GSE52600 projects/GSE52600/config/sra_dir.txt 
```
Once the `sra_dir.txt` file is prepared, run
```
python2 configure.py GSE52600
```
Then submit it by 
```
./prepare_individual_run.sh GSE52600
```