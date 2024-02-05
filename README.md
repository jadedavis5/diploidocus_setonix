# diploidocus_nextflow_HPC

1. Clone repo into your environment
2. Edit nextflow.config file to change container locations to where you have them stored (not necessary for Oceanomics users)
3. Run main.nf using the slurm template provided, the following params can be added as to the nextflow run command by adding '--param /path/to/file'. Without these parameters the script will search for the files in the default places and will make them if not found:

***Run parameters***
   
**--asssembly** Path to assembly file (necessary) 

**--projectDir** Where output directories will be made and the deafult directory where files will be searched for [default '.']

**--hifireads** Path to all reads files [default "$params.projectDir/reads/*"]

**--bam** Bam file for assembly. If running in diploid mode (using 2 assembly files) then set as the directory where both bam files can be found; don't use this parameter unless you have bam files already made. [default "$params.projectDir/bam"]

**--busco** Path to full_table.tsv output from busco [default "$params.projectDir/busco/**full_table.tsv"]

**--lineage** Lineage to use for BUSCO analysis (only necessary if busco file is not given) [default actinopterygii_odb10]

**--depthsizer** Path to folder with depthsizer output files [default search in ("${params.projectDir}/depthsizer/"]

**--kat** Path to folder with kat and selfkat output files  [default search in ("${params.projectDir}/kat/"]

**--purge_hap** Path to folder with purge_haplotigs output files  [default search in ("${params.projectDir}/purge_haplotigs/"]


