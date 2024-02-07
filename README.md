# diploidocus_setonix

Nextflow based workflow to parallelise https://github.com/slimsuite/diploidocus for HPC deployment


1. Clone repo into your environment
2. Edit nextflow.config file to change container locations to where you have them stored (not necessary for Oceanomics users). Highly recommend using a Seqera account so that the pipeline execution can be easily monitored (https://training.nextflow.io/basic_training/seqera_platform/), else delete the tower parameters in nextflow.config.
3. **Choose which script to run based on gzipped assembly file size: small (under 600Mb), medium (600Mb-4Gb), large (4-7Gb).** If running in diploid mode the assembly file size is based on one file.
5. Run chosen main.nf using the slurm template provided, the following params can be added to the nextflow run command by adding **'--param /path/to/file'**. Without these parameters the script will search for the files in the default places and will make them if not found:

## Run parameters
   
```
Assembly and reads are the only necessary input, all other files will be made by the workflow if needed
--asssembly Path to assembly file (necessary) 

--projectDir Where output directories will be made and the deafult directory where files will be searched for [default '.']

--hifireads Path to all reads files [default "$params.projectDir/reads/*"], you have to put "" around
the directory location if using '*' e.g. --hifireads "/path/to/*.fastq"

--bam Bam file for assembly. If running in diploid mode (using 2 assembly files) then set as
the directory where both bam files can be found. [default "$params.projectDir/bam"]

--busco Path to full_table.tsv output from busco [default "$params.projectDir/busco/**full_table.tsv"]

--lineage Lineage to use for BUSCO analysis (only necessary if busco file is not given) [default 'actinopterygii_odb10']

--depthsizer Path to folder with depthsizer output files [default search in ("${params.projectDir}/depthsizer/"]

--kat Path to folder with kat and selfkat output files  [default search in ("${params.projectDir}/kat/"]

--purge_hap Path to folder with purge_haplotigs output files  [default search in ("${params.projectDir}/purge_haplotigs/"]
```




## Diploid mode

Using diploid mode, the haplotype input will be pre-processed with its paired haplotype file. To use this mode the **--asssembly** input stays as one file and the **--PairedFactor** [default "'hap1', 'hap2'"] parameter can be set as the difference between the files. The other haplotype file will be searched for in the same directory and if found diploid pre-processing will be performed. If no file is found will proceed in haploid mode. For example, to run diploid mode using 2 haplotype files run: 

nextflow run main.nf --assembly /path/to/my_assembly.hap1.fa    #where my_assembly.hap2.fa file can be found in the same directory

Or to run using a different naming convention:

nextflow run main.nf --assembly /path/to/my_assembly.mat.fa --PairedFactor "'mat', 'pat'"   #where my_assembly.pat.fa file can be found in the same directory



**Diploid mode**
![image](https://github.com/jadedavis5/diploidocus_nextflow_HPC/assets/111946376/b1d16824-4b50-47ef-a8c1-42a958b53ed4)

**Haploid mode**
![image](https://github.com/jadedavis5/diploidocus_nextflow_HPC/assets/111946376/c894483d-dc88-4ebe-9673-3b3ed69d73bc)





