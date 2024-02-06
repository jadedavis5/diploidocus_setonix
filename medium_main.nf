#!/bin/bash


/*
 * pipeline input parameters
 */
params.projectDir = '.'

def determineFactor(fileName, factorsToCheck) {
    def matchingFactor = factorsToCheck.find { fileName.contains(it) }
    return matchingFactor
}

params.PairedFactor = "'hap1', 'hap2'"

params.assembly = "$params.projectDir/assembly/file"
params.hifireads = "$params.projectDir/reads/*"  //Input reads directory
params.bam = ''  //Give a bam file if primary mode and give directory if you want haplotype mode (if no bams exist leave this empty)
params.busco = ''  //Input full_table.tsv file
params.lineage = 'actinopterygii_odb10'
params.depthsizer = ''  //Input scdepth file
params.kat = ''  //Input directory with files
params.purge_hap = ''  //Input directory with files


log.info """
DIPLOIDOCUS-NF PIPELINE
==================
"""
.stripIndent()

workflow {
//Make directories
        ['busco', 'reads', 'bam', 'depthsizer', 'kat', 'purge_haplotigs'].each {
                "mkdir -p ${params.projectDir}/${it}".execute()
        }
//Check if Hifi Reads exist and exit if not
        hifireads_ch = Channel.fromPath(params.hifireads, checkIfExists: true).collect()

        assembly_ch = Channel.fromPath("$params.assembly", checkIfExists: true)

//Check if paired file naming convention is used and if so find matching assembly file and check for pre-made bams
        def pairedFactorString = params.PairedFactor
        def factorsToCheck = Eval.me("[${pairedFactorString}]")
        basename_ch = Channel.fromPath("$params.assembly", checkIfExists: true).map { file -> file.baseName }

        def assembly_file = new File(params.assembly)
        def assembly_basename = assembly_file.baseName
        def assemblydir = assembly_file.getAbsoluteFile().getParentFile().toString()
        def assemblyExistsMap = [:]
        def bamExistsMap = [:]
        def foundAssemblyFiles = [:]
        def foundBamFiles = [:]

        def haplotype_basename = factorsToCheck.inject(assembly_file.baseName, { result, factor -> result.split("\\.${factor}")[0] })
        for (assembly in factorsToCheck) {
                assemblyExistsMap[assembly] = false
        }
        if (determineFactor(assembly_file.toString(), factorsToCheck) != null) {
                println('Haplotype naming found in assembly name, searching for matching file')
                for (assembly in factorsToCheck) {
                def assemblies = new File(assemblydir).listFiles(new FilenameFilter() {
                        boolean accept(File dir, String name) {
                                return name.toString().contains(assembly) && name.toString().contains(haplotype_basename)
                }
        })

        def assemblyExists = assemblies?.size() > 0
        assemblyExistsMap[assembly] = assemblyExists

        if (assemblyExists) {
            foundAssemblyFiles[assembly] = assemblies[0]
        }
        def bamFile = new File(params.bam ?: "$params.projectDir/bam").listFiles()
                .sort { file -> file.name.contains('sort') ? 0 : 1 } //prioritize finding a bam file that is sorted
                .find { file -> file.name.endsWith(".bam") && file.name.contains(assembly) && file.name.contains(haplotype_basename) }

        def bamExists = bamFile != null
        bamExistsMap[assembly] = bamExists
        if (bamExistsMap) {
                foundBamFiles[assembly] = bamFile
                        }
                }
        }
//Create assembly_ch and bam_ch based on assembly and bam existances

        if (assemblyExistsMap[factorsToCheck[0]] && assemblyExistsMap[factorsToCheck[1]]) {
                println('Paired assemblies found\n========================\nRUNNING HAPLOTYPE PRE-PROCESSING\n========================')

                def assembly_type = determineFactor(assembly_file.toString(), factorsToCheck)
                hapbasename_ch = Channel.fromPath("$params.assembly", checkIfExists: true).map { file -> file.baseName.split("\\.${assembly_type}")[0] }

                assembly1_ch = Channel.fromPath(foundAssemblyFiles[factorsToCheck[0]]?.getCanonicalPath())
                assembly2_ch = Channel.fromPath(foundAssemblyFiles[factorsToCheck[1]]?.getCanonicalPath())

                if (!bamExistsMap[assembly_type]) {
                        println('Bam file for haplotype mode not found, will create')
                        diploid_file = haplotypemapping(assembly1_ch, assembly2_ch, hifireads_ch, hapbasename_ch)
                        haplotypesplit(diploid_file, assembly1_ch, assembly2_ch, hapbasename_ch, assembly_type, factorsToCheck)
                        sortedbam_ch = haplotypesplit.out.assembly_file
                } else {
                        bam_ch = Channel.fromPath(foundBamFiles[assembly_type]?.getCanonicalPath())
                        def usingbam = new File(foundBamFiles[assembly_type]?.getCanonicalPath())
                        if (!usingbam.name.contains("sort")) {
                                sortedbam_ch = bam_sortindex(bam_ch, hapbasename_ch, assembly_type)

                        } else {
                                println('Sorting bam file found/given')
                                sortedbam_ch = bam_ch
                        }
                }
                reads_ch = readpartition(sortedbam_ch, hapbasename_ch, assembly_type)

        } else {
                println('No haplotype pairs found- process assembly file in primary mode')
                def assembly_type = ''
                reads_ch = Channel.fromPath("$params.hifireads", checkIfExists: true)

        //Search for bam/create and sort if needed
                def bamDirectory = new File("$params.projectDir/bam/")
                def bamFile = bamDirectory.listFiles()
                        .sort { file -> file.name.contains('sort') ? 0 : 1 } //prioritize finding a bam file that is sorted
                        .find { file -> file.name =~ "${assembly_basename}.*\\.bam" }

                if (params.bam) {
                        println('Using user inputted bam')
                        bam_ch = Channel.fromPath("$params.bam")
                        def usingbam = new File(params.bam)

                        if (!usingbam.name.contains("sort")) {
                                println('Sorting user inputted bam')
                                sortedbam_ch = bam_sortindex(bam_ch, basename_ch, assembly_type)

                        } else {
                                sortedbam_ch = bam_ch
                        }

                } else if (bamFile?.exists()) {
                        println("Found bam in project directory: ${bamFile}")
                        bam_ch = Channel.fromPath(bamFile.absolutePath)

                        if (!bamFile.name.contains("sort")) {
                                println('Sorting the bam file found')
                                sortedbam_ch = bam_sortindex(bam_ch, basename_ch, assembly_type)

                        } else {
                                println('Found bam already sorted')
                                sortedbam_ch = bam_ch
                        }
                } else {
                        println('No bam found, will create')
                        sam_ch = mapping(assembly_ch, hifireads_ch, basename_ch)
                        sortedbam_ch = sam2bam(sam_ch, basename_ch)
                        }
                }

                //Find/create busco
        def fullTableFile = new File("$params.projectDir/busco").listFiles().find {
                it.isFile() && it.name.contains("full_table") && it.name.contains(assembly_basename) && it.extension == 'tsv'
        }
        if (params.busco) {
                println('Using user inputted busco table')
                busco_ch = Channel.fromPath("$params.busco")
        } else if (fullTableFile?.exists()) {
                println("Found busco table in project directory: ${fullTableFile.getCanonicalPath()}")
                busco_ch = Channel.fromPath(fullTableFile.getCanonicalPath())
        } else {
                println('No busco table given or found, will create table')
                busco_ch = busco('EMPTY', assembly_ch, basename_ch)
        }

        //Check if scdepth given or found and create if not
        def scdepthPath = file("${params.projectDir}/depthsizer")
        if (params.depthsizer) {
                scdepth_ch = Channel.fromPath("$params.depthsizer")
        } else {
                def scdepthFound = scdepthPath.listFiles()?.find { file ->
                file.name.contains(assembly_basename) && file.extension == 'scdepth'
        }

                if (scdepthFound) {
                        def scdepthFile = new File(scdepthFound.toAbsolutePath().toString())
                        println("Scdepth file found: ${scdepthFile.absolutePath}")
                        scdepth_ch = Channel.fromPath(scdepthFile.absolutePath)

                } else {
                        println('No scdepth file found, will create')
                        scdepth_ch = depthsizer('EMPTY', assembly_ch, sortedbam_ch, busco_ch, reads_ch, basename_ch)
                }
        }

        //Check for Kat input and run if not given (user should give folder with KAT analysis in it)
        //These are the file names diploidocus will search for hence input isn't accepted unless it meets these guidelines
        def KATfilesToCheck = [
                "${assembly_basename}.selfkat-stats.tsv",
                "${assembly_basename}.selfkat-counts.cvg",
                "${assembly_basename}.kat-stats.tsv",
                "${assembly_basename}.kat-counts.cvg"
        ]

        def katExistsMap = [:]
        def katGivenMap = [:]

        KATfilesToCheck.each { file ->
                katExistsMap[file] = new File("${params.projectDir}/kat/$file").canonicalFile.exists()  //Use canonical so symbolic links can be followed
        }

        KATfilesToCheck.each { file ->
                katGivenMap[file] = new File("${params.kat}/$file").canonicalFile.exists()
        }
        if (["selfkat-stats.tsv", "selfkat-counts.cvg"].every { katGivenMap["${assembly_basename}.$it"] }) {
                println('Correct selfkat files given')
                selfkat_ch = Channel.fromPath("$params.kat/${assembly_basename}.selfkat*").collect()
        } else if (["selfkat-stats.tsv", "selfkat-counts.cvg"].every { katExistsMap["${assembly_basename}.$it"] }) {
                println('Self kat files found')
                selfkat_ch = Channel.fromPath("$params.projectDir/kat/${assembly_basename}.selfkat*").collect()

        } else {
                println('All selfkat files not found or given, will create')
                selfkat_ch = selfkat('EMPTY', assembly_ch, basename_ch)
        }

        if (["kat-stats.tsv", "kat-counts.cvg"].every { katGivenMap["${assembly_basename}.$it"] }) {
                println('Correct kat files given')
                kat_ch = Channel.fromPath("$params.kat/${assembly_basename}.kat*").collect()
        } else if (["kat-stats.tsv", "kat-counts.cvg"].every { katExistsMap["${assembly_basename}.$it"] }) {
                println('Kat files found')
                kat_ch = Channel.fromPath("$params.projectDir/kat/${assembly_basename}.kat*").collect()
        } else {
                println('All kat files not found or given, will create')
                kat_ch = KAT('EMPTY', assembly_ch, reads_ch, basename_ch)
        }

        //Check for purge_haplotig input and run if not given  (user should give folder with purge hap analysis in it)
        def PURGEfilesToCheck = [
                "${assembly_basename}.purge.reassignments.tsv",
                "${assembly_basename}.purge.coverage_stats.csv",
                "${assembly_basename}.gencov"
        ]

        def purgeExistsMap = [:]
        def purgeGivenMap = [:]

        PURGEfilesToCheck.each { file ->
                purgeExistsMap[file] = new File("${params.projectDir}/purge_haplotigs/$file").canonicalFile.exists()
        }

        PURGEfilesToCheck.each { file ->
                purgeGivenMap[file] = new File("${params.purge_hap}/$file").canonicalFile.exists()
        }


        if (PURGEfilesToCheck.every { file -> purgeGivenMap[file] }) {
                println('Purge haplotig files given')
                purgehap_ch = Channel.fromPath("$params.purge_hap/${assembly_basename}*").collect()
        } else if (PURGEfilesToCheck.every { file -> purgeExistsMap[file] }) {
                println('Purge haplotig files found in directory')
                purgehap_ch = Channel.fromPath("$params.projectDir/purge_haplotigs/${assembly_basename}*").collect()
        } else {
                println "Some PURGE files are missing or not given, let's make them"
                purgehap_ch = PURGE_HAPLOTIGS('EMPTY', assembly_ch, sortedbam_ch, scdepth_ch, basename_ch)
        }

        //Run diploidocus
        diploidocus(reads_ch, assembly_ch, scdepth_ch, sortedbam_ch, busco_ch, purgehap_ch, kat_ch, selfkat_ch, basename_ch, factorsToCheck)

        }


process readpartition {
memory { 2.GB * task.attempt }
cpus 4
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
time '1h'
        input:
        path bam
        val basename
        val assembly_type

        output:
        path '*.fastq'

        script:
        """
        bedtools bamtofastq -i $bam -fq ${basename}.${assembly_type}.fastq
        """
}


process diploidocus {
publishDir "$params.projectDir/diploidocus", mode:'symlink'
memory { 15.GB * task.attempt }
time '1h'
cpus 8
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
        input:
        path reads
        path assembly
        path scdepth_file
        path bam
        path busco
        path purgehap
        path kat
        path selfkat
        val basename
        tuple val(factor1), val(factor2)

        output:
        path '*'

        script:
        """
        scdepth=\$(cat ${scdepth_file})
        rounded_scdepth=\$(printf "%.0f" \${scdepth})


        if [[ $basename == *"$factor1"* || $basename == *"$factor2"* ]]; then
                runmode=\$"diploidocus"
        else
                runmode=\$"dicycle"
        fi


        conda run --no-capture-output -n myenv python /app/diploidocus/code/diploidocus.py seqin=$assembly \
        basefile=${basename} scdepth=\$rounded_scdepth bam=$bam reads=$reads readtype=hifi \
        busco=$busco runmode=\$runmode
        """
}
process PURGE_HAPLOTIGS {
debug true
publishDir "$params.projectDir/purge_haplotigs", mode:'symlink'
memory { 75.GB * task.attempt }
time { 2.hour * task.attempt }
cpus 16
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
        input:
        val x
        path assembly
        path bam
        path scdepth_file
        val basename

        when:
        x == 'EMPTY'

        output:
        path '*'

        script:
        """
        scdepth=\$(cat ${scdepth_file})
        rounded_scdepth=\$(printf "%.0f" \${scdepth})

        low=\$((\${rounded_scdepth} / 4))
        high=\$((\${rounded_scdepth} * 2))

        purge_haplotigs hist -b $bam -g $assembly -t $task.cpus
        mv *.gencov ${bam}.gencov
        purge_haplotigs cov -i ${bam}.gencov -l \$low -m \$rounded_scdepth -h \$high -o ${basename}.purge.coverage_stats.csv
        purge_haplotigs purge -o ${basename}.purge -g $assembly -c ${basename}.purge.coverage_stats.csv -t $task.cpus

        """
}

process selfkat {
debug true
publishDir "$params.projectDir/kat", mode:'symlink'
memory { 62.GB * task.attempt }
cpus 32
time { 2.hour * task.attempt }
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
        input:
        val x
        path assembly
        val basename

        output:
        path '*'

        when:
        x == 'EMPTY'

        script:
        """
        kat sect -t $task.cpus -o ${basename}.selfkat $assembly $assembly
        """
}

process KAT {
debug true
tag "${assembly}"
publishDir "$params.projectDir/kat", mode:'symlink'
cpus 32
memory { 80.GB * task.attempt }
time { 4.hour * task.attempt }
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
        input:
        val x
        path assembly
        path reads
        val basename

        output:
        path '*'

        when:
        x == 'EMPTY'

        script:
        """
        kat sect -t $task.cpus -o ${basename}.kat $assembly $reads
        """
}
process depthsizer {
debug true
publishDir "$params.projectDir/depthsizer", mode:'symlink'
memory { 23.GB * task.attempt }
time '4h'
cpus 4
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
//For some reason nextflow thinks normal depthsizer output is errerous so had to force exit code as 0
        input:
        val x
        path assembly
        path bam
        path busco
        path reads
        val basename

        output:
        path '*.fastmp.scdepth'

        when:
        x == 'EMPTY'

        script:
        """
        python /opt/depthsizer_code/scripts/depthsizer.py seqin=$assembly bam=$bam reads=$reads busco=$busco i=-1 readtype=hifi basefile=${basename}.depth  > /dev/null 2>&1 || (exit 0)

        if [ -e *.fastmp.scdepth ]; then
                exit 0
        else
                exit 1
        fi
        """
}
process haplotypemapping {
debug true
publishDir "$params.projectDir/bam", mode:'symlink'
cpus 32
memory { 60.GB * task.attempt }
time { 5.hour * task.attempt }
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
        input:
        path hap1
        path hap2
        path reads
        val basename

        output:
        path "*.diploid.sam"

        script:
        """
        cat $hap1 $hap2 > ${basename}.diploid.fa
        minimap2 --secondary=no -L -ax map-hifi ${basename}.diploid.fa $reads > ${basename}.diploid.sam
        """
}

process haplotypesplit {
debug true
publishDir "$params.projectDir/bam", mode:'symlink'
cpus 4
memory { 10.GB * task.attempt }
time { 2.hour * task.attempt }
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
        input:
        path sam
        path hap1
        path hap2
        val basename
        val assembly_type
        tuple val(factor1), val(factor2)

        output:
        path "*.bam"
        path "*${assembly_type}*.bam", emit: assembly_file

        script:
        """
        awk 'sub(/^>/, "")' $hap1 > ${factor1}_sqlines.txt
        awk 'sub(/^>/, "")' $hap2 > ${factor2}_sqlines.txt
        samtools sort $sam -o ${basename}.sorted.diploid.bam
        samtools index ${basename}.sorted.diploid.bam

        #This filters out any of the supplementry mapping from the other file which mixes haplotype, samtools filtering by flags doesn't work 
        samtools view -h ${basename}.sorted.diploid.bam \$(cat ${factor1}_sqlines.txt) | grep -v -f ${factor2}_sqlines.txt | samtools view -b > ${basename}.${factor1}.sorted.bam
        samtools view -h ${basename}.sorted.diploid.bam \$(cat ${factor2}_sqlines.txt) | grep -v -f ${factor1}_sqlines.txt | samtools view -b > ${basename}.${factor2}.sorted.bam
        """
}

process busco {
debug true
publishDir "$params.projectDir/busco", mode:'symlink'
tag "BUSCO on ${assembly}"
cpus 32
memory { 100.GB * task.attempt }
time '10h'
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
        input:
        val x
        path assembly
        val basename

        output:
        path '**full_table.tsv'

        when:
        x == 'EMPTY'

        script:
        """
        busco -i $assembly -o ${basename} -m genome -l $params.lineage -c $task.cpus
        """
}
process mapping {
debug true
memory { 40.GB * task.attempt }
time '2h'
cpus 16
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
        input:
        path assembly
        path reads
        val basename

        output:
        path "*.sam"

        script:
        """
        echo mapping for $assembly
        minimap2 --secondary=no -L -ax map-hifi $assembly $reads > ${basename}.sam
        """
}

process sam2bam {
publishDir "$params.projectDir/bam", mode:'symlink'
memory { 5.GB * task.attempt }
time '2h'
cpus 2
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
        input:
        path sam
        val basename

        output:
        path '*sorted.bam'

        script:
        """
        samtools view -bS $sam > ${basename}.bam
        samtools sort -@ 16 -o ${basename}.sorted.bam ${basename}.bam
        samtools index -@ 16 ${basename}.sorted.bam
        """
}
process bam_sortindex {
debug true
memory { 27.GB * task.attempt }
time '2h'
cpus 16
publishDir "$params.projectDir/bam", mode:'symlink'
errorStrategy { (task.attempt <= 3) ? 'retry' : 'finish' }
        input:
        path bam
        val basename
        val assembly_type

        output:
        path '*sort*.bam'

        script:
        """
        echo 'Sorting BAM file'
        threads=\$(($task.cpus * 2))

        if [ -z $assembly_type ]; then
                output_basename=\$"${basename}"
        else
                output_basename=\$"${basename}.${assembly_type}"
        fi

        samtools sort -@ \$threads -o \${output_basename}.sorted.bam "$bam"
        samtools index \$output_file
        """
}
