version 1.0

struct GenomeResources {
    String refFai
    String refFasta
    String refDict
    String modules
    String mergeVcfModules
}

workflow vardict {
    input {
        File tumor_bam
        File normal_bam
        String tumor_sample_name
        String normal_sample_name
        String bed_file
        String reference
    }

    parameter_meta {
        tumor_bam: "tumor_bam file for analysis sample"
        normal_bam: "normal_bam file for analysis sample"
        tumor_sample_name:"Sample name for the tumor bam"
        normal_sample_name: "Sample name for the normal bam"
        bed_file: "BED files for specifying regions of interest"
        reference: "the reference genome for input sample"
    }

    Map[String, GenomeResources] resources = {
        "hg19": {
            "refFai" : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa.fai",
            "refFasta" : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa",
            "refDict" : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.dict",
            "modules" : "hg19/p13 rstats/4.2 java/9 perl/5.30 vardict/1.8.3 bcftools/1.9",
            "mergeVcfModules": "picard/2.19.2 hg19/p13"
        },
        "hg38": {
            "refFai" : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa.fai",
            "refFasta" : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa",
            "refDict" : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.dict",
            "modules" : "hg38/p12 rstats/4.2 java/9 perl/5.30 vardict/1.8.3 bcftools/1.9",
            "mergeVcfModules": "picard/2.19.2  hg38/p12"
        }
    }

    call splitBedByChromosome {
        input:
        bed_file = bed_file
    }

    # run vardict

    scatter (i in range(length(splitBedByChromosome.bed_files))) {
        call runVardict { 
            input: 
                tumor_bam = tumor_bam,
                normal_bam = normal_bam,
                tumor_sample_name = tumor_sample_name,
                normal_sample_name = normal_sample_name,
                bed_file = splitBedByChromosome.bed_files[i],
                memory_coefficient = splitBedByChromosome.memory_coefficients[i],
                modules = resources [ reference ].modules,
                refFai = resources[reference].refFai,
                refFasta = resources[reference].refFasta,
        }
    }
    Array[File] vardictVcfs = runVardict.vcf_file
    Array[File] vardictVcfIndexes = runVardict.vcf_index

    call mergeVCFs {
    input:
      vcfs = vardictVcfs,
      vcfIndexes = vardictVcfIndexes,
      refDict = resources[reference].refDict,
      modules = resources[reference].mergeVcfModules
    }

    meta {
        author: "Gavin Peng"
        email: "gpeng@oicr.on.ca"
        description: "VarDict is an ultra sensitive variant caller for both single and paired sample variant calling from BAM files. VarDict implements several novel features such as amplicon bias aware variant calling from targeted sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability than many Java based variant callers."
        dependencies: [
            {
                name: "vardict/1.8.3",
                url: "https://github.com/pachterlab/vardict"
            },
            {
                name: "java",
                url: "https://www.java.com/en/"
            }
        ]
        output_meta: {
            vardict_vcf: {
                description: "VCF file for variant calling from vardict",
                vidarr_label: "vardict_vcf"
            }
        }
    }
    
    output {
        File vardictVcf = mergeVCFs.mergedVcf
        File vardictVcfIndex = mergeVCFs.mergedVcfIdx
    }

}

task splitBedByChromosome {
    input {
        File bed_file  
        Int memory = 1
        Int timeout = 1
    }
    parameter_meta {
        bed_file: "Input BED file to split"
        memory: "Memory allocated for task in GB"
        timeout: "Timeout in hours"
    }

    command <<<
        set -euo pipefail
        
        mkdir split_beds
        CHROMS=($(seq 1 22) X Y)
        
        for chr in "${CHROMS[@]}"; do
            grep -E "^(chr)?${chr}[[:space:]]" ~{bed_file} > split_beds/chr${chr}.bed || true
            if [ -s split_beds/chr${chr}.bed ]; then
                echo "split_beds/chr${chr}.bed" >> split_beds.list
                # Calculate range size for this bed
                range_size=$(awk '{sum += ($3 - $2)} END {print sum}' split_beds/chr${chr}.bed)
                echo "${range_size}" >> range_sizes.txt
            fi
        done
        max_size=$(sort -n range_sizes.txt | tail -n1)
        while IFS= read -r size; do
            awk -v size="$size" -v max="$max_size" 'BEGIN {printf "%.1f\n", size/max}' 
        done < range_sizes.txt > memory_coefficients.txt
    >>>

    output {
        Array[File] bed_files = read_lines("split_beds.list")
        Array[Float] memory_coefficients = read_lines("memory_coefficients.txt")
    }

    runtime {
        memory:  "~{memory} GB"
        timeout: "~{timeout}"
    }
}   

# ==========================
#  configure and run vardict
# ==========================
task runVardict {
    input {
        File tumor_bam
        File normal_bam
        String tumor_sample_name
        String normal_sample_name
        String refFasta
        String refFai
        String AF_THR = 0.01
        String MAP_QUAL = 10
        String READ_POSTION_FILTER = 5
        String modules
        String bed_file
        Int timeout = 96
        Int memory = 32
        Int minMemory = 24
        Float memory_coefficient
    }
    parameter_meta {
        tumor_bam: "tumor_bam file for analysis sample"
        normal_bam: "normal_bam file for analysis sample"
        tumor_sample_name:"Sample name for the tumor bam"
        normal_sample_name: "Sample name for the normal bam"
        refFai: "Reference fasta fai index"
        refFasta: "The reference fasta"
        AF_THR: "The threshold for allele frequency, default: 0.01 or 1%"
        MAP_QUAL: " Mapping quality. If set, reads with mapping quality less than the number will be filtered and ignored"
        READ_POSTION_FILTER: "The read position filter. If the mean variants position is less that specified, it is considered false positive. Default: 5"
        bed_file: "BED files for specifying regions of interest"
        modules: "Names and versions of modules"
        timeout: "Timeout in hours, needed to override imposed limits"
    } 
    Int allocatedMemory = if minMemory > round(memory * memory_coefficient) then minMemory else round(memory * memory_coefficient)

    command <<<
        set -euo pipefail
        cp ~{refFai} .
        
        export JAVA_OPTS="-Xmx$(echo "scale=0; ~{memory} * 0.9 / 1" | bc)G"
        $VARDICT_ROOT/bin/VarDict \
            -G ~{refFasta} \
            -f ~{AF_THR} \
            -N ~{tumor_sample_name} \
            -b "~{tumor_bam}|~{normal_bam}" \
            -Q ~{MAP_QUAL} \
            -P ~{READ_POSTION_FILTER} \
            -c 1 -S 2 -E 3 -g 4 \
            ~{bed_file} | \
            $RSTATS_ROOT/bin/Rscript $VARDICT_ROOT/bin/testsomatic.R | \
            $PERL_ROOT/bin/perl $VARDICT_ROOT/bin/var2vcf_paired.pl \
            -N "~{tumor_sample_name}|~{normal_sample_name}" \
            -f ~{AF_THR} | gzip  > vardict.vcf.gz

            # the vardict generated vcf header missing contig name, need extract contig lines from refFai
            bcftools view -h vardict.vcf.gz > header.txt     
            while read -r name length rest; do 
                echo "##contig=<ID=$name,length=$length>" 
            done < ~{refFai} >> header.txt

            bcftools reheader -h header.txt -o ~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf.gz vardict.vcf.gz
            bcftools index --tbi ~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf.gz
    >>>

    runtime {
        modules: "~{modules}"
        timeout: "~{timeout}"
        memory:  "~{allocatedMemory} GB"
    }

    output {
        File vcf_file = "~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf.gz"
        File vcf_index = "~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf.gz.tbi"
    }
}

task mergeVCFs {
  input {
    String modules
    Array[File] vcfs
    Array[File] vcfIndexes 
    String refDict
    Int memory = 4
    Int timeout = 12
  }

  parameter_meta {
    modules: "Environment module names and version to load (space separated) before command execution"
    vcfs: "Vcf's from scatter to merge together"
    memory: "Memory allocated for job"
    timeout: "Hours before task timeout"
  }

  meta {
    output_meta: {
      mergedVcf: "Merged vcf, unfiltered.",
      mergedVcfIdx: "Merged vcf index, unfiltered."
    }
  }

  String outputName = basename(vcfs[0])

  command <<<
    set -euo pipefail

    java "-Xmx~{memory-3}g" -jar $PICARD_ROOT/picard.jar MergeVcfs \
    I=~{sep=" I=" vcfs} \
    O=~{outputName} \
    SEQUENCE_DICTIONARY=~{refDict}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedVcf = "~{outputName}"
    File mergedVcfIdx = "~{outputName}.tbi"
  }
}