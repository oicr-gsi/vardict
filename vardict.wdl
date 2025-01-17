version 1.0

struct GenomeResources {
    String refFai
    String refFasta
    String refDict
    String modules
    String splitRegionModule
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
        Int mem_coefficient
    }

    parameter_meta {
        tumor_bam: "tumor_bam file for analysis sample"
        normal_bam: "normal_bam file for analysis sample"
        tumor_sample_name:"Sample name for the tumor bam"
        normal_sample_name: "Sample name for the normal bam"
        bed_file: "BED files for specifying regions of interest"
        reference: "the reference genome for input sample"
        intervalsToParallelizeBy: "Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4)"
    }

    Map[String, GenomeResources] resources = {
        "hg19": {
            "refFai" : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa.fai",
            "refFasta" : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa",
            "refDict" : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.dict",
            "modules" : "hg19/p13 rstats/4.2 java/9 perl/5.30 vardict/1.8.3",
            "splitRegionModule": "hg19/p13",
            "mergeVcfModules": "gatk hg19/p13"
        },
        "hg38": {
            "refFai" : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa.fai",
            "refFasta" : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa",
            "refDict" : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.dict",
            "modules" : "hg38/p12 rstats/4.2 java/9 perl/5.30 vardict/1.8.3",
            "splitRegionModule": "hg38/p12",
            "mergeVcfModules": "gatk hg38/p12"
        }
    }

    call splitRegion {
        input:
        refFai = resources[reference].refFai,
        mem_coefficient = mem_coefficient,
        modules = resources [ reference ].splitRegionModule
    }

    # run vardict
    scatter (region in splitRegion.regions) {
        call runVardict { 
            input: 
                tumor_bam = tumor_bam,
                normal_bam = normal_bam,
                tumor_sample_name = tumor_sample_name,
                normal_sample_name = normal_sample_name,
                bed_file = bed_file,
                modules = resources [ reference ].modules,
                refFai = resources[reference].refFai,
                refFasta = resources[reference].refFasta,
                region = region,
                memory = splitRegion.region_memory_map[region]
        }
    }
    Array[File] vardictVcfs = runVardict.vcf_file

    call mergeVCFs {
    input:
      vcfs = vardictVcfs,
      refDict = resources[reference].refDict,
      refFasta = resources[reference].refFasta,
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

task splitRegion {
    input {
        File refFai  
        Int mem_coefficient = 500  # Memory coefficient (GB per Gb)
        Int overhead = 4       # Overhead memory (GB)
        String modules
        Int memory = 1
        Int timeout = 1
    }

    command <<<
        set -euo pipefail
        grep -E '^chr[0-9XY]{1,2}\s' ~{refFai} > filtered_genome.fasta.fai

        awk -v coef=~{mem_coefficient} -v over=~{overhead} '{
            size_gb = $2 / 1e9
            memory = int(size_gb * coef + over + 0.5)
            region = $1 "_1_" $2
            print region "\t" memory > "region_memory_map.tsv"  # Tab-separated key-value map
            print region > "regions.txt"  # List of regions only
        }' filtered_genome.fasta.fai 
    >>>

    runtime {
        memory:  "~{memory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }
    output {
        Array[String] regions = read_lines("regions.txt")  # Array of regions only
        Map[String, Int] region_memory_map = read_map("region_memory_map.tsv")  # Map of region to memory
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
        String region
        Int memory
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

    command <<<
        set -euo pipefail
        cp ~{refFai} .
        original_region=$(echo ~{region} | sed 's/_/:/; s/_/-/')
        
        export JAVA_OPTS="-Xmx$(echo "scale=0; ~{memory} * 0.9 / 1" | bc)G"

        $VARDICT_ROOT/bin/VarDict \
            -G ~{refFasta} \
            -f ~{AF_THR} \
            -N ~{tumor_sample_name} \
            -b "~{tumor_bam}|~{normal_bam}" \
            -Q ~{MAP_QUAL} \
            -P ~{READ_POSTION_FILTER} \
            -c 1 -S 2 -E 3 -g 4 \
            -R ${original_region} \
            ~{bed_file} | \
            $RSTATS_ROOT/bin/Rscript $VARDICT_ROOT/bin/testsomatic.R | \
            $PERL_ROOT/bin/perl $VARDICT_ROOT/bin/var2vcf_paired.pl \
            -N "~{tumor_sample_name}|~{normal_sample_name}" \
            -f ~{AF_THR} | gzip  > ~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf.gz
    >>>

    runtime {
        modules: "~{modules}"
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
    }

    output {
        File vcf_file = "~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf.gz"
    }
}

task mergeVCFs {
  input {
    String modules
    Array[File] vcfs
    String refDict
    String refFasta
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

    cp ~{refDict} .

    gatk --java-options "-Xmx~{memory-3}g" MergeVcfs \
    -I ~{sep=" -I " vcfs} \
    -O ~{outputName}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedVcf = "~{outputName}"
    File mergedVcfIdx = "~{outputName}.idx"
  }
}