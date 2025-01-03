version 1.0

struct GenomeResources {
    String refFai
    String refFasta
    String modules
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
            "refFai" : "$HG19_ROOT/hg19_random.fa.fai",
            "refFasta" : "$HG19_ROOT/hg19_random.fa",
            "modules" : "hg19/p13 rstats/4.2 java/9 perl/5.30 vardict/1.8.3 samtools/1.16.1"
        },
        "hg38": {
            "refFai" : "$HG38_ROOT/hg38_random.fa.fai",
            "refFasta" : "$HG38_ROOT/hg38_random.fa",
            "modules" : "hg38/p12 rstats/4.2 java/9 perl/5.30 vardict/1.8.3 samtools/1.16.1"
        }
    }

    # run vardict
    call runVardict 
        { 
            input: 
                tumor_bam = tumor_bam,
                normal_bam = normal_bam,
                tumor_sample_name = tumor_sample_name,
                normal_sample_name = normal_sample_name,
                bed_file = bed_file,
                modules = resources [ reference ].modules,
                refFai = resources[reference].refFai,
                refFasta = resources[reference].refFasta
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
                name: "samtools/1.16.1",
                url: "https://www.htslib.org/"
            },
            {
                name: "perl",
                url: "https://www.perl.org/"
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
        File vardict_vcf = runVardict.vcf_file
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
        Int timeout = 48
        Int jobMemory = 48
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
        jobMemory: "Memory in Gb for this job"
        modules: "Names and versions of modules"
        timeout: "Timeout in hours, needed to override imposed limits"
    }

    command <<<
        set -euo pipefail
        $PERL_ROOT/bin/perl /.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/workflow/vardict/test/vardict.pl \
            -G ~{refFasta} \
            -f ~{AF_THR} \
            -N ~{tumor_sample_name} \
            -b "~{tumor_bam}|~{normal_bam}" \
            -c 1 -S 2 -E 3 -g 4 \
            -Q ~{MAP_QUAL} \
            ~{bed_file} | \
            $RSTATS_ROOT/bin/Rscript $VARDICT_ROOT/bin/testsomatic.R | \
            $PERL_ROOT/bin/perl $VARDICT_ROOT/bin/var2vcf_paired.pl \
            -N "~{tumor_sample_name}|~{normal_sample_name}" \
            -f ~{AF_THR} \
            -P 0.05   > ~{tumor_sample_name}_~{normal_sample_name}_perl.vardict.vcf
  
    >>>

    runtime {
        memory:  "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }

    output {
        File vcf_file = "~{tumor_sample_name}_~{normal_sample_name}_perl.vardict.vcf"

    }
}