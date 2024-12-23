version 1.0


workflow vardict {
    input {
        File tumor_bam
        File normal_bam
        String tumor_sample_name
        String normal_sample_name
    }

    parameter_meta {
        tumor_bam: "tumor_bam file for analysis sample"
        normal_bam: "normal_bam file for analysis sample"
        tumor_sample_name:"Sample name for the tumor bam"
        normal_sample_name: "Sample name for the normal bam"
    }

    # run vardict
    call runVardict 
        { 
            input: 
                tumor_bam = tumor_bam,
                normal_bam = normal_bam,
                tumor_sample_name = tumor_sample_name,
                normal_sample_name = normal_sample_name
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
                name: "java",
                url: "https://www.java.com/en/"
            }
        ]
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
        String refFasta = "$HG38_ROOT/hg38_random.fa"
        String AF_THR = 0.01
        String MAP_QUAL = 10
        String modules = "samtools/1.16.1 rstats/4.2 java/9 perl/5.30 vardict/1.8.3 hg38/p12"
        String bed_file = "/.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/workflow/vardict/test/renamed_hg38.bed"
        Int timeout = 48
        Int jobMemory = 24
    }
    parameter_meta {
        tumor_bam: "tumor_bam file for analysis sample"
        normal_bam: "normal_bam file for analysis sample"
        tumor_sample_name:"Sample name for the tumor bam"
        normal_sample_name: "Sample name for the normal bam"
        refFasta: "The reference fasta"
        AF_THR: "The threshold for allele frequency, default: 0.01 or 1%"
        MAP_QUAL: " Mapping quality. If set, reads with mapping quality less than the number will be filtered and ignored"
        bed_file: "BED files for specifying regions of interest"
        jobMemory: "Memory in Gb for this job"
        modules: "Names and versions of modules"
        timeout: "Timeout in hours, needed to override imposed limits"
    }

    command <<<
        /.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/workflow/vardict/test/vardict.pl \
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
            -f ~{AF_THR}  > ~{tumor_sample_name}_~{normal_sample_name}_perl.vardict.vcf
  
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