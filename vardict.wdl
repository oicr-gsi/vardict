version 1.0


workflow vardict {
    input {
        File tumor_bam
        File normal_bam
        String tumor_sample_name
        String normal_sample_name
    }

    parameter_meta {
        tumor_bam: "Input fastqR1 file for analysis sample"
        normal_bam: "Input fastqR2 file for analysis sample"
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
        File vcf_out = runVardict.vcf_file
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
        String refFasta = "$HG19_ROOT/hg19_random.fa"
        String AF_THR = 0.01
        String modules = "samtools/1.16.1 rstats/4.2 java/9 perl/5.30 vardict/1.8.3 hg19/p13"
        String bed_file = "/.mounts/labs/gsi/testdata/mutect2/input_data/PCSI0022.val.bed"
        Int timeout = 48
        Int jobMemory = 24
    }

    command <<<
        $VARDICT_ROOT/bin/VarDict \
            -G ~{refFasta} \
            -f ~{AF_THR} \
            -N "~{tumor_sample_name}" \
            -b ~{tumor_bam} | ~{normal_bam} \
            -Q 10 \
            -c 1 -S 2 -E 3 -g 4 \
            ~{bed_file} | \
            $VARDICT_ROOT/bin/testsomatic.R | \
            $VARDICT_ROOT/bin/var2vcf_paired.pl -N "~{tumor_sample_name}| ~{normal_sample_name}" -f 0.03 > ~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf
  
    >>>

    runtime {
        memory:  "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }

    output {
        File vcf_file = "~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf"

    }
}