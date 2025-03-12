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
        call normalizeVardictVcf {
            input:
                vcf_file = runVardict.vcf_file
        }
    }
    Array[File] vardictVcfs = normalizeVardictVcf.normalized_vcf
    Array[File] vardictVcfIndexes = normalizeVardictVcf.vcf_index

    call mergeVcfs {
    input:
      vcfs = vardictVcfs,
      vcfIndexes = vardictVcfIndexes,
      tumor_sample_name = tumor_sample_name,
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
            },
            vardictVcfIndex: {
                description: "VCF index for variant calling from vardict",
                vidarr_label: "vardcit_index"
            }
        }
    }
    
    output {
        File vardictVcf = mergeVcfs.mergedVcf
        File vardictVcfIndex = mergeVcfs.mergedVcfIdx
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
        min_coff=0.8
        p=0.8
        max_size=$(sort -n range_sizes.txt | tail -n1)

        awk -v max="$max_size" -v min_coff="$min_coff" -v p="$p" '
        {
            size = $1
            ratio = size / max
            coefficient = min_coff + (1 - min_coff) * (ratio ^ p)
            printf "%.2f\n", coefficient
        }' range_sizes.txt > memory_coefficients.txt
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
        String READ_POSITION_FILTER = 5
        String modules
        String bed_file
        Int timeout = 120
        Int memory = 32
        Int minMemory = 24
        Int numThreads = 8
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
        READ_POSITION_FILTER : "The read position filter. If the mean variants position is less that specified, it is considered false positive. Default: 5"
        bed_file: "BED files for specifying regions of interest"
        modules: "Names and versions of modules"
        timeout: "Timeout in hours, needed to override imposed limits"
        memory: "base memory for this job"
        memory_coefficient: "coefficient for calculating allocated memory"
        minMemory: "The minimum value for allocated memory"
        numThreads: "Number of threads"
    } 
    Int allocatedMemory = if minMemory > round(memory * memory_coefficient) then minMemory else round(memory * memory_coefficient)

    command <<<
        set -euo pipefail
        cp ~{refFai} .
        
        export JAVA_OPTS="-Xmx$(echo "scale=0; ~{allocatedMemory} * 0.8 / 1" | bc)G"
        $VARDICT_ROOT/bin/VarDict \
            -th ~{numThreads} \
            -G ~{refFasta} \
            -f ~{AF_THR} \
            -N ~{tumor_sample_name} \
            -b "~{tumor_bam}|~{normal_bam}" \
            -Q ~{MAP_QUAL} \
            -P ~{READ_POSITION_FILTER} \
            -c 1 -S 2 -E 3 -g 4 \
             ~{bed_file} | \
            $RSTATS_ROOT/bin/Rscript $VARDICT_ROOT/bin/testsomatic.R | \
            $PERL_ROOT/bin/perl $VARDICT_ROOT/bin/var2vcf_paired.pl \
            -N "~{tumor_sample_name}|~{normal_sample_name}" \
            -f ~{AF_THR} | bgzip  > vardict.vcf.gz

            # the vardict generated vcf header missing contig name, need extract contig lines from refFai
            bcftools view -h vardict.vcf.gz > header.txt     
            while read -r name length rest; do 
                echo "##contig=<ID=$name,length=$length>" 
            done < ~{refFai} >> header.txt

            bcftools reheader -h header.txt -o ~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf.gz vardict.vcf.gz
            
    >>>

    runtime {
        modules: "~{modules}"
        timeout: "~{timeout}"
        memory:  "~{allocatedMemory} GB"
    }

    output {
        File vcf_file = "~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf.gz"
    }
}

task normalizeVardictVcf {
  input {
    File vcf_file
    Int memory = 4
    String modules = "python/3.6 htslib/1.9"
    Int timeout = 12
  }
  parameter_meta {
    modules: "Environment module names and version to load (space separated) before command execution"
    vcf_file: "Vcf that's need be normalized because vardict generates non-standard notation vcf"
    memory: "Memory allocated for job"
    timeout: "Hours before task timeout"
  }

  String output_name = basename(vcf_file, ".vcf.gz") + ".normalized.vcf"

command <<<
    set -euo pipefail
    # normalize vcf file from vardict since vardict may generate vcf with non-standard notation
    
    python3 <<CODE
    import sys
    import re
    import gzip
    import os

    input_file = "~{vcf_file}"
    output_file = "~{output_name}"

    is_gzipped = input_file.endswith('.gz')
    if is_gzipped:
        in_fh = gzip.open(input_file, 'rt')
    else:
        in_fh = open(input_file, 'r')
    out_fh = open(output_file, 'w')
    dup_pattern = re.compile(r'<dup-(\d+)>')

    header_added = False
    for line in in_fh:
        # Add custom field definitions to the header
        if line.startswith('#') and not header_added and "##INFO=<ID=STATUS" in line:
            # Insert our custom header lines before the first INFO line
            header_lines = [
                '##INFO=<ID=DUP_SIZE,Number=1,Type=Integer,Description="Size of the duplication in base pairs">\n',
                '##INFO=<ID=IS_DUP,Number=0,Type=Flag,Description="Variant was originally annotated as a duplication">\n',
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n',
                '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">\n',
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n'
            ]
            for header_line in header_lines:
                out_fh.write(header_line)
            header_added = True
        
        # Insert headers at the end of the header section if not added earlier
        if line.startswith('#CHROM') and not header_added:
            header_lines = [
                '##INFO=<ID=DUP_SIZE,Number=1,Type=Integer,Description="Size of the duplication in base pairs">\n',
                '##INFO=<ID=IS_DUP,Number=0,Type=Flag,Description="Variant was originally annotated as a duplication">\n',
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n',
                '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">\n',
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n'
            ]
            for header_line in header_lines:
                out_fh.write(header_line)
            header_added = True
        
        # Output header lines unchanged
        if line.startswith('#'):
            out_fh.write(line)
            continue
        
        # Parse fields
        fields = line.strip().split('\t')
        
        # Look for <dup-XX> pattern in ALT field (field 5)
        if len(fields) >= 5 and len(fields[4]) > 0:
            match = dup_pattern.search(fields[4])
            if match:
                # Remove the <dup-XX> notation
                fields[4] = dup_pattern.sub("", fields[4])
                
                # Update INFO field (field 8)
                if len(fields) >= 8:
                    dup_size = match.group(1)
                    pos = int(fields[1])
                    end = pos + int(dup_size)
                    
                    # Add structural variant information to INFO field
                    info = fields[7]
                    
                    # If it was already marked as an insertion, add duplication info
                    if "TYPE=Insertion" in info:
                        # Convert to a more standard SVTYPE format
                        info = info.replace("TYPE=Insertion", "SVTYPE=INS")
                        # Add SVLEN and END for downstream compatibility
                        info = info + ";DUP_SIZE=" + dup_size + ";SVLEN=" + dup_size + ";END=" + str(end)
                        # Add flag indicating this was originally a duplication
                        info = info + ";IS_DUP=1"
                    
                    fields[7] = info
        
        # Write the modified or unmodified line
        out_fh.write('\t'.join(fields) + '\n')

    in_fh.close()
    out_fh.close()
    CODE
    bgzip ~{output_name}
    tabix -p vcf ~{output_name}.gz
>>>

  output {
    File normalized_vcf = "~{output_name}.gz"
    File vcf_index ="~{output_name}.gz.tbi"
  }

  runtime {
    memory: "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
}

task mergeVcfs {
  input {
    String modules
    String tumor_sample_name
    Array[File] vcfs
    Array[File] vcfIndexes 
    String refDict
    Int memory = 4
    Int timeout = 12
  }

  parameter_meta {
    modules: "Environment module names and version to load (space separated) before command execution"
    vcfs: "Vcf's from scatter to merge together"
    vcfIndexes: "The indices for the input vcfs"
    memory: "Memory allocated for job"
    timeout: "Hours before task timeout"
    refDict: "reference sequence dictionary"
  }

  meta {
    output_meta: {
      mergedVcf: "Merged vcf, unfiltered.",
      mergedVcfIdx: "Merged vcf index, unfiltered."
    }
  }

   command <<<
    set -euo pipefail

    java "-Xmx~{memory-3}g" -jar $PICARD_ROOT/picard.jar MergeVcfs \
    ~{sep=" " prefix("I=", vcfs)} \
    O=~{tumor_sample_name}.vardict.vcf.gz \
    SEQUENCE_DICTIONARY=~{refDict}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedVcf = "~{tumor_sample_name}.vardict.vcf.gz"
    File mergedVcfIdx = "~{tumor_sample_name}.vardict.vcf.gz.tbi"
  }
}