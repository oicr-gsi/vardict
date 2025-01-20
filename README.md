# vardict

VarDict is an ultra sensitive variant caller for both single and paired sample variant calling from BAM files. VarDict implements several novel features such as amplicon bias aware variant calling from targeted sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability than many Java based variant callers.

## Overview

## Dependencies

* [vardict 1.8.3](https://github.com/pachterlab/vardict)
* [java](https://www.java.com/en/)


## Usage

### Cromwell
```
java -jar cromwell.jar run vardict.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumor_bam`|File|tumor_bam file for analysis sample
`normal_bam`|File|normal_bam file for analysis sample
`tumor_sample_name`|String|Sample name for the tumor bam
`normal_sample_name`|String|Sample name for the normal bam
`bed_file`|String|BED files for specifying regions of interest
`reference`|String|the reference genome for input sample


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`splitBedByChromosome.memory`|Int|1|Memory allocated for task in GB
`splitBedByChromosome.timeout`|Int|1|Timeout in hours
`runVardict.AF_THR`|String|0.01|The threshold for allele frequency, default: 0.01 or 1%
`runVardict.MAP_QUAL`|String|10| Mapping quality. If set, reads with mapping quality less than the number will be filtered and ignored
`runVardict.READ_POSTION_FILTER`|String|5|The read position filter. If the mean variants position is less that specified, it is considered false positive. Default: 5
`runVardict.timeout`|Int|96|Timeout in hours, needed to override imposed limits
`runVardict.memory`|Int|32|base memory for this job
`runVardict.minMemory`|Int|24|The minimum value for allocated memory
`mergeVCFs.memory`|Int|4|Memory allocated for job
`mergeVCFs.timeout`|Int|12|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`vardictVcf`|File|Merged vcf, unfiltered.|
`vardictVcfIndex`|File|VCF index for variant calling from vardict|vidarr_label: vardcit_index


## Commands
 This section lists command(s) run by vardict workflow
 
 * Running vardict
 
 ```
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
```
 ```
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
```
 ```
     set -euo pipefail
 
     java "-Xmx~{memory-3}g" -jar $PICARD_ROOT/picard.jar MergeVcfs \
     I=~{sep=" I=" vcfs} \
     O=~{outputName} \
     SEQUENCE_DICTIONARY=~{refDict}
```

 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
