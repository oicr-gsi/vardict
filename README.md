# vardict

VarDict is an ultra sensitive variant caller for both single and paired sample variant calling from BAM files. VarDict implements several novel features such as amplicon bias aware variant calling from targeted sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability than many Java based variant callers.

## Overview

## Dependencies

* [vardict 1.8.3](https://github.com/pachterlab/vardict)
* [tabix 0.2.6](http://www.htslib.org)
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
`tumor_bai`|File|index for tumor bam file
`normal_bam`|File|normal_bam file for analysis sample
`normal_bai`|File|index for normal bam file
`tumor_sample_name`|String|Sample name for the tumor bam
`normal_sample_name`|String|Sample name for the normal bam
`bed_file`|String|BED files for specifying regions of interest
`reference`|String|the reference genome for input sample


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`splitBedByChromosome.memory`|Int|1|Memory allocated for task in GB
`splitBedByChromosome.timeout`|Int|1|Timeout in hours
`runVardict.AF_THR`|String|0.01|The threshold for allele frequency, default: 0.01 or 1%
`runVardict.MAP_QUAL`|String|10| Mapping quality. If set, reads with mapping quality less than the number will be filtered and ignored
`runVardict.READ_POSITION_FILTER`|String|5|The read position filter. If the mean variants position is less that specified, it is considered false positive. Default: 5
`runVardict.timeout`|Int|120|Timeout in hours, needed to override imposed limits
`runVardict.memory`|Int|32|base memory for this job
`runVardict.minMemory`|Int|24|The minimum value for allocated memory
`runVardict.numThreads`|Int|8|Number of threads
`mergeVcfs.memory`|Int|4|Memory allocated for job
`mergeVcfs.timeout`|Int|12|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`vardictVcf`|File|Merged vcf, unfiltered.|
`vardictVcfIndex`|File|VCF index for variant calling from vardict|vidarr_label: vardcit_index


## Commands
 This section lists command(s) run by vardict workflow
 
 * Running vardict
 
 
 ### Calculate memory allocation coefficients
 
 Depending on the size, each chromosome-specific shard gets a specific amount of RAM
 
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
 ```
 
 ### Run vardict
 
 This is a task which runs in a scatter, one shard per chromosome
 
 ```
         set -euo pipefail
         cp ~{refFai} .
         
         export JAVA_OPTS="-Xmx$(echo "scale=0; ~{allocatedMemory} * 0.8 / 1" | bc)G"
         $VARDICT_ROOT/bin/VarDict \
             -th ~{numThreads} \
             -G ~{refFasta} \
             -f ~{AF_THR} \
             -N ~{tumor_sample_name} \
             -b "~{normal_bam}|~{tumor_bam}" \
             -Q ~{MAP_QUAL} \
             --nosv \
             -P ~{READ_POSITION_FILTER} \
             -c 1 -S 2 -E 3 -g 4 \
              ~{bed_file} | \
             $RSTATS_ROOT/bin/Rscript $VARDICT_ROOT/bin/testsomatic.R | \
             $PERL_ROOT/bin/perl $VARDICT_ROOT/bin/var2vcf_paired.pl \
             -N "~{normal_sample_name}|~{tumor_sample_name}" \
             -f ~{AF_THR} | bgzip  > vardict.vcf.gz
 
             # the vardict generated vcf header missing contig name, need extract contig lines from refFai
             bcftools view -h vardict.vcf.gz > header.txt     
             while read -r name length rest; do 
                 echo "##contig=<ID=$name,length=$length>" 
             done < ~{refFai} >> header.txt
 
             bcftools reheader -h header.txt -o ~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf.gz vardict.vcf.gz
             tabix -p vcf ~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf.gz
             
 ```
 ### Merge chromosome-specific vcf files
 
 Using vcftools we merge chromosome-specific chunks into one vcf. Additional post-processing required
 to remove extra header lines from chunk-specific vcf files
 
 ```
     set -euo pipefail
     $VCFTOOLS_ROOT/bin/vcf-concat ~{sep=" " vcfs} > temp.vcf
     grep ^# test.vcf | perl -ne 'BEGIN{$end=0}{print if !$end;if(/CHROM/){$end = 1;}}' > ~{tumor_sample_name}.vardict.vcf
     grep -v ^# temp.vcf >> ~{tumor_sample_name}.vardict.vcf
     bgzip ~{tumor_sample_name}.vardict.vcf
     tabix -p vcf ~{tumor_sample_name}.vardict.vcf.gz
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
