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


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`runVardict.refFasta`|String|"$HG38_ROOT/hg38_random.fa"|The reference fasta
`runVardict.AF_THR`|String|0.01|The threshold for allele frequency, default: 0.01 or 1%
`runVardict.MAP_QUAL`|String|10| Mapping quality. If set, reads with mapping quality less than the number will be filtered and ignored
`runVardict.READ_POSTION_FILTER`|String|5|The read position filter. If the mean variants position is less that specified, it is considered false positive. Default: 5
`runVardict.modules`|String|"rstats/4.2 java/9 perl/5.30 vardict/1.8.3 hg38/p12"|Names and versions of modules
`runVardict.timeout`|Int|48|Timeout in hours, needed to override imposed limits
`runVardict.jobMemory`|Int|24|Memory in Gb for this job


### Outputs

Output | Type | Description
---|---|---
`vardict_vcf`|File|{'description': 'VCF file for variant calling from vardict', 'vidarr_label': 'vardict_vcf'}


## Commands
 This section lists command(s) run by WORKFLOW workflow
 
 * Running WORKFLOW
 
 === Description here ===.
 
 <<<
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
             -f ~{AF_THR}  > ~{tumor_sample_name}_~{normal_sample_name}.vardict.vcf
   
     >>>
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
