#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

for f in $(find . -xtype f -name "*.vcf.gz" | sort -V);do zcat $f | grep -v ^# | md5sum