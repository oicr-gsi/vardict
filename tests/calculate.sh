#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

find . -name '*.vcf'  | xargs md5sum 