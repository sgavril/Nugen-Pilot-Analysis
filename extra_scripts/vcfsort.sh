#!/bin/bash
# Slower, but handles streams.
[ $# -eq 0 ] && { echo "Sorts a VCF file in natural chromosome order";\
                  echo "Usage: $0 [my.vcf | my.vcf.gz]"; \
                  echo "Usage: cat [my.vcf | my.vcf.gz] | $0"; \
                  exit 1;
                 }
# cheers, @colbychiang
if zless $1 | awk '$0~"^#" { print $0; next } { print $0 | "LC_ALL=C sort -k1,1V -k2,2n" }';
then
    exit 0
else
    printf 'sort failed. Does your version of sort support the -V option?\n'
    printf 'If not, you should update sort with the latest from GNU coreutils:\n'
    printf 'git clone git://git.sv.gnu.org/coreutils'
fi
