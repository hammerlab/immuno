#!/bin/sh

java -Xmx4g -jar snpEff.jar eff -v GRCh37.66 ${1} > ${1}.eff.vcf
