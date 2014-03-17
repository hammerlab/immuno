#!/bin/sh

wget http://sourceforge.net/projects/snpeff/files/snpEff_v3_0_core.zip
unzip snpEff_v3_0_core.zip
mv snpEff_3_0 snpEff 
cd snpEff
java -jar snpEff.jar download GRCh37.66
