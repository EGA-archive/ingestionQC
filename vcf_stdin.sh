#!/bin/bash

currentDir=$1; #KK
fileDir=$2; #outputDIR
fileName=$3; #filename to extract extensions/format
filePath=$4; #KK
pipeErrorFolder=$5;
KeyFile=$6;


######################################################################################
###########Paths to the tools used in this pipeline###################################
######################################################################################
vcfToolCommand="/homes/angel/vcftools_0.1.13/bin/vcftools";
vcfToolCommand="/usr/local/vcftools-0.1.13/bin/vcftools";
vcfValidatorCommand="/homes/angel/vcf_validator/vcf_validator-v0.5";
vcfValidatorCommand="/usr/local/vcf_validator";
bcfToolCommand="/homes/angel/bcftools-1.4-INST-EBI-comp-RH7.4/bin/bcftools";
bcfToolCommand="bcftools";
bgzipCommand="/homes/angel/samtoolsINST-1.3.1-EBI-comp-RH6.3/bin/bgzip";
bgzipCommand="bgzip";
#tabixCommand="/homes/angel/samtoolsINST-1.3.1-EBI-comp-RH6.3/bin/tabix";
#tabixCommand="tabix";
#plotsStatCommand="/homes/angel/bcftools-1.4-INST-EBI-comp-RH7.4/bin/plot-vcfstats"; (echo "$plotsStatCommand;$($plotsStatCommand -h)") >> $SWversions 2>&1
######################################################################################
######################################################################################

######################################################################################
############Printing the above tools version to a file################################
######################################################################################
SWversions="$fileDir/SW_versions"       ###file where to store the versions of the tools used in this pipeline
> $SWversions
(echo "$vcfToolCommand;$($vcfToolCommand -h)") >>$SWversions 2>&1
(echo "$vcfValidatorCommand;$($vcfValidatorCommand -h)") >> $SWversions 2>&1
(echo "$bcfToolCommand;$($bcfToolCommand -v)") >> $SWversions 2>&1
(echo "$bgzipCommand;$($bgzipCommand -h)") >> $SWversions 2>&1
#(echo "$tabixCommand;$($tabixCommand -h)") >> $SWversions 2>&1
######################################################################################
######################################################################################

######################################################################################
###############Additional resouces needed to run the pipeline#########################
######################################################################################
teePath="/homes/angel/coreutils-8.25-INST/bin/tee";
teePath="tee";
bunzip2Path="/homes/angel/bunzip2";
bunzip2Path="bunzip2";
export PERL5LIB=/usr/local/vcftools-0.1.13/perl/:${PERL5LIB}
#export JAVA_HOME="/nfs/ega/private/ega/production/java/latest"
#export PATH="$JAVA_HOME/bin/:$PATH"
######################################################################################
######################################################################################

######################################################################################
#Output, standard output and standard error files generated during pipeline execution#
######################################################################################
#StatsPlotsFolder="$fileDir/PlotsStats";                mkdir -p $StatsPlotsFolder;
                                                        viewErr="$pipeErrorFolder/view.e";                              uncompressErr="$pipeErrorFolder/uncompress.e";
headerFile="$fileDir/header.txt.openssl.gz";            headerErr="$pipeErrorFolder/header.e";
randomSNPfile="$fileDir/randomSNPs.txt";                randomSNPerr="$pipeErrorFolder/randomSNP.e";
validationFile="$fileDir/validation.txt.gz";            validationErr="$pipeErrorFolder/validation.e";
vcfVersionFile="$fileDir/vcfVersion.txt";               vcfVersionErr="$pipeErrorFolder/vcfVersion.e";
TsTvFile="$fileDir/vcf";                                TsTvErr="$pipeErrorFolder/TsTv.e";
#gb17File="$fileDir/gb17.txt.gz";                       gb17Err="$pipeErrorFolder/gb17.e";
#gb37File="$fileDir/gb37.txt.gz";                       gb37Err="$pipeErrorFolder/gb37.e";
#gb38File="$fileDir/gb38.txt.gz";                       gb38Err="$pipeErrorFolder/gb38.e";
densityFile="$fileDir/density.snpden.gz";               densityErr="$pipeErrorFolder/density.e";
statFile="$fileDir/stats.txt.openssl.gz";               statErr="$pipeErrorFolder/stats.e";
qualityFile="$fileDir/quality.txt.openssl.gz";          qualityErr="$pipeErrorFolder/quality.e";
frqFile="$fileDir/frequency.txt.openssl.gz";            frqErr="$pipeErrorFolder/frequency.e";
#plotErr="$pipeErrorFolder/plots.e";
#mdCheckSumFile="$fileDir/$fileName.md5.txt";
indexFile="$fileDir/out.index";                         indexErr="$pipeErrorFolder/index.e"
#                                                       avroErr="$pipeErrorFolder/avroTransform.e"
#avroConvertMetaFile="temp.vcf"
#avroConvertMetaFileSuffix="file.json.gz"
#avroConvertFile="variants.avro.gz"
#avroThreads="16"
######################################################################################
######################################################################################



EncryptionKey=$(sed '1q;d' $KeyFile);
ACTpattern='vcf.gz.[0-9]*.(cip|gpg)$'
if [[ "$fileName" =~ $ACTpattern || "$fileName" == *"gz.gpg"  ||  "$fileName" == *"gz.cip" ]];then
        read_input="cat - | gunzip -d - 2>$uncompressErr "
elif [[ "$fileName" == *"bz2.gpg" || "$fileName" == *"bz2.cip" ]]; then
        read_input="cat - | $bunzip2Path - 2>$uncompressErr "
else
        if [[ "$fileName" == *"bcf.gpg" || "$fileName" == *"bcf.cip" ]];then
                read_input="cat - | $bcfToolCommand view - 2>$viewErr"
        else
                read_input="cat - "
        fi;
fi;

echo "ACT DEBUG $EncryptionKey"

cwd="$(dirname $0)"
tmpdir="$(mktemp -d)"
{
        {
                cat - | $teePath -p >($bcfToolCommand index - -o $indexFile >$indexErr 2>&1) | eval $read_input | $teePath -p \
                >((cut -f 1,2,3 | awk '$3 ~ /rs/' | shuf -n 1000) > $randomSNPfile 2>$randomSNPerr) \
                >((sed -n '/^#/p' | grep -vi "^#CHROM" | openssl aes-256-cbc -a -salt -k $EncryptionKey | gzip) > $headerFile 2>$headerErr) \
                >((head -n 1 - | awk -F'VCF' '{print $2}') > $vcfVersionFile 2>$vcfVersionErr) \
                >($vcfToolCommand --vcf - --TsTv-summary --out $TsTvFile 2>$TsTvErr) \
                >(($vcfToolCommand --vcf - --SNPdensity 1000 --stdout -c | gzip) >$densityFile 2>$densityErr) \
                >(($vcfToolCommand --vcf - --site-quality --stdout -c | sed '/^#/d' | grep -Pv '\t-1$' | openssl aes-256-cbc -a -salt -k $EncryptionKey | gzip) > $qualityFile 2>$qualityErr) \
                >(($bcfToolCommand stats - | openssl aes-256-cbc -a -salt -k $EncryptionKey | gzip) > $statFile 2>$statErr) \
                >(($vcfValidatorCommand -i stdin -l error -r stdout -o "$fileDir" | grep -v '^Lines read' | cut -d':' -f2 | sort -T /slgpfs/projects/slc00/slc00474/tmp/ -u | gzip) >$validationFile 2>$validationErr) \
                >(($vcfToolCommand --vcf - --freq --stdout -c | sed '/^#/d' | awk '{if (NF!=6) print $0}' | openssl aes-256-cbc -a -salt -k $EncryptionKey | gzip) > $frqFile 2>$frqErr) \
                > /dev/null;
        } 3>&1 >&4 4>&- | cat;
} 4>&1
#               >((/homes/angel/opencga-1.4.0-storage/bin/opencga-storage.sh variant index --transform --input $avroConvertMetaFile --stdin --outdir ${tmpdir} --stdout -Dtransform.threads=$avroThreads | openssl aes-256-cbc -a -salt -k $EncryptionKey | gzip) >$fileDir/${avroConvertFile}.openssl.gz 2>$avroErr) \
#(cat ${tmpdir}/${avroConvertMetaFile}.${avroConvertMetaFileSuffix} | openssl aes-256-cbc -a -salt -k $EncryptionKey) >$fileDir/${avroConvertMetaFileSuffix}.openssl 2>$pipeErrorFolder/transform.e
#               >(($bcfToolCommand query -T $cwd/dic_bcftools_gb17 -f'%CHROM\t%POS\t%REF\t%ALT\n' - | gzip) >$gb17File 2>$gb17Err) \
#               >(($bcfToolCommand query -T $cwd/dic_bcftools_gb37 -f'%CHROM\t%POS\t%REF\t%ALT\n' - | gzip) >$gb37File 2>$gb37Err) \
#               >(($bcfToolCommand query -T $cwd/dic_bcftools_gb38 -f'%CHROM\t%POS\t%REF\t%ALT\n' - | gzip) >$gb38File 2>$gb38Err) \
rm -rf $tmpdir