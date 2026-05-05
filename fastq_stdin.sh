
#!/bin/bash

currentDir=$1; #KK
fileDir=$2; #output
fileName=$(echo $3 | sed 's/\.[0-9]\+\././g');
filePath=$4; #KK
pipeErrorFolder=$5;
KeyFile=$6;
typeFile=$7


SWversions="$fileDir/SW_versions"
> $SWversions
teePath="/homes/angel/coreutils-8.25-INST/bin/tee";
teePath="tee";
bunzip2Path="/homes/angel/bunzip2";
bunzip2Path="bunzip2";
FastqcCommand="$currentDir/PipelineScripts/FastQC/fastqc";
FastqcCommand="/usr/local/FastQC/fastqc";
FastqscreenCommand="/usr/local/FastQ-Screen/fastq_screen";
#FastqcCommand="/homes/angel/FastQC_v0.11.8/fastqc";
echo "$FastqcCommand;$($FastqcCommand -v 2>&1)" >> $SWversions
errorFile="$pipeErrorFolder/fastqc.e";          outFile="$pipeErrorFolder/fastqc.o";
errorFile2="$pipeErrorFolder/unzip.e";
errorFile3="$pipeErrorFolder/fastqscreen.e";    outFile3="$pipeErrorFolder/fastqscreen.o";
uncompressErrorFile="$pipeErrorFolder/uncompress.e";

fastq_mem="$((SLURM_CPUS_PER_TASK*3900))"
fastq_mem_under_9000="$((fastq_mem<9000 ? fastq_mem : 9000))"

EncryptionKey=$(sed '1q;d' $KeyFile);
if [[ "$fileName" == *"gz.gpg" || "$fileName" == *"gz.cip" ]];then
        sec="gunzip -d - 2>$uncompressErrorFile"
elif [[ "$fileName" == *"bz2.gpg" || "$fileName" == *"bz2.cip" ]]; then
        sec="$bunzip2Path - 2>$uncompressErrorFile"
else
        sec="cat"
fi;

{
        {
                cat - | eval $sec | $teePath -p >(awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\n");} else { printf("TAB-MARK");} }' | awk -v k=100000 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | awk -F"TAB-MARK" '{print $1"\n"$2"\n"$3"\n"$4}' | gzip > $fileDir/input-100K.fastq.gz) >(awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\n");} else { printf("TAB-MARK");} }' | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .01) print $0}' | awk -F"TAB-MARK" '{print $1"\n"$2"\n"$3"\n"$4}' | gzip > $fileDir/input-1percent.fastq.gz) >(wc -l >$fileDir/howmany.lines) >($FastqcCommand /dev/stdin -o "$fileDir" -f fastq --memory $fastq_mem_under_9000 >$outFile 2>$errorFile) >/dev/null
        } 3>&1 >&4 4>&- | cat;
} 4>&1

fastqcZip="$fileDir/stdin_fastqc.zip";
#unzip $fastqcZip -d $fileDir 1>/dev/null 2>$errorFile2;


#$FastqscreenCommand --quiet --conf ~/qc/fastq_screen.conf $fileDir/input.fastq.gz --outdir $fileDir >$outFile3 2>$errorFile3

#rm $fileDir/input.fastq.gz $fileDir/input_screen.html