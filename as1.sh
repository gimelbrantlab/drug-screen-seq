#!/bin/bash
## This performs ASYMMETRIC alignment:
# first, aligns only long reads from R2 using bowtie,
# then adds short reads from R1 (with appropriate FLAG setting) to SAM file
#
#~ bash s1.sh config.yaml
###################################################

module load gcc/6.2.0
module load samtools/1.3.1
module load bowtie2/2.2.9
module load perl/5.24.0 

parse_yaml() {
## from https://gist.github.com/pkuczynski/8665367
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\)\($w\)$s:$s\"\(.*\)\"$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}


# read yaml file
eval $(parse_yaml $1 "config_")

# check version
ver_expected="6.0"
if [[ $ver_expected = $config_version ]]
then 
	echo ""
else
	echo "ERROR $0: wrong version of YAML file:"
	echo ">$ver_expected< expected, >$config_version< actually"
	echo
	exit 2
fi	

# access yaml content
where=$config_alignment_fastq_dir
f1=$config_alignment_R1
f2=$config_alignment_R2

ref=$config_alignment_reference
#### the following is a kludge for o2 index location
ref=/n/groups/shared_databases/bowtie2_indexes/$ref

sam=$config_alignment_SAM_name_base
workdir=$config_alignment_SAM_location

loc=$config_scripts_location
asym=$config_scripts_asymagic
##### now just do work

if [ -d $workdir ] 
then
    echo "Using existing working directory [$workdir]"
else
	mkdir $workdir
fi
cd $workdir
h=`pwd`
here="$h/"

echo "Aligning..."

bowtie2 --threads 8 -q --fast-local --reorder -x $ref -U ${where}$f2 -S ${sam}_r2only.sam

echo "Asymmetric magic..."
perl $loc/$asym ${sam}_r2only.sam ${where}$f1 > ${sam}.sam

echo -n "Samtools transforms... sam->bam... "

samtools view -ubSh $here$sam.sam > $here${sam}.bam
echo -n "sort... "
samtools sort -l 0 -o $here${sam}.sorted.bam $here${sam}.bam 
echo -n "index... "
samtools index $here${sam}.sorted.bam $here${sam}.sorted.bai
echo  "done "
rm $here${sam}.bam

cd ..
