### ONLY run on O2 as an interactive job, 

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
yamlfile=$1
eval $(parse_yaml $yamlfile "config_")

loc=$config_scripts_location
s2=$config_scripts_s2

perl $loc/$s2 $yamlfile
