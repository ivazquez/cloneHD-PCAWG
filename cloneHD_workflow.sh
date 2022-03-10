#!/bin/bash

# set an initial value for the flag
input_vcf=false
input_cna=false
sample_name=false
output_dir=false
trials=false
restarts=false
seed=false
snv_fprate=false
snv_fpfreq=false
show_help=false
debug=false

if [ $# -eq 0 ];
then
    show_help=true
fi

# read the options
TEMP=`getopt -o v:c:s:o:t:r:x:R:F:hd --long vcf:,cna:,sample:,output:,trials:,restarts:,seed:,snv-fprate:,snv-fpfreq:,help,debug -n 'cloneHD_workflow.sh' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables
while true ; do
    case "$1" in
        -v|--vcf) input_vcf=$2 ; shift 2 ;;
        -c|--cna) input_cna=$2 ; shift 2 ;;
        -s|--sample) sample_name=$2 ; shift 2 ;;
        -o|--output) output_dir=$2 ; shift 2 ;;
        -t|--trials) trials=$2 ; shift 2 ;;
        -r|--restarts) restarts=$2 ; shift 2 ;;
        -x|--seed) seed=$2 ; shift 2 ;;
        -R|--snv-fprate) snv_fprate=$2 ; shift 2 ;;
        -F|--snv-fpfreq) snv_fpfreq=$2 ; shift 2 ;;
        -h|--help) show_help=true ; shift ;;
        -d|--debug) debug=true ; shift ;;
        --) shift ; break ;;
        *) echo $2 "Internal error!" ; exit 1 ;;
    esac
done

if [[    $input_vcf == false
      || $input_cna == false
      || $sample_name == false
      || $output_dir == false
    ]];
then
    show_help=true
fi

if [ $show_help == true ];
then
    echo "usage: $0 [options]"
    echo "options:"
    echo "    -v, --vcf arg       VCF file describing the SNV data"
    echo "    -c, --cna arg       Battenberg-format segmentation file describing the CNA data"
    echo "    -s, --sample arg    sample name"
    echo "    -o, --output arg    write output into this directory named dir/sample.cloneHD.gz"
    echo "    -t, --trials arg    number of independent optimizations [default: 10]"
    echo "    -r, --restarts arg  number of perturbations in local random search mode [default: 10]"
    echo "    -x, --seed arg      seed [default: 123]"
    echo "    -d, --debug         turns on debugging"
    echo "    -h, --help          this text"
    echo "    add --debug for debugging output"
    echo "Run cloneHD and save results in the output directory."
    exit
fi

prefix=$sample_name
snv=$prefix.snv.txt
mean_tcn=$prefix.mean_tcn.txt
avail_cn=$prefix.avail_cn.txt

### SNV parser ###
python3 /opt/cloneHD-PCAWG/cloneHD_snv_parser.py \
	--vcf $input_vcf \
	--vcf-type 'mutect-smchet' \
	--sample 'tumor' \
	--snv $snv

### CNA parser ###
gender=`awk '{if($1==24){sum++}}END{if(sum>5){print "male"}else{print "female"}}' ${snv}`

perl /opt/cloneHD-PCAWG/cloneHD_cna_parser.pl \
	-g $gender \
	-m "mean-tcn" \
	-c $input_cna \
	-o $mean_tcn

perl /opt/cloneHD-PCAWG/cloneHD_cna_parser.pl \
	-g $gender \
	-m "avail-cn" \
	-c $input_cna \
	-o $avail_cn

### cloneHD ###
declare -A clones_to_clusters=( [1]=0 [2]=1 [3]=3 )
n_clones=1

while [ $n_clones -le 3 ]
do
	
	summary[$n_clones]=$prefix.Nc$n_clones.summary.txt
	snv_posterior[$n_clones]=$prefix.Nc$n_clones.snv.posterior.txt
	
	n_clusters=${clones_to_clusters[$n_clones]}
	n_snvs=`wc -l $snv | awk '{print $1}'`	
	rescaled_snv_fprate=`echo "$snv_fprate $n_snvs" | awk '{print $1*3E9/$2}'`
	
	# fixed number of trials and restarts for large number of SNVs
	if [[ "$n_snvs" -ge 25000 &&  "$n_snvs" -le 50000 ]]; then
		trials=8
		restarts=8
	elif [[ "$n_snvs" -ge 50000 &&  "$n_snvs" -le 100000 ]]; then
		trials=6
		restarts=6
	elif [[ "$n_snvs" -ge 100000 &&  "$n_snvs" -le 500000 ]]; then
		trials=4
		restarts=4
	elif [[ "$n_snvs" -ge 500000 ]]; then
		trials=2
		restarts=2
	else
		:
	fi
	
	/opt/cloneHD/build/cloneHD \
		--pre $prefix.Nc$n_clones \
		--snv $snv \
		--seed $seed \
		--force $n_clones \
		--trials $trials \
		--restarts $restarts \
		--max-tcn 8 \
		--mean-tcn $mean_tcn \
		--avail-cn $avail_cn \
		--snv-rnd 1E-4 \
		--snv-fpfreq $snv_fpfreq \
		--snv-fprate $rescaled_snv_fprate \
		--learn-cluster-w $n_clusters \
		--snv-pen-high 3E-1 \
		--print-all 0

	n_clones=`expr $n_clones + 1`
	
done

## Model selection ###
perl /opt/cloneHD-PCAWG/cloneHD_model_selection.pl \
    -i ${summary[1]} -j ${snv_posterior[1]} \
    -k ${summary[2]} -l ${snv_posterior[2]} \
    -m ${summary[3]} -n ${snv_posterior[3]} \
    -a 10.0 -s 50.0 \
    -o $prefix