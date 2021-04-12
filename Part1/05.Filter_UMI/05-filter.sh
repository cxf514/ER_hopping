#!/bin/bash
time1=`date +"%s"`  # record start time
thread_num=1        # number of threads prepare to be uesd(default=1)
####################### Process the passing arguments ############################
while getopts "t:c:s:r:" opt       
do
        case ${opt} in
                t)  # pass the number of thread prepared to used
                        if [ ${OPTARG:0:1} == "-" ]
                        then
                                echo "$0: Must provide an argument for -${opt}"  
                                exit
                        fi
                        thread_num=${OPTARG}
                        ;;
		c)  # pass the path of the predict cutfile
                        if [ ${OPTARG:0:1} == "-" ]
                        then
                                echo "$0: Must provide an argument for -${opt}"
                                exit
                        fi
                        cutfile=${OPTARG}
                        ;;
		s)  # pass the directory where the sam file is located
                        if [ ${OPTARG:0:1} == "-" ]
                        then
                                echo "$0: Must provide an argument for -${opt}"
                                exit
                        fi
                        dir_sam=${OPTARG}
                        ;;
		r)  # pass the directory where the "RNCU" file is located
                        if [ ${OPTARG:0:1} == "-" ]
                        then
                                echo "$0: Must provide an argument for -${opt}"
                                exit
                        fi
                        dir_RNCU=${OPTARG}
                        ;;
                *)
                        echo "passing arguments error"
                        exit
                        ;;
        esac
done
######################### Multi-threaded control ##################################
tmp_fifo_name="/tmp/$$.fifo"
mkfifo ${tmp_fifo_name}
exec 4<>${tmp_fifo_name}
rm -f ${tmp_fifo_name}
for ((i=0;i<thread_num;i++))
do
        echo "" >&4
done
########################## Arguments and Directory ################################
dir_RNCUcut="RNCUcut"
dir_filter="filter_sam"
file="filter.sam"
mkdir -p ${dir_RNCUcut}
mkdir -p ${dir_filter}
################################ Main process #####################################
for line in $(tail -n +2 ${cutfile})
do
	read -u 4
	{
	sample=$(echo ${line}| awk -F "," '{print $1}')
	cutoff=$(echo ${line}| awk -F "," '{print $3}')
	UMI_length=$(awk -F "[_\t]" 'NR==2{print length($3)}' ${dir_RNCU}/${sample}.RNCU.txt)
	awk -F "\t" '$2<='${cutoff}'{print $0}' ${dir_RNCU}/${sample}.RNCU.txt > ${dir_RNCUcut}/${sample}.RNCUcut.txt
	awk -F "\t" 'NR==FNR{C_P_U[$1]=$1;next};{UMI=substr($1,length($1)-"'${UMI_length}'"+1);if(!($3"_"$4"_"UMI in C_P_U)){print $0}}' ${dir_RNCUcut}/${sample}.RNCUcut.txt ${dir_sam}/${sample}.sam > ${dir_filter}/${sample}_${file}
	echo "" >&4
        } &
done
wait >/dev/null 2>&1
time2=`date +"%s"` # record end time
echo "run time = "$[${time2}-${time1}]" seconds"
exec 4>&-  # Delete file descriptor
