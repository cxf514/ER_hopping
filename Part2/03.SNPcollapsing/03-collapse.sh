#!/bin/bash
time1=`date +"%s"`  # record start time
thread_num=1  # number of threads prepare to be uesd(default=1)
collapsed_percent=0.9  # consistency of one UMI;below this value the UMI will be discarded(default=0.9) 
extension=".SNPcounting.txt"  # the extension of file2 generated in the script "02-SNPcounting.sh"
####################### Process the passing arguments ############################
while getopts "t:p:c:s:S:" opt       
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
		p)  # pass the collapsed percent(0-1,default=0.9) 
			if [ ${OPTARG:0:1} == "-" ]
			then
				echo "$0: Must provide an argument for -${opt}"
				exit
			fi
			collapsed_percent=${OPTARG}
			;;
		c)  # decide whether to cutoff and pass the path of "cutoff" file(default=False)
			if [ ${OPTARG:0:1} == "-" ]
			then
				echo "$0: Must provide an argument for -${opt}"
				exit
			fi
			cutoff=${OPTARG}
			cutoff_need=1			
			;;
		 s) #  pass the path of a SNPcounting file to be processed
                        if [ ${OPTARG:0:1} == "-" ]
                        then
                                echo "$0: Must provide an argument for -${opt}"
                                exit
                        fi
                        if [ ${S_done} ]
                        then
                                echo "$0: argument \"-s\" and \"-S\" can not be used together"
                                exit
                        fi
			file_dir=$(dirname ${OPTARG})
			file_list=$(basename ${OPTARG})
			s_done=1
                        ;;
                S)  # pass the path of a directory contain the SNPcounting files to be processed
                        if [ ${OPTARG:0:1} == "-" ]
                        then
                                echo "$0: Must provide an argument for -${opt}"
                                exit
                        fi
                        if [ ${s_done} ]
                        then
                                echo "$0: argument \"-s\" and \"-S\" can not be used together"
                                exit
                        fi
			file_dir=${OPTARG}
			file_list=$(ls ${file_dir}|grep "${extension}$")
                        S_done=1
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
dir3="collapsed"
file3="collapsed.txt"
if [ ${cutoff_need} ] && [ ! -f ${cutoff} ]
then
        echo -e "cutoff file:\"${cutoff}\" not exist\nuse [-c value] to perform cutoff, but the cutoff file must be provided"
        exit
fi

mkdir -p ${dir3} # make the directory 3
################################ Main process #####################################
for sample in ${file_list}
do
	read -u 4
	{
	name=${sample%${extension}}
        if [ ${cutoff_need} ]
        then
		cut=`awk -F, '$1 ~/^'${name}'$/{found=1;print $2};END{if(found!=1){print 0}}' ${cutoff}`
        else
                cut=0
        fi
# read the "cufoff" file to get the cutoffs of each samples(cut=0 if not match); then discard the line(one line store the imformation of reads with same UMI) which "max_percent(column 17)" lower than "collapsed_percent" or "length(column 6)" lower than cutoff; use the "max_genotype" as representative of one line and count each SNP number(ATCGindel) in one locus
	for each in ${cut}
	do	
		gawk -F, 'BEGIN{
				OFS=",";
				print "rs_number,chr,read_start_pos,SNP_pos,total_UMI_num,A,T,C,G,INS,[ins_type],DEL,[del_type]";
	  		  };
			  (NR!=1) && ($17>='${collapsed_percent}') && ($6>'${each}'){
				if(len[$2","$3","$4","$5]==""){n[j++]=$2","$3","$4","$5};
				len[$2","$3","$4","$5]++;
				if(a[$2","$3","$4","$5]==""){a[$2","$3","$4","$5]=0};
				if($15=="A"){a[$2","$3","$4","$5]++};
				if(t[$2","$3","$4","$5]==""){t[$2","$3","$4","$5]=0};
				if($15=="T"){t[$2","$3","$4","$5]++};
				if(c[$2","$3","$4","$5]==""){c[$2","$3","$4","$5]=0};
				if($15=="C"){c[$2","$3","$4","$5]++};
				if(g[$2","$3","$4","$5]==""){g[$2","$3","$4","$5]=0};
	    			if($15=="G"){g[$2","$3","$4","$5]++};
				if(ins[$2","$3","$4","$5]==""){ins[$2","$3","$4","$5]=0};
				if($15 ~/+/){ins_arr[$2","$3","$4","$5"|"$15]++};
				if(del[$2","$3","$4","$5]==""){del[$2","$3","$4","$5]=0};
	 			if($15 ~/-/){del_arr[$2","$3","$4","$5"|"$15]++};
			   };
			   END{
			     	l1=asorti(ins_arr,sort1);
				l2=asorti(del_arr,sort2);
				for(i=1;i<=l1;i++){
					split(sort1[i],arr1,"|");
					ins_type[arr1[1]]=ins_type[arr1[1]]";"arr1[2]"="ins_arr[sort1[i]];
					if(ins[arr1[1]]<ins_arr[sort1[i]]){ins[arr1[1]]=ins_arr[sort1[i]]};
				};
				for(i=1;i<=l2;i++){
					split(sort2[i],arr2,"|");
					del_type[arr2[1]]=del_type[arr2[1]]";"arr2[2]"="del_arr[sort2[i]];
					if(del[arr2[1]]<del_arr[sort2[i]]){del[arr2[1]]=del_arr[sort2[i]]};
				};
				for(i=0;i<j;i++){
					print n[i],len[n[i]],a[n[i]],t[n[i]],c[n[i]],g[n[i]],ins[n[i]],"["substr(ins_type[n[i]],2)"]",del[n[i]],"["substr(del_type[n[i]],2)"]";
				}	
			   }'  ${file_dir}/${sample} > ${dir3}/${name}.${each}.${file3}
	done 
	echo "" >&4
	} &
done
wait >/dev/null 2>&1
time2=`date +"%s"` # record end time
echo "running time = "$[${time2}-${time1}]" seconds"  # output running time
exec 4>&-  # Delete file descriptor


