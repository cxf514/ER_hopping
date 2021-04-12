#!/bin/bash
time1=`date +"%s"`  # record start time
thread_num=1        # number of threads prepare to be uesd(default=1)
####################### Process the passing arguments ############################
while getopts "u:t:s:" opt       
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
		u)  # pass the length of UMI
                        if [ ${OPTARG:0:1} == "-" ]
                        then
                                echo "$0: Must provide an argument for -${opt}"
                                exit
                        fi
                        UMI_length=${OPTARG}
                        ;;
		s)  # pass the directory where the sam file is located
                        if [ ${OPTARG:0:1} == "-" ]
                        then
                                echo "$0: Must provide an argument for -${opt}"
                                exit
                        fi
                        sam_dir=${OPTARG}
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
if [ ! ${UMI_length} ]  # determine whether the necessary argument "-f" and its file exist
then
        echo "$0: Argument [-u value] is necessary, you must provide the length of UMI uesd in your data" 
        exit
fi

dir_RNCU="RNCU"
file="machine_X"
mkdir -p ${dir_RNCU} 
################################ Main process #####################################
file_list=$(ls ${sam_dir}|grep "\.sam$")
for sample in ${file_list}

do
	read -u 4
	{
	name=${sample%".sam"}
	gawk -F "\t" 'BEGIN{
			max=0;
			printf("%s","'${name}'") > "'${name}.${file}_1.txt'";
			printf("%s","'${name}'") > "'${name}.${file}_2.txt'";
		  };
		  $1!~/^@/{
		  	UMI=substr($1,length($1)-"'${UMI_length}'"+1);
		  	C_P_U[$3"_"$4"_"UMI]++;
		  };
		  END{
		  	for(i in C_P_U){
				RNCU[C_P_U[i]]++;
				if(max<C_P_U[i]){max=C_P_U[i]};
				print i"\t"C_P_U[i] >> "'${dir_RNCU}/${name}.RNCU1.txt'";
			};
		  	for(i=1;i<=max;i++){
				if(RNCU[i]){
					printf(",%d",RNCU[i]) > "'${name}.${file}_1.txt'";
					printf(",%d",i) > "'${name}.${file}_2.txt'";
				};
			};
			print "" > "'${name}.${file}_1.txt'";
			print "" > "'${name}.${file}_2.txt'";
			delete C_P_U;
			delete RNCU;
		  }' ${sam_dir}/${sample}
	sort -t "_" -k 1,1 -k 2n,2n  ${dir_RNCU}/${name}.RNCU1.txt > ${dir_RNCU}/${name}.RNCU.txt
	rm ${dir_RNCU}/${name}.RNCU1.txt
	echo "" >&4
	} &
done
wait >/dev/null 2>&1
cat *${file}_1.txt >> ${file}_values.txt
cat *${file}_2.txt >> ${file}_index.txt
rm *${file}_[12].txt

time2=`date +"%s"` # record end time
echo "run time = "$[${time2}-${time1}]" seconds"
exec 4>&-  # Delete file descriptor
