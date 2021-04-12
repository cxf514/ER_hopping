#!/bin/bash
time1=`date +"%s"`  # record start time
thread_num=1        # number of threads prepare to be uesd(default=1)
extension=".SNPcalling_summarize.txt"  # the extension of file1 generated in the script "01-SNPcalling.sh"
####################### Process the passing arguments ############################
while getopts "t:s:S:" opt       
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
		s)  # pass the path of a SNPcalling_summarize file to be processed
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
		S)  # pass the path of a directory contain the SNPcalling_summarize files to be processed
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
dir2="SNP_counting"  # the directory contain the output file "file2" 
file2="SNPcounting.txt"  # the common part of file2'name
mkdir -p ${dir2}  # make the directory 2
################################ Main process #####################################
for sample in ${file_list}
do
	read -u 4
	{
	name=${sample%${extension}}
        column_seq=`head -1 ${file_dir}/${sample} |awk '{print NF-1}'`  # get the column number of SNP information

# count the number of each types of SNP and indel in one line(one line store the information of reads with same UMI); the column named "INS" recording the number of the max insertion; the column named "[ins_type]" shows all the types of ins as well as their number, such as[+1=19;+3=43]; so does "DEL" and "[del_type]"; meanwhile, find out the SNP or indel with max number in one line and calculate its percentage.
	gawk 'BEGIN{
		OFS=",";
		if('${column_seq}'==6){print "rs_number,chr,read_start_pos,SNP_pos,total_num,A,T,C,G,INS,[ins_type],DEL,[del_type],max_genotype,max_num,max_percent"}
		else{print "UMI,rs_number,chr,read_start_pos,SNP_pos,total_num,A,T,C,G,INS,[ins_type],DEL,[del_type],max_genotype,max_num,max_percent"};
	      };
	      NR!=1{
	     	gsub(/d/,"",$'${column_seq}')
	      }; 
	      (NR!=1) && ($'${column_seq}' !~/[+-]/){
		str1=toupper($'${column_seq}');
		len=length(str1);
		A=gsub(/A/,"&",str1);
		T=gsub(/T/,"&",str1); 
		C=gsub(/C/,"&",str1);
		G=gsub(/G/,"&",str1);
		INS=0;
		DEL=0;
		if (A < T){max=T;max_genotype="T"}else{max=A;max_genotype="A"};
		if (max < C){max=C;max_genotype="C"};
		if (max < G){max=G;max_genotype="G"};
		ins_type=""
		del_type=""
      	      };
	      (NR!=1) && ($'${column_seq}' ~/[+-]/){
		str2=$'${column_seq}';
		gsub(/[^0-9+-]/," ",str2);
		split(str2,indel," ");
		str3=toupper($'${column_seq}');
		for (i in indel){
			a=" "indel[i]" ";
			b="[A-Z]["substr(indel[i],1,1)"]"substr(indel[i],2)"[A-Z]{"substr(indel[i],2)"}";
			sub(b,a,str3)
		};
		str4=str3;	       
		gsub(/[0-9 ]/,"",str4);	      
		len=length(str4);
		A=gsub(/A/,"&",str4);
		T=gsub(/T/,"&",str4); 
		C=gsub(/C/,"&",str4);
		G=gsub(/G/,"&",str4);
		INS=0;
		DEL=0;
		if (A < T){max=T;max_genotype="T"}else{max=A;max_genotype="A"};
		if (max < C){max=C;max_genotype="C"};
		if (max < G){max=G;max_genotype="G"};
		str5=str3;
		gsub(/[A-Z]/,"",str5);
		split(str5,indel2," ");
		for (i in indel2){
			if(indel2[i] ~/+/){ins[indel2[i]]++};
			if(indel2[i] ~/-/){del[indel2[i]]++};
		};
		for (i in ins){
			if(INS < ins[i]){INS=ins[i]};
			if(max < ins[i]){max=ins[i];max_genotype=i};
		};
		for (i in del){
			if(DEL < del[i]){DEL=del[i]};
			if(max < del[i]){max=del[i];max_genotype=i};
		};
		l1=asorti(ins,sort1);
		l2=asorti(del,sort2);
		for (i=1;i<=l1;i++){ins_type=ins_type";"sort1[i]"="ins[sort1[i]]};
		for (i=1;i<=l2;i++){del_type=del_type";"sort2[i]"="del[sort2[i]]};
		ins_type=substr(ins_type,2);
		del_type=substr(del_type,2);
	      };
	      (NR!=1) && (len!=0){
	      	if('${column_seq}'==6){print $1,$2,$3,$4,len,A,T,C,G,INS,"["ins_type"]",DEL,"["del_type"]",max_genotype,max,max/len}
		else{print $1,$2,$3,$4,$5,len,A,T,C,G,INS,"["ins_type"]",DEL,"["del_type"]",max_genotype,max,max/len};
	      };
	      NR!=1{
		ins_type="";
		del_type="";
		delete ins;
		delete del;
	      }' ${file_dir}/${sample} > ${dir2}/${name}.${file2}
	echo "" >&4
	}&
done
wait >/dev/null 2>&1
time2=`date +"%s"`  # record end time
echo "running time = "$[${time2}-${time1}]" seconds"  # output running time
exec 4>&-  # Delete file descriptor

