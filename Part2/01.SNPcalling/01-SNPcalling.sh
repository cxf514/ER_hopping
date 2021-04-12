#!/bin/bash
time1=`date +"%s"`  # record start time
thread_num=1  # number of threads prepare to be uesd(default=1)
Qscore=30  # Phred Quality Score (default=30)
UMI_length=0  # The length of UMI used in the sam file, 0 represent not used (default=0)
extension=".sam" # extension of sam file
####################### Process the passing arguments ############################
while getopts "t:Q:f:u:s:S:" opt       
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
		Q)  # pass the minimum Qscore threshold
			if [ ${OPTARG:0:1} == "-" ]                             
			then                                                                    
				echo "$0: Must provide an argument for -${opt}"
				exit    
			fi
			Qscore=${OPTARG}
			;;
		f)  # pass the path of target loci file
			if [ ${OPTARG:0:1} == "-" ]                             
			then                                                                    
				echo "$0: Must provide an argument for -${opt}"
				exit    
			fi
			loci_file=${OPTARG}
			;;		
		u)  # pass the length of UMI
			if [ ${OPTARG:0:1} == "-" ]                             
			then
				echo "$0: Must provide an argument for -${opt}"
				exit
			fi
			UMI_length=${OPTARG}
			;;
		s)  # pass the path of a sam file to be processed
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
			sam_dir=$(dirname ${OPTARG})
			sam_list=$(basename ${OPTARG})
			s_done=1
			;;
		S)  # pass the path of a directory contain the sam files to be processed
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
			sam_dir=${OPTARG}
			sam_list=$(ls ${sam_dir}|grep "${extension}$")
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
dir0="SNP_calling"  # the directory contain the output file "file0", which record the extract SNP of each read 
file0="SNPcalling.txt"  # the common part of file0'name
dir1="SNP_calling_summarize"  # the directory contain the summarize output of "file0"
file1="SNPcalling_summarize.txt"  # the common part of file1'name
Phred="z!""\"""#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"  # ascii 33~74, the "z" is used to help the system recognize "!" as a string
Qscore=$((${Qscore}+1))
Phred=${Phred:0:${Qscore}}  # select the corresponding ascii according to the Qscore
if [ ! ${loci_file} ]  # determine whether the necessary argument "-f" and its file exist
then
	echo "$0: Argument [-f value] is necessary" 
	exit
else
	if [ ! -f ${loci_file} ]
	then	
		echo -e "$0:\nfile \"${loci_file}\" not exist\nyou need to provide a file like this:\nchr1 1050069 rs2799067\nchr6 87472986 rs9450704\nchr12 4662516 rs2267549\n......\n"
		exit
	fi
fi
mkdir -p ${dir0}
mkdir -p ${dir1} # make the directory 0 and 1
################################ Main process #####################################
for sample in ${sam_list}
do
	read -u 4
	{
	name=${sample%${extension}}  # get the file name without extension

# Since the 'gawk' cannot comment inside the code, so each part are commented at the beginning.
# "BEGIN": output the column name
# "NR==FNR": read the "${file_name}" and store its second and third columns in the array "pos" and "rs" respectively. The array "chr" is used to reduce the number of subsequent searches.
# "(NR>FNR) && ($1!~/^@/) && ($5>=30) && ($5<=60)": read the sam file and process the row without "@" at the beginning and 30 <= MAPing Quality(column 5) <= 60
# SNP calling for each reads, "[M=XHP]" means the sequence has no indel at this base; "[S]" means there is an insertion at the start of sequence; "[I]" means insertion; "[D]" means deletion; "[N]" means intron, usuallly appear in RNAseq. Restore the sequence and its sequencing quality(column 10 and 11) according to the reference so that we can extract the base by chromosome location; for insertion[I], its sequence will be store in the nearest base before it, and can get it by calling this base; so is deletion[D], except its sequence will replace with "D" because the sam file does not provide its sequence. If the chromosome location of read exist in the target loci file, the base information will be extract and output; otherwise, its chromosome and start location (column 3 and 4) will be output.
	gawk 'BEGIN{
		if('${UMI_length}'==0){print "QNAME rs_number chr read_start_pos SNP_pos SNP qual" > "'${dir0}/${name}.${file0}'"}
		else{print "QNAME UMI rs_number chr read_start_pos SNP_pos SNP qual" > "'${dir0}/${name}.${file0}'"};
	     };

	     NR==FNR{
		chr[$1]++;
		pos[$1","chr[$1]]=$2;
		rs[$1","chr[$1]]=$3;
	     };
			
	     (NR>FNR) && ($1!~/^@/) && ($5>=30) && ($5<=60){
		str1=$6;
		gsub(/[MIDNSHPX=]/," ",str1);
		length_len=split(str1,len," ");
		str1=$6;
		gsub(/[0-9]/," ",str1);
		split(str1,type," ");
		n1=1;
		n2=1;
		for(i1=1;i1<=length_len;i1++){
			if(type[i1] ~ /[M=XHP]/){
				seq1=seq1""substr($10,n1,len[i1]);
				seq2=seq2""substr($11,n1,len[i1]);
				n1+=len[i1];
				n2+=len[i1];
			}
			else if(type[i1] ~ /[S]/){
				n1+=len[i1];
			}
			else if(type[i1] ~ /[I]/){
				indel[$4+n2-2]="+"len[i1]""substr($10,n1,len[i1]);
				n1+=len[i1];
			}     
			else if(type[i1] ~ /[DN]/){
				for(i2=0;i2<len[i1];i2++){
					D1=D1"d";
					D2=D2"D";
				};
				seq1=seq1""D1;
				seq2=seq2""D1;
				indel[$4+n2-2]="-"len[i1]""D2;
				n2+=len[i1];
			};
		};
		end=$4+length(seq1);
		for(i3=1;i3<=chr[$3];i3++){
			if(pos[$3","i3]>=$4 && pos[$3","i3]<end){
				SNP=substr(seq1,pos[$3","i3]-$4+1,1)""indel[pos[$3","i3]];
				quality=substr(seq2,pos[$3","i3]-$4+1,1);
				rs_number=rs[$3","i3];
				SNP_pos=pos[$3","i3];
				break;
			};
			if(i3==chr[$3]){
				notfound[$3"-"$4]++;
				SNP="";
				quality="";
				rs_number="";
				SNP_pos="";
			};
		};
		if(chr[$3]==""){
			notfound[$3"-"$4]++;
			SNP="";
			quality="";
			rs_number="";
			SNP_pos="";
		}
		if(SNP!="" && quality!~/['${Phred}']/){
			if('${UMI_length}'==0){print $1,rs_number,$3,$4,SNP_pos,SNP,quality > "'${dir0}/${name}.${file0}'"}
			else{print $1,substr($1,length($1)-"'${UMI_length}'"+1),rs_number,$3,$4,SNP_pos,SNP,quality > "'${dir0}/${name}.${file0}'"};
		};
		seq1="";
		seq2="";
		D1="";
		D2="";
		delete indel;
	     };
	     END{
		for(i4 in notfound){print "'${name}:' "i4,notfound[i4],"reads not match" >> "'${dir0}/reads_not_match1.txt'"};
	     }' ${loci_file} ${sam_dir}/${sample}   
	
	

# gather together the SNP information of reads with same UMI into one line
	gawk 'BEGIN{
	      if('${UMI_length}'==0){print "rs_number chr read_start_pos SNP_pos number SNP qual" > "'${dir1}/${name}.${file1}'"}
		else{print "UMI rs_number chr read_start_pos SNP_pos number SNP qual" > "'${dir1}/${name}.${file1}'"}
	      };
	      NR!=1{
		if('${UMI_length}'==0){
	     		if(b[$2" "$3" "$4" "$5]==""){a[j++]=$2" "$3" "$4" "$5};
			b[$2" "$3" "$4" "$5]++;
			c[$2" "$3" "$4" "$5]=c[$2" "$3" "$4" "$5]""$6;
			d[$2" "$3" "$4" "$5]=d[$2" "$3" "$4" "$5]""$7;
		}else{
			if(b[$2" "$3" "$4" "$5" "$6]==""){a[j++]=$2" "$3" "$4" "$5" "$6};
			b[$2" "$3" "$4" "$5" "$6]++;
			c[$2" "$3" "$4" "$5" "$6]=c[$2" "$3" "$4" "$5" "$6]""$7;
			d[$2" "$3" "$4" "$5" "$6]=d[$2" "$3" "$4" "$5" "$6]""$8;
		};
       	      };
	      END{
	     	for(i=0;i<j;i++){print a[i],b[a[i]],c[a[i]],d[a[i]] > "'${dir1}/${name}.${file1}'"};
		print "'${name}' done";
	      }' ${dir0}/${name}.${file0}  
	echo "" >&4
	} &
done
wait >/dev/null 2>&1  # Do not output the prompt content
echo "sample_name:chr-pos,read_number" > ${dir0}/reads_not_match.txt
if [ -f ${dir0}/reads_not_match1.txt ]
then
	sort ${dir0}/reads_not_match1.txt >> ${dir0}/reads_not_match.txt
	rm ${dir0}/reads_not_match1.txt
fi
time2=`date +"%s"`  # record end time
echo "running time = $[${time2}-${time1}] seconds"  # output running time
exec 4>&-  # Delete file descriptor


