
########### Here is the one line that needs to be modified
Illumeta=/cluster/project8/vyp/sofia/IlluMeta_new


########## everything else below should be automated
Software=${Illumeta}/exec
fastqCollapse=${Software}/fastqCollapse/fastqCollapse
fastqCollapse2=${Software}/fastqCollapse2/fastqCollapse
concatReads=${Software}/SHERAc/concatReads
blastn=${Software}/ncbi-blast-2.2.29+/bin/blastn
blastx=${Software}/ncbi-blast-2.2.29+/bin/blastx
novoalign=${Software}/novoalign
velveth=${Software}/velvet_1.2.10/velveth
velvetg=${Software}/velvet_1.2.10/velvetg
perlExtractContigs=${Illumeta}/scripts/extract_largeContigs.pl
perlExtractReads=${Software}/velvet_1.2.10/contrib/extractContigReads/extractContigReads.pl 

  
############ default values
local=interactive
contigLengthCutoff=150
script=${output}/cluster/submission/default.sh

step1=FALSE
step2=FALSE
step3=FALSE
step4=FALSE
step5=FALSE
step6=FALSE
step7=FALSE
blastnDone=FALSE


until [ -z "$1" ]; do
	# use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--inputFiles)
	    shift
	    i=0
	    for fileloc in $@; do 
		inputFiles[ $i ]=$fileloc
		((i=i+1))
	    done;;
	--local )
	    shift
	    local=$1;;
	--outDir )
	    shift
	    output=$1;;
	--type )
	     shift
	     type=$1;;
	--delimiter )
	    shift
	    delimiter=$1;;
	--script)
	    shift
	    script=$1;;
	 --reference)
            shift
            reference=$1;;
	 --db )
	   shift
	   db=$1;;
	--dbViral)
	    shift
	    dbViral=$1;;
	--dbnr)
	    shift
	    dbnr=$1;;
         --sample )
           shift
           sample=$1;;
	--step1 )
	    shift
	    step1=$1;;
	--step2 )
	    shift
	    step2=$1;;
	--step3 )   
           shift
           step3=$1;;
	--stepMerge )   
           shift
           stepMerge=$1;;
	--step4 )
	    shift
	    step4=$1;;
	--step5 )
	    shift
	    step5=$1;;
	--step6 )
	    shift
	    step6=$1;;
        --step7 )
            shift
            step7=$1;;
	--step8 )
	    shift
	    step8=$1;;
	--step9 )
	    shift
	    step9=$1;;
        --step10 )
            shift
            step10=$1;;
	--step11 )
	    shift
	    step11=$1;;
	--merged )
	    shift
	    merged=$1;;
	--blastnDone )
	    shift
	    blastnDone=$1;;	
	--paired )
	    shift
	    paired=$1;;	

	-* )
	    echo "Unrecognized option: $1"
	    exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done 


################ creating all the output folders
echo -e "Output folder: $output\n"



clusterDir=${output}/cluster_noexpCov
cluster_out=${clusterDir}/out
cluster_error=${clusterDir}/error
cluster_submission=${clusterDir}/submission

output_qc=${output}/QC
output_collapsed=${output}/collapsed
output_merged=${output}/merged
output_novoalign=${output}/novoalign

output_blastn=${output}/blastn
output_blastnrRNA=${output}/blastn_rRNA

output_velvet=${output}/velvet_noexpCov
output_velvet_input=${output_velvet}/input
contigs_subsetted=${output_velvet}/contigs_subsetted
readsContigs=${output_velvet}/readsContribContigs




myFolders="$output $clusterDir $cluster_out $cluster_error $cluster_submission $output_velvet $output_velvet_input $output_merged $output_qc $output_collapsed  $output_blastn $output_blastnrRNA   $output_novoalign  $contigs_subsetted $readsContigs"

for folder in $myFolders; do
    if [ ! -e $folder ]; then 
	echo "Creating $folder"
	mkdir $folder
    fi
done




################################################################### Now writing the script

echo "Output script:  $script"

echo "
#!/bin/bash
#$ -S /bin/bash

date ##to measure the duration


" > $script



######################################################################################### STEP1 (collapsing)   #####################################################                                                                     
if [[ "$step1" == "TRUE" ]]; then

    if [[ "$paired" == "TRUE" ]]; then    ########################## step1 for PAIRED END data                                                                                                               

	
        nfiles=${inputFiles[0]}
	
        if [[ "$nfiles" != "2" ]]; then
            echo "You specified $nfiles input files".
            echo "Error: currently the input data MUST be paired end."
            exit;
        fi
	
	
        seq1=${inputFiles[ 1 ]}
        seq2=${inputFiles[ 2 ]}
	
	
###############  check that raw data files & reference exist                                                                                                                                   
        for file in $seq1 $seq2 $reference; do
            ls -lh $file
            if [ ! -e "$file" ]; then
                echo "Error, file $file does not exist"
                exit
            fi
        done



	
#################################################################  Collapsing                                                                                                                  
        echo "1) Collapse paired reads"
	
        summaryfile=${sample}_summary.txt
	
###change the time depending on the size of dataset

       echo "
#!/bin/bash -l
#$ -S /bin/bash
#$ -o $cluster_out 
#$ -e $cluster_error
#$ -l h_rt=23:00:00
#$ -l tmem=14G,h_vmem=14G
#$ -N step1
#$ -wd  ${output}
#$ -V

${fastqCollapse2}  -i $seq1 $seq2 -o ${output_collapsed}/${sample} -sigLength 40 -summary ${output_collapsed}/${summaryfile}                                                                                  

                                                                                                                                                                                               
" >> $script

qsub $script



   else   ##################check we have single end files 

echo "check single files"
    fi ####paired

fi  ################## end of step1








############################################### STEP2 #########################################################################################

if [[ "$step2" == "TRUE" ]]; then

    if [[ "$paired" == "TRUE" ]]; then    ########################## step2 for PAIRED END data                                                                                                               


	echo "2) PrinSeq reads"

	echo "
#!/bin/bash -l
#$ -S /bin/bash
#$ -o $cluster_out 
#$ -e $cluster_error
#$ -l h_rt=24:00:00
#$ -l tmem=14G,h_vmem=14G
#$ -N step2
#$ -wd  ${output}
#$ -V

                                                                                                      
export PERL5LIB=\${PERL5LIB}:${Software}/bioperl-live            
                                                                                                                                                                                                              
perl ${Software}/prinseq-lite-0.20.3/prinseq-lite.pl -fastq ${output_collapsed}/${sample}_1.fq -fastq2 ${output_collapsed}/${sample}_2.fq  -min_qual_mean 15 -out_good ${output_qc}/Sample${sample}  -trim_qual_right 10 -trim_qual_left 10  -lc_method dust -lc_threshold 7 -no_qual_header -qual_noscale -graph_data ${output_qc}/Sample${sample}_quality.plot  -out_bad null                                                                                                                                             
cat ${output_qc}/Sample${sample}_1_singletons.fastq  ${output_qc}/Sample${sample}_2_singletons.fastq  > ${output_qc}/Sample${sample}_singletons.fq 
                                                                                                                                                                                                             
" >> $script


	qsub $script

    else   ##############  here I will put QC step for single reads



        nfiles=${inputFiles[0]}
	
        if [[ "$nfiles" != "1" ]]; then
            echo "You specified $nfiles input files".
            echo "Error: currently the input data MUST be single end."
            exit;
        fi
	
	
        seq1=${inputFiles[ 1 ]}
	
	
###############  check that raw data files & reference exist                                                                                                                                   
        for file in $seq1 $reference; do
            ls -lh $file
            if [ ! -e "$file" ]; then
                echo "Error, file $file does not exist"
                exit
            fi
        done



	echo "
#!/bin/bash -l
#$ -S /bin/bash
#$ -o $cluster_out 
#$ -e $cluster_error
#$ -l h_rt=12:00:00
#$ -l tmem=14G,h_vmem=14G
#$ -N step2
#$ -wd  ${output}
#$ -V

                                                                                                      
export PERL5LIB=\${PERL5LIB}:${Software}/bioperl-live            
                                                                                                                                                                                                              
zcat ${seq1} | perl ${Software}/prinseq-lite-0.20.3/prinseq-lite.pl -fastq stdin -min_qual_mean 15 -out_good ${output_qc}/Sample${sample}  -trim_qual_right 10 -trim_qual_left 10  -lc_method dust -lc_threshold 7 -no_qual_header -qual_noscale -graph_data ${output_qc}/Sample${sample}_quality.plot  -out_bad null                                                                                                                                                                                                                                                         
" >> $script


	qsub $script



    fi   ###end of paired

fi   ############ end of step2
	
	







######################################################################## STEP_merge ##############################################################################

if [[ "$stepMerge" == "TRUE" ]]; then

    if [[ "$paired" == "TRUE" ]]; then    ########################## step_merge  for PAIRED END data                                                                                                               


##########################################################################Merging                                                                                                                       
       echo "2) Merge Reads"
	
echo "

#!/bin/bash -l
#$ -S /bin/bash
#$ -o $cluster_out 
#$ -e $cluster_error
#$ -l h_rt=03:00:00
#$ -l tmem=20G,h_vmem=20G
#$ -N stepMerge
#$ -wd  ${output}
#$ -V

${concatReads} --adaptersFile ${Illumeta}/support/adapters.fa ${output_qc}/Sample${sample}_1.fastq  ${output_qc}/Sample${sample}_2.fastq  ${output_merged}/${sample}                               

" >> $script
	
		
################################################################ Check whether merging necessary                                                                                               	
        echo -e "3) Is merging necessary? Print merged percentage."	
        echo "                                                                                                                                                                                 
                                                                                                                                                                                               
awk -F":" '{if ((\$1~/Sequences Processed/) || (\$1~/Sequences Overlapping/)) print \$2}' ${output_merged}/$sample.summary  | awk '{ a[\$0] } END {for (i in a){for (j in a){if (i+0 < j+0)  print (i/j)*100} }}' > ${output}/percentage_merged.txt                                                                                                                                              
                                                                                                                                                                                               
" >> $script
	
qsub $script	

    fi
fi   ####### End of stepMerge
	











########################################################################### STEP3 - Novoalign ###################################
if [[ "$step3" == "TRUE" ]]; then

    if [[ "$paired" == "TRUE" ]]; then    ########################## step 3 for PAIRED END data                                                                                                               

        echo -e "4a) Merged data more than 50% -> Align merged (single end) and remaining (paired end) against host."
        echo -e "4b) Merged data less than 50% -> Align original collapsed files (only paired end data) against host."

sed -i 's/ /_/g' ${output_qc}/Sample${sample}_1.fastq
sed -i 's/ /_/g' ${output_qc}/Sample${sample}_2.fastq       
sed -i 's/ /_/g' ${output_qc}/Sample${sample}_singletons.fq       

	
 	echo "

#!/bin/bash -l
#$ -S /bin/bash
#$ -o $cluster_out 
#$ -e $cluster_error
#$ -l h_rt=06:00:00
#$ -pe smp 12
#$ -l tmem=1.2G,h_vmem=1.2G
#$ -N step3
#$ -wd  ${output}
#$ -V

	
# ################################################################# Novoalign                                                                                                                    
	
percentage=\`awk -F\".\" '{print \$1}' ${output}/percentage_merged.txt\`                                                                                                                       
                                                                                                                                                                                               
if [[ \$percentage -ge 50 ]]; then                                                                                                                                                             
############## Hardclipping but no quality calibration                                                                                                                                         
${novoalign} -c 12 -t180  -H -a -e 500  -F STDFQ -f ${output_merged}/$sample.fq -d $reference > ${output_novoalign}/${sample}_overlapping.novo                                                              
${novoalign} -c 12 -t250 -H -a -e 500 -F STDFQ -f ${output_merged}/${sample}_single_1.fq ${output_merged}/${sample}_single_2.fq -d $reference > ${output_novoalign}/${sample}_nonOverlapping.novo         
${novoalign} -c 12 -t180  -H -a -e 500 -F STDFQ -f ${output_qc}/Sample${sample}_singletons.fq -d $reference > ${output_novoalign}/${sample}_singletons.novo     
                                                                                                                                                                                               
else                                                                                                                                                                                           
############# Hardclipping and quality calibration                                                                                                                                             
${novoalign} -c 12 -t250 -H -a -k -e 500 -F STDFQ -f ${output_qc}/Sample${sample}_1.fastq ${output_qc}/Sample${sample}_2.fastq -d $reference > ${output_novoalign}/${sample}.novo                              
${novoalign} -c 12 -t180  -H -a -e 500  -F STDFQ -f ${output_qc}/Sample${sample}_singletons.fq -d $reference > ${output_novoalign}/${sample}_singletons.novo                                                          
   
                                                                                                                                                                                               
fi                                                                                                                                                                                             

" >> $script                                                                                                                                                                                               

        echo -e "5) Select pairs of reads that are both NM ->make txt and fasta files. When merged, select all NM ones."


echo "
	
novofiles="${output_novoalign}/*.novo"                                                                                                                                                         
echo "Novoalign files are '$novofiles'"                                                                                                                                                        
                                                                                                                                                                                               
                                                                                                                                                                                               
for i in \$novofiles                                                                                                                                                                           
  do                                                                                                                                                                                             
  echo "Print '$i'"                                                                                                                                                                          
  filename=\`basename \$i .novo\`                                                                                                                                                            
  echo "Filename is '$filename'"                                                                                                                                                             
                                                                                                                                                                                               
                                                                                                                                                                                               
  if [[ \$filename == *_overlapping* ]]; then
      awk -F\"\\t\" 'BEGIN {OFS=\"\\t\"} {if ((\$5~/NM/) && (length(\$3)>50)) print \$1,\$3}'  \$i |  uniq > ${output_novoalign}/NM_\$filename.txt
      wc -l ${output_novoalign}/NM_\$filename.txt > ${output_novoalign}/NM.number
      awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' ${output_novoalign}/NM_\$filename.txt > ${output_novoalign}/NM_\$filename.fasta                                                                        
                                                                                                                                                                                               
  elif    [[ \$filename == *_singletons* ]]; then
      awk -F\"\\t\" 'BEGIN {OFS=\"\\t\"} {if ((\$5~/NM/) && (length(\$3)>50)) print \$1,\$3}'  \$i |  uniq > ${output_novoalign}/NM_\$filename.txt
      wc -l ${output_novoalign}/NM_\$filename.txt > ${output_novoalign}/NM_singletons.number
      awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' ${output_novoalign}/NM_\$filename.txt > ${output_novoalign}/NM_\$filename.fasta                                                                             
          
  else                                                                                                                                                                              
      awk 'BEGIN {FS=OFS=\"\\t\"} {if ((\$5~/NM/) && (length(\$3)>50))  print \$1, \$3}'  \$i |  uniq > ${output_novoalign}/NM_\$filename.txt                                                       
############keep only the NM  pairs                                                  
      awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$2}' ${output_novoalign}/NM_\$filename.txt > ${output_novoalign}/NM_\${filename}_tmp.txt

      awk 'BEGIN {FS=OFS=\"\\t\"} FNR==NR{a[\$1]++;next} {if (a[\$1] > 1) print \$2, \$3}' ${output_novoalign}/NM_\${filename}_tmp.txt ${output_novoalign}/NM_\${filename}_tmp.txt > ${output_novoalign}/NM_paired_\$filename.txt
     wc -l ${output_novoalign}/NM_paired_\$filename.txt > ${output_novoalign}/NM_paired.number

     rm ${output_novoalign}/NM_\$filename.txt 
     rm ${output_novoalign}/NM_\${filename}_tmp.txt
     awk 'BEGIN {FS=OFS=\"\\t\"} {print \">\"\$1\"\n\"\$2}' ${output_novoalign}/NM_paired_\$filename.txt > ${output_novoalign}/NM_paired_\$filename.fasta                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                                     
  fi

done                      
" >> $script

	qsub   $script




    else   ######################################################### step 3 for SINGLE END data                                                                                                




################################################################# Novoalign                                                                                                                    

        echo "                                                                                                                                                                                 

#!/bin/bash -l
#$ -S /bin/bash
#$ -o $cluster_out 
#$ -e $cluster_error
#$ -l h_rt=06:00:00
#$ -pe smp 12
#$ -l tmem=1.2G,h_vmem=1.2G
#$ -N step3
#$ -wd  ${output}
#$ -V

                                                                                                                                                                                               
############## Hardclipping and quality calibration                                                                                                                                            
${novoalign} -c 12 -t180  -H -a  -k  -f  ${output_qc}/Sample${sample}.fastq -d $reference > ${output_novoalign}/${sample}.novo                                                                                                     
                                                                                                                                                                                               
" >> $script


        echo -e "5) Select all NM reads."

        echo "                                                                                                                                                                                 
novofiles="${output_novoalign}/*.novo"                                                                                                                                                         
echo "Novoalign file is '$novofiles'"                                                                                                                                                          
                                                                                                                                                                                               
                                                                                                                                                                                               
for i in \$novofiles   
do                                                                                                                                                                                             
    echo "Print '$i'"                                                                                                                                                                          
    filename=\`basename \$i .novo\`                                                                                                                                                            
    echo "Filename is '$filename'"                                                                                                                                                             
                                                                                                                                                                                               
                                                                                                                                                                                               
        awk -F\"\\t\" 'BEGIN {OFS=\"\\t\"} {if (\$5~/NM/) print \$1,\$3}'  \$i |  uniq > ${output_novoalign}/NM_\$filename.txt                                                                 
        wc -l ${output_novoalign}/NM_\$filename.txt > ${output_novoalign}/NM.number                                                                                                            
        awk '{print \">\"\$1\"\n\"\$2}' ${output_novoalign}/NM_\$filename.txt > ${output_novoalign}/NM_\$filename.fasta                                                                        
                                                                                                                                                                                               
                                                                                                                                                                                               
done                                                                                                                                                                                           
" >> $script



            qsub   $script



    fi
fi

############## end of the STEP3                                                                                                                                                               









########################################## STEP4  blastn against host - optional step                                                        

if [[ "$step4" == "TRUE" ]]; then

    
######New portion of code
    if [ ! -e $dbnr ]; then echo "File $dbnr does not exist."; exit; fi

    echo -e "12) Blastn NM reads against human_genomic"
                                                                                                                                                                                     
####################USING BLAST+

    if [[ "$paired" == "TRUE" ]]; then   ########################### step 2 for PAIRED END data                                                                                                

          echo -e "6) Submit blastn job against human reference"

        if [[ "$merged" == "TRUE" ]]; then

            blastn_input_files="${output_novoalign}/NM_${sample}_overlapping.fasta  ${output_novoalign}/NM_paired_${sample}_nonOverlapping.fasta  ${output_novoalign}/NM_${sample}_singletons.fasta"
            echo "Input files for blastn are in  '$blastn_input_files'"


        else

            blastn_input_files="${output_novoalign}/NM_paired_${sample}.fasta ${output_novoalign}/NM_${sample}_singletons.fasta"
            echo "Input files for blastn are in  '$blastn_input_files'"
        fi

	echo "                                                                                                                                                                         

#!/bin/bash -l
#$ -S /bin/bash
#$ -N step4
#$ -wd  ${output}
#$ -l h_rt=06:00:00
#$ -o $cluster_out
#$ -e $cluster_error
#$ -pe smp 12
#$ -l tmem=1.2,h_vmem=1.2G
#$ -V
" >> $script



        for i in $blastn_input_files; do


            ls -lh $i
            if [ ! -e "$i" ]; then
                echo "Error, input file for blastn $i does not exist"    ##########Check input for blastn (created in previous step) exists                                                    
                exit


            else

                filename=`basename $i .fasta`
                outputblastn=${output_blastn}/$filename.ncbiBLASTn
                echo $i
                echo $filename
                echo $output


echo "                                                                                                                                                                                               

$blastn -db $db -query  $i  -outfmt \"6 qacc sacc evalue pident qstart qend sstart send\"  -num_alignments 1   -evalue 0.1  -culling_limit 1 -num_threads 12  > $outputblastn                  
" >> $script

           fi

	done




        echo -e "7) filter out blastn hits and keep filtered dataset"

        echo "                                                                                                                                                                                                                                                                                                                                                                                
blastnfiles="${output_blastn}/*.ncbiBLASTn"                                                                                                                                                    
                                                                                                                                                                                               
for i in \$blastnfiles                                                                                                                                                                         
   do                                                                                                                                                                                          
   filename=\`basename \$i .ncbiBLASTn\`                                                                                                                                                       
       echo \$i                                                                                                                                                                                
        echo \$filename                                               
         
                                                                                                                              
        if [[ \$i == *_overlapping* ]]; then                                                                                                                                                   
####### Make filtered files for merged reads (single end)       
                                                                                                                                                                                     
            echo \$i                                                                                                                                                                           
        
                                                                                                                                                                                       
            initialNM=${output_novoalign}/NM_${sample}_overlapping.txt                                                                                                                         
            filteredNM=${output_blastn}/\${filename}_filtered.txt                                                                                                                              
            echo \$initialNM                                                                                                                                                                   
            echo \$filteredNM                                                                                                                                                                  
            

            awk -F\"\\t\" 'BEGIN {OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$1, \$2}' \${i} \${initialNM} > \${filteredNM}                                                      
            wc -l \${filteredNM} > ${output_blastn}/NM_single.number                                                                                                                           

        elif    [[ \$filename == *_singletons* ]]; then       
       
             initialNM=${output_novoalign}/NM_${sample}_singletons.txt
             filteredNM=${output_blastn}/\${filename}_filtered.txt

             awk -F\"\\t\" 'BEGIN {OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$1, \$2}' \${i} \${initialNM} > \${filteredNM}                                                                      
             wc -l \${filteredNM} > ${output_blastn}/NM_singletons.number                                                                                              
            awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' \${filteredNM} > ${output_blastn}/NM_${sample}_singletons_filtered.fasta

                                                                                                                                                                                               
        else                                                                                                                                                                                   
                                                                                                                                                                                               
####### Make fasta files for non overlapping reads (paired end)                                                                                                                                
            echo \$i                                                                                                                                                                           
                                                                                                                                                                                               
            file1=\${i}                                                                                                                                                                        
            file1tmp=${output_blastn}/\${filename}_ncbi_tmp.txt                                                                                                                                
                                                                                                                                                                                               
                                                                                                                                                                                               
            if [[ "$merged" == "TRUE" ]]; then      

                    file2=${output_novoalign}/NM_paired_${sample}_nonOverlapping.txt                                                                                                           
                    file2tmp=${output_blastn}/NM_paired_${sample}_nonOverlapping_tmp.txt                                                                                                       
          else                                                                                                                                                                                 
                   file2=${output_novoalign}/NM_paired_${sample}.txt                                                                                                                           
                    file2tmp=${output_blastn}/NM_paired_${sample}_tmp.txt                                                                                                                      
                                                                                                                                                                                               
        fi                                                                                                                                                                                     
            filtered=${output_blastn}/\${filename}_filtered.txt                                                                                                                                

                   
            awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$3}' \${file1} | uniq > \${file1tmp}
            awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$2}' \${file2} | uniq > \${file2tmp}                                                                                                                                                                                                                                                                                     
            echo \$file1tmp                                                                                                                                                                    
            echo \$file2tmp                                                                                                                                                                                                                                                                                                                                                                                 
            awk  'BEGIN {FS=OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$2,\$3}' \${file1tmp} \${file2tmp} > \${filtered}                                                     
            wc -l \${filtered} > ${output_blastn}/NM_paired.number
            awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' \${filtered} > ${output_blastn}/NM_paired_${sample}_filtered.fasta

                                                                                                                                                                                               
        fi                                                                                                                                                                                     
                                                                                                                                                                                               
    done                                                                                                                                                                                       
                                                                                                                                                                                               
" >> $script

                                                                                                                                                                                                              #rm \${file1tmp}                                                                                                                                                                    
#rm \${file2tmp}                                                                                                                                                                    

            qsub  $script


    else  ############################################################## step 4 for  SINGLE END data         

 echo -e "6) Submit blastn job against human reference"


        blastn_input_files="${output_novoalign}/NM_${sample}.fasta"
        echo "Input files for blastn are in  '$blastn_input_files'"


    for i in $blastn_input_files; do


    ls -lh $i
    if [ ! -e "$i" ]; then
        echo "Error, input file for blastn $i does not exist"    ##########Check input for blastn (created in previous step) exists                                                            
        exit


    else

        filename=`basename $i .fasta`
        outputblastn=${output_blastn}/$filename.ncbiBLASTn	
        echo $i
        echo $filename
        echo $output

        echo "
#!/bin/bash -l
#$ -S /bin/bash
#$ -N step4
#$ -wd  ${output}
#$ -l h_rt=03:00:00
#$ -o $cluster_out
#$ -e $cluster_error
#$ -pe smp 4
#$ -l tmem=6G,h_vmem=6G
#$ -V

$blastn -db $db -query  $i  -outfmt \"6 qacc sacc evalue pident qstart qend sstart send\"  -num_alignments 1   -evalue 0.1  -culling_limit 1 -num_threads 4 > $outputblastn                  
" >> $script
	
 fi
    done





    echo -e "7) filter out blastn hits and keep filtered dataset"

    echo "     

 blastnfiles="${output_blastn}/*.ncbiBLASTn"

for i in \$blastnfiles
   do
   filename=\`basename \$i .ncbiBLASTn\`
       echo \$i
        echo \$filename
####### Make filtered files for single end reads


            initialNM=${output_novoalign}/NM_${sample}.txt
            filteredNM=${output_blastn}/\${filename}_filtered.txt
            echo \$initialNM
            echo \$filteredNM
            awk -F\"\\t\" 'BEGIN {OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$1, \$2}' \${i} \${initialNM} > \${filteredNM}
            wc -l \${filteredNM} > ${output_blastn}/NM_single.number
            awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' \${filteredNM} > ${output_blastn}/${filename}_filtered.fasta

    done    

" >> $script


        qsub   $script



    fi

fi


############ end of STEP4  







########################################### Step5
#blastn against ribosomal RNA (silva)                                                        

if [[ "$step5" == "TRUE" ]]; then

    if [[ "$paired" == "TRUE" ]]; then   ########################### step 2 for PAIRED END data                                                                                                

          echo -e "6) Submit blastn job against rRNA reference"

        if [[ "$merged" == "TRUE" ]]; then

            blastn_input_files="${output_blastn}/NM_${sample}_overlapping_filtered.fasta  ${output_blastn}/NM_paired_${sample}_nonOverlapping_filtered.fasta  ${output_blastn}/NM_${sample}_singletons_filtered.fasta"
            echo "Input files for blastn are in  '$blastn_input_files'"


        else

            blastn_input_files="${output_blastn}/NM_paired_${sample}_filtered.fasta ${output_blastn}/NM_${sample}_singletons_filtered.fasta"
            echo "Input files for blastn are in  '$blastn_input_files'"
        fi




	echo "                                                                                                                                                                         
#!/bin/bash -l
#$ -S /bin/bash
#$ -N step5
#$ -wd  ${output}
#$ -l h_rt=06:00:00
#$ -o $cluster_out
#$ -e $cluster_error
#$ -pe smp 12
#$ -l tmem=1.2G,h_vmem=1.2G
#$ -V
" >> $script



        for i in $blastn_input_files; do


            ls -lh $i
            if [ ! -e "$i" ]; then
                echo "Error, input file for blastn $i does not exist"    ##########Check input for blastn (created in previous step) exists                                                    
                exit


            else

                filename=`basename $i .fasta`
                outputblastn=${output_blastnrRNA}/$filename.ncbiBLASTn
                echo $i
                echo $filename
                echo $output


echo "                                                                                                                                                                                               

$blastn -db $db -query  $i  -outfmt \"6 qacc sacc evalue pident qstart qend sstart send\"  -num_alignments 1   -evalue 0.1  -culling_limit 1 -num_threads 12  > $outputblastn                  
" >> $script

           fi

	done




        echo -e "7) filter out blastn hits and keep filtered dataset"

        echo "                                                                                                                                                                                                                                                                                                                                                                                
blastnfiles="${output_blastnrRNA}/*.ncbiBLASTn"                                                                                                                                                    
                                                                                                                                                                                               
for i in \$blastnfiles                                                                                                                                                                         
   do                                                                                                                                                                                          
   filename=\`basename \$i .ncbiBLASTn\`                                                                                                                                                       
       echo \$i                                                                                                                                                                                
        echo \$filename                                               
         
                                                                                                                              
        if [[ \$i == *_overlapping* ]]; then                                                                                                                                                   
####### Make filtered files for merged reads (single end)       
                                                                                                                                                                                     
            echo \$i                                                                                                                                                                           
        
                                                                                                                                                                                       
            initialNM=${output_blastn}/NM_${sample}_overlapping_filtered.txt                                                                                                                         
             initialNMtmp=${output_blastn}/NM_${sample}_singletons_filtered_tmp.txt
            filteredNM=${output_blastnrRNA}/\${filename}_filtered.txt                                                                                                                              
            echo \$initialNM                                                                                                                                                                   
            echo \$filteredNM                                                                                                                                                                  
            
            awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$2}' \${initialNM} | uniq > \${initialNMtmp}                                                                        
            awk -F\"\\t\" 'BEGIN {OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$2, \$3}' \${i} \${initialNMtmp} > \${filteredNM}                                                                     
            wc -l \${filteredNM} > ${output_blastnrRNA}/NM_single.number                                                                                                                           
            awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' \${filteredNM} > ${output_blastnrRNA}/NM_${sample}_singletons_filtered.fasta


        elif    [[ \$filename == *_singletons* ]]; then       
       
             initialNM=${output_blastn}/NM_${sample}_singletons_filtered.txt
             initialNMtmp=${output_blastn}/NM_${sample}_singletons_filtered_tmp.txt

             filteredNM=${output_blastnrRNA}/NM_${sample}_singletons_filtered.txt
                                                                                                                                                                                                     
#            awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$2}' \${initialNM} | uniq > \${initialNMtmp}                                                                        
             awk -F\"\\t\" 'BEGIN {OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$1, \$2}' \${i} \${initialNM} > \${filteredNM}                                                                     
             wc -l \${filteredNM} > ${output_blastnrRNA}/NM_singletons.number                                                                                              
             awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' \${filteredNM} > ${output_blastnrRNA}/NM_${sample}_singletons_filtered.fasta


                                                                                                                                                                                               
        else                                                                                                                                                                                   
                                                                                                                                                                                               
####### Make fasta files for non overlapping reads (paired end)                                                                                                                                
            echo \$i                                                                                                                                                                           
                                                                                                                                                                                               
            file1=\${i}                                                                                                                                                                        
            file1tmp=${output_blastnrRNA}/\${filename}_ncbi_tmp.txt                                                                                                                                
                                                                                                                                                                                               
                                                                                                                                                                                               
            if [[ "$merged" == "TRUE" ]]; then      

                    file2=${output_blastn}/NM_paired_${sample}_nonOverlapping_filtered.txt                                                                                                           
                    file2tmp=${output_blastnrRNA}/NM_paired_${sample}_nonOverlapping_filtered_tmp.txt                                                                                                       
          else                                                                                                                                                                                 
                   file2=${output_blastn}/NM_paired_${sample}_filtered.txt                                                                                                                           
                    file2tmp=${output_blastnrRNA}/NM_paired_${sample}_filtered_tmp.txt                                                                                                                      
                                                                                                                                                                                               
        fi                                                                                                                                                                                     
            filtered=${output_blastnrRNA}/NM_paired_${sample}_filtered.txt                                                                                                                                

                   
            awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$3}' \${file1} | uniq > \${file1tmp}
            awk  -v x=\"${delimiter}\"  'BEGIN {FS=OFS=\"\\t\"} {split(\$1, a, x);  print a[1], \$1, \$2}' \${file2} | uniq > \${file2tmp}                                                                                                                                                                                                                                                                                     
            echo \$file1tmp                                                                                                                                                                    
            echo \$file2tmp                                                                                                                                                                                                                                                                                                                                                                                 
            awk  'BEGIN {FS=OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$2,\$3}' \${file1tmp} \${file2tmp} > \${filtered}                                                     
            wc -l \${filtered} > ${output_blastnrRNA}/NM_paired_filtered.number
            awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' \${filtered} > ${output_blastnrRNA}/NM_paired_${sample}_filtered.fasta

                                                                                                                                                                                               
        fi                                                                                                                                                                                     
                                                                                                                                                                                               
    done                                                                                                                                                                                       
                                                                                                                                                                                               
" >> $script

                                                                                                                                                                                                              #rm \${file1tmp}                                                                                                                                                                    
#rm \${file2tmp}                                                                                                                                                                    

            qsub  $script


    else  ############################################################## step 5 for  SINGLE END data         

 echo -e "6) Submit blastn job against rRNA reference"


        blastn_input_files="${output_blastn}/NM_${sample}_filtered.fasta"
        echo "Input files for blastn are in  '$blastn_input_files'"


    for i in $blastn_input_files; do


    ls -lh $i
    if [ ! -e "$i" ]; then
        echo "Error, input file for blastn $i does not exist"    ##########Check input for blastn (created in previous step) exists                                                            
        exit


    else

        filename=`basename $i .fasta`
        outputblastn=${output_blastnrRNA}/$filename.ncbiBLASTn	
        echo $i
        echo $filename
        echo $output

        echo "                                                                                                                                                                                 
#!/bin/bash -l
#$ -S /bin/bash
#$ -N step5
#$ -wd  ${output}
#$ -l h_rt=06:00:00
#$ -o $cluster_out
#$ -e $cluster_error
#$ -pe smp 4
#$ -l tmem=5.9G,h_vmem=5.9G
#$ -V

$blastn -db $db -query  $i  -outfmt \"6 qacc sacc evalue pident qstart qend sstart send\"  -num_alignments 1   -evalue 0.1  -culling_limit 1 -num_threads 4  > $outputblastn                  
" >> $script
	
 fi
    done





    echo -e "7) filter out blastn hits and keep filtered dataset"

    echo "     

 blastnfiles="${output_blastnrRNA}/*.ncbiBLASTn"

for i in \$blastnfiles
   do
   filename=\`basename \$i .ncbiBLASTn\`
       echo \$i
        echo \$filename
####### Make filtered files for single end reads


            initialNM=${output_blastn}/NM_${sample}_filtered.txt
            filteredNM=${output_blastnrRNA}/NM_${sample}_filtered.txt
            echo \$initialNM
            echo \$filteredNM
            awk -F\"\\t\" 'BEGIN {OFS=\"\\t\"} NR==FNR{a[\$1]=\$1;next} a[\$1]!=\$1{print \$1, \$2}' \${i} \${initialNM} > \${filteredNM}
            wc -l \${filteredNM} > ${output_blastnrRNA}/NM_single.number
             awk -F\"\\t\" '{print \">\"\$1\"\n\"\$2}' \${filteredNM} > ${output_blastnrRNA}/NM_${sample}_filtered.fasta

    done    

" >> $script


        qsub   $script



    fi

fi


########### end of STEP5                                                                                                                                                                       







###################################################### STEP6 with queue1-Velvet

if [[ "$step6" == "TRUE" ]]; then

echo "
#!/bin/bash -l
#$ -S /bin/bash
#$ -N step6
#$ -wd  ${output}
#$ -l h_rt=12:00:00
#$ -o $cluster_out
#$ -e $cluster_error
#$ -l tmem=14G,h_vmem=14G
#$ -V
#$ -l s_stack=20M



 " >> $script                        

    if [[ "$paired" == "TRUE" ]]; then     ################################ step4 for PAIRED END data

        if [[ "$merged" == "TRUE" ]]; then
            if [[ "$blastnDone" == "TRUE" ]]; then                              ########## use filtered datasets


                velvet_input_files="${output_blastnrRNA}/NM_${sample}_overlapping_filtered_all.txt ${output_blastnrRNA}/NM_paired_${sample}_nonOverlapping_filtered_all.txt ${output_blastnrRNA}/NM_${sample}_singletons_filtered_all.txt"
                echo "Files to create input files for velvet are in '$velvet_input_files'"

            else                                                                        ######### use datasets from novoalign

                velvet_input_files="${output_novoalign}/NM_${sample}_overlapping.txt  ${output_novoalign}/NM_paired_${sample}_nonOverlapping.txt ${output_novoalign}/NM_${sample}_singletons.txt"
                echo "Files to create input files for velvet are in  '$velvet_input_files'"
            fi   ###end of blastndone  loop


        else     #### don't have merged files

           if [[ "$blastnDone" == "TRUE" ]]; then
                velvet_input_files="${output_blastnrRNA}/NM_paired_${sample}_filtered.txt ${output_blastnrRNA}/NM_${sample}_singletons_filtered.txt"
            else
                velvet_input_files="${output_novoalign}/NM_paired_${sample}.txt ${output_novoalign}/NM_${sample}_singletons.txt"
                echo "Files to create input files for velvet are in  '$velvet_input_files'"
 

           fi   ###end of blastndone loop

        fi        ####end of merged loop


       for i in $velvet_input_files; do

           ls -lh $i
           if [ ! -e "$i" ]; then

 echo "Error, input file for velvet $i does not exist"    ##########Check that input for velvet (created in previous step) exists
              exit


           else

                velvetInputPaired=${output_velvet_input}/velvet_input_paired.fasta
                velvetInputSingle=${output_velvet_input}/velvet_input_single.fasta



               if [[ $i =~ .*_paired.* ]]; then   ################## if the file consists of paired reads
                   awk -F"\t" '{print ">"$1"\n"$2}' $i >> $velvetInputPaired

               else

                   awk -F"\t" '{print ">"$1"\n"$2}' $i >> $velvetInputSingle

             fi

           fi   ####end of $i loop

       done

        echo -e "10) Velvet. Try several parameters."


        if [[ "$merged" == "TRUE" ]]; then

            velvetInputFiles="${velvetInputSingle} ${velvetInputPaired}"
            echo "Input files for velvet are in  '$velvetInputFiles'"

        

###Preprocess no HASH

minK=13


		echo "
${velveth} ${output_velvet}/output_k${minK}  ${minK} -fasta -shortPaired $velvetInputPaired -short $velvetInputSingle -noHash

" >> $script


	    for id in `seq 25 2 51`;

	    do

		outdir=${output_velvet}/output_k${id}/

		mkdir $outdir
		echo $outdir

		echo "
ln -s ${output_velvet}/output_k${minK}/Sequences 	 ${output_velvet}/output_k${id}/Sequences

${velveth}  ${outdir} $id -fasta -shortPaired $velvetInputPaired -short $velvetInputSingle -reuse_Sequences

" >> $script

	    done
	

	else
minK=13
	    velvetInputFiles="$velvetInputPaired $velvetInputSingle"
	    echo "Input files for velvet are in  '$velvetInputFiles'"

		echo "
${velveth} ${output_velvet}/output_k${minK}  ${minK} -fasta -shortPaired $velvetInputPaired -short $velvetInputSingle -noHash

" >> $script


	    for id in `seq 25 2 51`;

	    do

		outdir=${output_velvet}/output_k${id}/

		mkdir $outdir
		echo $outdir

		echo "
ln -s ${output_velvet}/output_k${minK}/Sequences 	 ${output_velvet}/output_k${id}/Sequences

${velveth}  ${outdir} $id -fasta -shortPaired $velvetInputPaired -short $velvetInputSingle -reuse_Sequences

" >> $script

	    done
	   
       
fi

	for id in `seq 25 2 51`;

	do

	    dir=${output_velvet}/output_k$id

	    echo "
${velvetg} $dir  -read_trkg yes -amos_file yes -unused_reads yes 

" >> $script

	done


	    qsub    $script


    else  #####################################step 6 for SINGLE END data

	echo -e "9) Take non human reads and make fasta file for Velvet."   

	echo -e "10) Velvet. Try several parameters."

	velvetInputSingle="${output_blastnrRNA}/NM_${sample}_filtered.fasta"

####### Make fasta files for merged reads (single end)
	         
	echo "Input files for velvet are in  '$velvetInputSingle'"



###Preprocess no HASH
      
minK=13
		echo "

${velveth} ${output_velvet}/output_k$minK  $minK -fasta -short $velvetInputSingle -noHash
" >> $script


	    for id in `seq 15 2 39`;
	    do

		outdir=${output_velvet}/output_k$id/

		mkdir $outdir
		echo $outdir

		echo "
ln -s ${output_velvet}/output_k${minK}/Sequences 	 ${output_velvet}/output_k$id/Sequences

${velveth}  $outdir $id -fasta -short $velvetInputSingle -reuse_Sequences

" >> $script

	    done		   

	    for id in `seq 15 2 39`;
	do

	    dir=${output_velvet}/output_k$id

	    echo "
${velvetg} $dir  -read_trkg yes -amos_file yes -unused_reads yes 

" >> $script

	done

   


	    qsub   $script
fi

    fi



############# end of STEP6




################################################# STEP7  - Extract contigs greater than 150bp for all k and print in a file their identifiers


if [[ "$step7" == "TRUE" ]]; then

echo "
#!/bin/bash -l
#$ -S /bin/bash
#$ -N step7
#$ -wd  ${output}
#$ -l h_rt=0:10:00
#$ -o $cluster_out
#$ -e $cluster_error
#$ -l tmem=1.2G,h_vmem=1.2G
#$ -V


" >> $script


   for k in `seq 25 2 51`; do 

#   for k in `seq 15 2 39`; do 

    inputContigs=${output_velvet}/output_k$k/contigs.fa

    echo $k
    echo "Extract Long contigs from  $inputContigs"
    
    echo "

export PERL5LIB=\${PERL5LIB}:${Software}/bioperl-live

perl $perlExtractContigs  $inputContigs ${output_velvet} k$k $contigLengthCutoff
" >> $script
	    
done

qsub     $script
 

	       
fi

################# end of STEP7






###################################  STEP8: Exctract number of reads forming each contig (need it for MCMC to adjust for contigs weights)

if [[ "$step8" == "TRUE" ]]; then


   awk '{print $0, FILENAME}'  ${output_velvet}/k* | sort -k 2 -n  | tail -1 |awk -F" " '{print $1}' | awk -F"/" '{print $NF}'> ${output_velvet}/kmer_params.txt

for k in `cat ${output_velvet}/kmer_params.txt`; do 

echo $k


        contigID="${output_velvet}/output${k}_node_ids.txt"

	awk -F"\t" '{print $1}' ${output_velvet}/output${k}_contigs_grt150.txt | awk 'BEGIN{OFS="\t"} {split($1, a, "_"); print a[2], $1}' > ${contigID}



############dokimh gia array
                wc -l $contigID > ${output_velvet}/contigs_${k}.number  

                numExtractJobs=`awk '{print $1}' ${output_velvet}/contigs_${k}.number | awk '{print $1}'`

                goodNum=$(( $numExtractJobs + 1 ))
                echo "Contigs for blastx against nr are in $i and number of jobs $goodNum"
                    extractReads=${cluster_submission}/extractReads_$k.sh
                    echo $extractReads





echo "
#!/bin/bash -l
#$ -S /bin/bash
#$ -N step8
#$ -wd  ${output}
#$ -l h_rt=02:10:00
#$ -t 1-${goodNum}
#$ -o $cluster_out
#$ -e $cluster_error
#$ -l tmem=3G,h_vmem=3G
#$ -V




# Change index from 1-2000 to 0-1999.
index=\$((\$SGE_TASK_ID-1))

id=\`awk -F\"\\t\" -v i=\$index '{if (NR==i) print \$1}' ${contigID}\`    

contigId=\`awk -F\"\\t\" -v i=\$index '{if (NR==i) print \$2}' ${contigID}\`    

export PERL5LIB=\${PERL5LIB}:${Software}/bioperl-live

perl $perlExtractReads \$id ${output_velvet}/output_$k/ > ${readsContigs}/\$contigId.fa

grep '>' ${readsContigs}/\$contigId.fa | wc -l > ${readsContigs}/\$contigId.number

" >> $extractReads


	    qsub  $extractReads

	    echo "sub script: $extractReads"
	  
	done              
#done


fi

##########end of STEP8




################### START of STEP9

if [[ "$step9" == "TRUE" ]]; then



for k in `cat ${output_velvet}/kmer_params.txt`; do 
    echo $k


    usedReads=${output_velvet}/output${k}_reads_contrib_longContigs.fa
    usedReads1Line=${output_velvet}/output${k}_reads_contrib_longContigs1Line.fa
    usedReadsTab=${output_velvet}/output${k}_reads_contrib_longContigs.tab
    contigsnumberReads=${output_velvet}/contigs_numberReads.tab
    unusedReadsTab=${output_velvet}/output${k}_reads_NOTcontrib_longContigs.tab

    if [ -e $usedReads ]; then rm $usedReads; fi
    if [ -e $usedReadsTab ]; then rm $usedReadsTab; fi
    if [ -e $contigsnumberReads ]; then rm $contigsnumberReads; fi
    if [ -e $usedReads1Line ]; then rm $usedReads1Line; fi
    if [ -e $unusedReadsTab ]; then rm $unusedReadsTab; fi

    


    find ${readsContigs}  -name \*.fa -execdir cat {} \; > $usedReads

    contigID="${output_velvet}/output${k}_node_ids.txt"

                                                                                                                                                                                             
    

    for id in  `awk -F"\t" '{print $2}' $contigID`; do
	
###print in a file the contig id plus the number of reads that have created it

	awk 'BEGIN{OFS="\t"} {nbr=split(FILENAME, a, "/");  if (i=nbr) print $0, a[i]}' ${readsContigs}/${id}.number | awk -F".number" '{print $1}'  >> ${contigsnumberReads}
 
    done
done


#####Make fasta tab-delimited

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $usedReads > $usedReads1Line   ###first put multi-line sequence in one line


awk 'BEGIN{RS=">"}{gsub("\n","\t",$0); print $0}' $usedReads1Line >> $usedReadsTab

numbReadsInfo="${output_velvet}/all_contigs${k}.number"
wc -l $usedReadsTab > ${numbReadsInfo}


if [[ "$paired" == "TRUE" ]]; then
####### Find unused reads                                                                                                                                
                          
###cat the txt files in novoalign directory so that we only look into one file
if [[ "$merged" == "TRUE" ]]; then
     
    if [[ "$blastnDone" == "TRUE" ]]; then
         
	fileOverlap=${output_blastnrRNA}/NM_${sample}_overlapping_filtered.txt
        filenonOverlap=${output_blastnrRNA}/NM_paired_${sample}_nonOverlapping_filtered.txt
	fileSingletons=${output_blastnrRNA}/NM_${sample}_singletons_filtered.txt

        file=${output_velvet}/NM_${sample}.txt

	cat $fileOverlap $filenonOverlap $fileSingletons > $file

############################ need to do this !!!! #######################     
    else

	fileOverlap=${output_novoalign}/NM_${sample}_overlapping.txt                                                                                                         
	filenonOverlap=${output_novoalign}/NM_paired_${sample}_nonOverlapping.txt                                                                                            
	fileSingletons=${output_novoalign}/NM_${sample}_singletons.txt

	file=${output_velvet}/NM_${sample}.txt
 
	cat $fileOverlap $filenonOverlap $fileSingletons > $file
    fi    #### end of blastnDone loop

else    #### not merged

    if [[ "$blastnDone" == "TRUE" ]]; then
         
	filePaired=${output_blastnrRNA}/NM_paired_${sample}_filtered.txt
        fileSingletons=${output_blastnrRNA}/NM_${sample}_singletons_filtered.txt

	file=${output_velvet}/NM_${sample}.txt

	cat ${filePaired} ${fileSingletons} > ${file}

     
    else

	filePaired=${output_novoalign}/NM_paired_${sample}.txt
        fileSingletons=${output_novoalign}/NM_${sample}_singletons.txt
	file=${output_velvet}/NM_${sample}.txt


	cat $filePaired $fileSingletons > $file

    fi   ###end of blastndone loop                                                     


fi   ###end of merged loop

          

awk -F" " 'BEGIN {FS=OFS="\t"} NR==FNR{a[$1]=$1;next} a[$1]!=$1{print $1,$2}' ${usedReadsTab} ${file} > ${unusedReadsTab}



    else  #####################################step 9 for SINGLE END data

       
        fileSingle=${output_blastnrRNA}/NM_${sample}_filtered.txt

	file=${output_velvet}/NM_${sample}.txt

	cat $fileSingle > $file
	awk -F" " 'BEGIN {FS=OFS="\t"} NR==FNR{a[$1]=$1;next} a[$1]!=$1{print $1,$2}' ${usedReadsTab} ${file} > ${unusedReadsTab}

fi


fi
##########end of STEP9



#####################################################  STEP10 - BLASTx ###############################################################
if [[ "$step10" == "TRUE" ]]; then
   
    echo -e "12) Blastx contigs against custom"

    for k in `cat ${output_velvet}/kmer_params.txt`; do 

	contigs_input=${output_velvet}/output${k}_contigs_grt150.txt


	unusedReadsTab=${output_velvet}/output${k}_reads_NOTcontrib_longContigs.tab
	blastx_input_temp=${output_velvet}/output${k}_contigs_unusedReads.tab


	cat ${contigs_input} ${unusedReadsTab} > ${blastx_input_temp}
	ls -lh ${blastx_input_temp}
	blastx_input=${output_velvet}/output${k}_contigs_unusedReads.fa
	

	awk -F"\t" '{print ">"$1"\n"$2}' ${blastx_input_temp} > ${blastx_input}

	ls -lh ${blastx_input}

	echo "Input files for blastx are in  '$blastx_input'"
	
	for i in $blastx_input; do

	    ls -lh $i
 	    if [ ! -e "$i" ]; then 
 		echo "Error, input contigs for blastx against nr $i do not exist"    ##########Check that input for blastx (created in previous step) exist
 		exit
    
	    else
		filename=`basename $i .fa`
		outputblastn=${contigs_subsetted}/$filename.tab
             
                echo $filename
                echo $outputblastn

		


		echo "                                                                                                                                                                         


#!/bin/bash -l
#$ -S /bin/bash
#$ -N step10
#$ -wd  ${output}
#$ -l h_rt=24:00:00
#$ -o $cluster_out
#$ -e $cluster_error
#$ -pe smp 12
#$ -l tmem=1.1G,h_vmem=1.1G
#$ -V


$blastx -db $dbnr -query  $blastx_input  -outfmt \"6 qacc sacc stitle score bitscore length pident evalue qstart qend sstart send\"  -word_size 3   -evalue 10  -max_target_seqs 10  -num_threads 12  > $outputblastn                  
" >> $script

		qsub $script

	    fi
	done
    done
fi   ##########end of step10

       






############ STEP11 - calculate mismatches - move it to the 1st step wuth R so that last step of bioinf pipeline is the BLAST

if [[ "$step11" == "TRUE" ]]; then

#!/bin/bash
#$ -S /bin/bash
    

    for k in `cat ${output_velvet}/kmer_params.txt`; do 
	

mismatchesRead=${output_velvet}/output${k}_peptMismatches_read.txt
mismatchesContig=${output_velvet}/output${k}_peptMismatches_contig.txt
mismatches=${output_velvet}/output${k}_peptMismatches.txt

mismatchesallInfo=${output_velvet}/output${k}_allInfo.txt
mismatchesallInfoWeights=${output_velvet}/output${k}_allInfo_weights.txt

refInfo=/home/ucbtsmo/scratchDir/sofia/sequence_database/custom_db/extractHMfromNR/gi_taxon_protLength.tab


awk 'BEGIN {FS=OFS="\t"} FNR==NR { g[$1] = $0; next } {if ($1 in g)  print  g[$1], $0 }'   ${output_velvet}/output${k}_contigs_unusedReads.tab  ${contigs_subsetted}/output${k}_contigs_unusedReads.tab |  awk 'BEGIN {FS=OFS="\t"} {print $1, $2, 10, $4, $5, $6, $7, $8, $9, $10,1,  $11, $12, $13, $14, length($2)}' > ${output_velvet}/output${k}_contigs_vs_nr.txt


awk  'BEGIN{FS=OFS="\t"} {if ($3!~/Failed_LC/) print $0, length($2)+0}' ${output_velvet}/output${k}_contigs_vs_nr.txt | awk  'BEGIN{OFS=FS="\t"} !($16+0) || !($8+0) || !($9+0) {print $0; next}{print $0, ($16/3) -(($8) * ($9))/100}' > ${mismatches}


awk 'BEGIN {FS=OFS="\t"} FNR==NR { g[$1] = $0; next } {if ($4 in g) print  $0, g[$4] }' $refInfo  $mismatches  >  $mismatchesallInfo

awk  'BEGIN {FS=OFS="\t"} NR==FNR {a[$2]=$1;next}{if ($1 in a) print $0, a[$1]; else print $0, 1}' ${output_velvet}/contigs_numberReads.tab $mismatchesallInfo > $mismatchesallInfoWeights



    done
fi
