#!/bin/bash

firstDir=${PWD}
for i in *TM*/
do
	echo $i
	cp consensus_haplotype_V5.pl $i${i%/}_results
	cp change_process_to_amb.pl $i${i%/}_results
	cd $i
	cd ${i%/}_results
	for f in process*
	do
		#echo $f
		fas=${f#process}
		fas=${fas%.haps.nex.*}
		nex=$fas
		fas=$fas".fasta"
		nex=$nex."haps.nex"
		
		echo $fas
		perl change_process_to_amb.pl $f
		perl consensus_haplotype_V5.pl -proc_file "$f" -fas_file "$fas" -indv "$i" -nex_file "$nex"
	done
	for f in 3process*
	do
		#echo $f
		fas=${f#3process}
		fas=${fas%.haps.nex.*}
		nex=$fas
		fas=$fas".fasta"
		nex=$nex."haps.nex"
		echo $fas
		perl change_process_to_amb.pl $f
		perl consensus_haplotype_V5.pl -proc_file "$f" -fas_file "$fas" -indv "$i" -nex_file "$nex"
	done
	#rm consensus_haplotype_V5.pl
	cd $firstDir
	
done
#for i in *ABE*/
#do
#	echo $i
#	cp consensus_haplotype_V5.pl $i${i%/}_results
#	cp change_process_to_amb.pl $i${i%/}_results
#	cd $i
#	cd ${i%/}_results
#	for f in process*
#	do
#		#echo $f
#		fas=${f#process}
#		fas=${fas%.haps.nex.*}
#		nex=$fas
#		fas=$fas".fasta"
#		nex=$nex."haps.nex"
#		
#		echo $fas
#		perl change_process_to_amb.pl $f
#		perl consensus_haplotype_V5.pl -proc_file "$f" -fas_file "$fas" -indv "$i" -nex_file "$nex"
#	done
#	for f in 3process*
#	do
#		#echo $f
#		fas=${f#3process}
#		fas=${fas%.haps.nex.*}
#		nex=$fas
#		fas=$fas".fasta"
#		nex=$nex."haps.nex"
#		echo $fas
#		perl change_process_to_amb.pl $f
#		perl consensus_haplotype_V5.pl -proc_file "$f" -fas_file "$fas" -indv "$i" -nex_file "$nex"
#	done
#	#rm consensus_haplotype_V5.pl
#	cd $firstDir
#	
#done
#for i in *SW1*/
#do
#	echo $i
#	cp consensus_haplotype_V5.pl $i${i%/}_results
#	cp change_process_to_amb.pl $i${i%/}_results
#	cd $i
#	cd ${i%/}_results
#	for f in process*
#	do
#		#echo $f
#		fas=${f#process}
#		fas=${fas%.haps.nex.*}
#		nex=$fas
#		fas=$fas".fasta"
#		nex=$nex."haps.nex"
#		
#		echo $fas
#		perl change_process_to_amb.pl $f
#		perl consensus_haplotype_V5.pl -proc_file "$f" -fas_file "$fas" -indv "$i" -nex_file "$nex"
#	done
#	for f in 3process*
#	do
#		#echo $f
#		fas=${f#3process}
#		fas=${fas%.haps.nex.*}
#		nex=$fas
#		fas=$fas".fasta"
#		nex=$nex."haps.nex"
#		echo $fas
#		perl change_process_to_amb.pl $f
#		perl consensus_haplotype_V5.pl -proc_file "$f" -fas_file "$fas" -indv "$i" -nex_file "$nex"
#	done
#	#rm consensus_haplotype_V5.pl
#	cd $firstDir
#	
#done
#for i in *SW8*/
#do
#	echo $i
#	
#	cp consensus_haplotype_V5.pl $i${i%/}_results
#	cp change_process_to_amb.pl $i${i%/}_results
#	cd $i
#	cd ${i%/}_results
#	for f in process*
#	do
#		#echo $f
#		fas=${f#process}
#		fas=${fas%.haps.nex.*}
#		nex=$fas
#		fas=$fas".fasta"
#		nex=$nex."haps.nex"
#		
#		echo $fas
#		perl change_process_to_amb.pl $f
#		perl consensus_haplotype_V5.pl -proc_file "$f" -fas_file "$fas" -indv "$i" -nex_file "$nex"
#	done
#	for f in 3process*
#	do
#		#echo $f
#		fas=${f#3process}
#		fas=${fas%.haps.nex.*}
#		nex=$fas
#		fas=$fas".fasta"
#		nex=$nex."haps.nex"
#		echo $fas
#		perl change_process_to_amb.pl $f
#		perl consensus_haplotype_V5.pl -proc_file "$f" -fas_file "$fas" -indv "$i" -nex_file "$nex"
#	done
#	#rm consensus_haplotype_V5.pl
#	cd $firstDir
#	
#done
#