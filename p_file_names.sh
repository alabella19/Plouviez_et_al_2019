#!/bin/bash

firstDir=${PWD}
for i in *TM*/
do
	echo $i
	cd $i
	cd ${i%/}.assembly
	for f in T_*
	do
		mv $f ${f/T_/}
	done
#	for f in 3process*
#	do
#		mv $f ${f/_T/}
#	done
#	
	cd $firstDir

done
#for i in *ABE*/
#do
#	echo $i
#	cd $i
#	cd ${i%/}.assembly
#	for f in T_*
#	do
#		mv $f ${f/T_/}
#	done
##	for f in 3process*
##	do
##		mv $f ${f/_T/}
##	done
##	
#	cd $firstDir
#
#done
#
#for i in *SW1*/
#do
#	echo $i
#	cd $i
#	cd ${i%/}.assembly
#	for f in T_*
#	do
#		mv $f ${f/T_/}
#	done
##	for f in 3process*
##	do
##		mv $f ${f/_T/}
##	done
##	
#	cd $firstDir
#
#done
#for i in *SW8*/
#do
#	echo $i
#	cd $i
#	cd ${i%/}.assembly
#	for f in T_*
#	do
#		mv $f ${f/T_/}
#	done
##	for f in 3process*
##	do
##		mv $f ${f/_T/}
##	done
##	
#	cd $firstDir
#done