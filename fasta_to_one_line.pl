#!/usr/bin/perl
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

#This file turns a FASTA file into a fasta file withall sequences on one line 
if($#ARGV<0){
	print "*******************************************\nSyntax: big_gap_closer aligned_file.fasta";
	exit;
}

$file = $ARGV[0];

open (INPUT, $file);
@fasta = <INPUT>;
close(INPUT);
$first = 0;

foreach $line (@fasta){
	if($line =~ />/ && $first == 0){
		$line =~ s/\n//g;
		$line = $line."\n";
		$first = 1;
	}elsif($line =~ />/){
		$line =~ s/\n//g;
		$line = "\n".$line."\n";
	}else{
		$line =~ s/\n//g;
	}
}

	open FILE, ">" , "$file";
	print FILE @fasta;
	close FILE;
