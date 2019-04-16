#!/usr/bin/perl
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

#This file turns a FASTA file into a fasta file withall sequences on one line 
if($#ARGV<0){
	print "*******************************************\nSyntax: split_fasta.pl IMGC_output.fasta";
	exit;
}

$file = $ARGV[0];

open (INPUT, $file);
@fasta = <INPUT>;
close(INPUT);
@filename = split('\.',$file);
$filename = $filename[0];
$block = 0;
$block_num = 0;
@current_block =();
foreach $line (@fasta){
	if($line !~ /^\n/ && $block == 0){
		push(@current_block, $line);
	}elsif($line =~ /^\n/ && $block ==0){
		$block_num = $block_num + 1;
	}elsif($block == 1 && $line =~ />/){
		#PRINT previous block 
		$prev_block_num = $block_num - 1;
		open OUTPUT, '>', "$filename.b$prev_block_num.fasta";
		print OUTPUT @current_block; 
		close OUTPUT;
		#reset the current block
		@current_block
	}
}

	open FILE, ">" , "$file";
	print FILE @fasta;
	close FILE;
