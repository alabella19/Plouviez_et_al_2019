#!/usr/bin/perl
use Data::Dumper;

$file = $ARGV[0];

open(INPUT, $file) || die $!;
@fas_array = <INPUT>;
close (INPUT);

$new_array = ();

foreach $line (@fas_array){
	@chars = ();
	if($line !~ />/){
		@chars = split('',$line);
		$l = scalar (@chars);
		$n = 0; 
		while($n < $l){
			if($chars[$n] !~ /\-/){
				last
			}else{
				$chars[$n] = "?";
			}
			$n = $n + 1;
		}
		$n = $l;
		$n = $n-2;
		
		while($n > 0) {
			if($chars[$n] !~ /\-/){
				last;
			}else{
				$chars[$n] = "?";
			}
			$n = $n - 1; 
		}
		$newSeq = join('',@chars);
		push(@new_array,$newSeq);
	}else{
		push(@new_array, $line);
	}
}

open FILE, ">", $file;
print FILE @new_array;
close FILE;
