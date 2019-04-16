#!/usr/bin/perl
if($#ARGV<0){
	print "*******************************************\nSyntax: split_frags.pl input.fasta";
	exit;
}



$file = $ARGV[0];
open (INPUT, $file);
@fasta = <INPUT>;
close(INPUT);

@frags = ();
@not_frags = ();
$this_read;

foreach $line (@fasta){
	if($line =~ />/){
		$this_read = $line;
	}else{
		$line_test = $line;
		chomp($line_test);
		$line_rev = reverse $line_test;
		chomp($line_rev);
		$front = chop($line_test);
		$end = chop ($line_rev);
		#print ("Testing $front and $end\n");
		
		if($front eq "?" || $end eq "?"){
			push(@frags, $this_read);
			push(@frags, $line);
		}else{
			push(@not_frags, $this_read);
			push(@not_frags, $line);
		}
	}
}

open FILE, ">" , "$file";
print FILE @not_frags;
close FILE;

#$f_z = -z $file;
#print "sites is $f_z\n";

open FILE, ">" , "frags_$file";
print FILE @frags;
close FILE;

#$f_z = -z "frags_$file";
#print "frags is $f_z\n";