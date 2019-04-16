#!/usr/bin/perl
use Cwd;
use Data::Dumper;
$original_dir = getcwd;
@dir_contents = grep -d, <$original_dir/*>;
%loci_hash = ();

foreach $indv_dir(@dir_contents){
	@assembly_contents = grep -d, <$indv_dir/*>;
	#Thes two lines to get into the assemlby file 
	$assembly = @assembly_contents[0];
	@indv_contents = grep -f, <$assembly/*>;
	
	foreach $loci_file (@indv_contents){
		print $loci_file;
		if($loci_file =~ /\.fasta/ && $loci_file =~ /haps/){
			#print "$loci_file\n";
			open(INPUT, $loci_file);
			@loci_array = <INPUT>;
			close(INPUT);
			
			@name = split('/',$loci_file);
			$name = pop @name; 
			#print "$name\n";
			if( $name =~ /\.haps\..\.fasta/){
				@hap_num = split('\.',$name);
				$hap_num = $hap_num[2];
				#print "hap num: $hap_num\n";
			}
			
			$name=~ s/\.haps\.fasta//;
			$name =~ s/\.haps\..\.fasta//;
			
			
			if($name =~ /^3/){
				
				$name =~ s/^3//;
				
			}
			#print "final name: $name\n";
			
			
			#print "$name\n";
			
			open(INPUT, $loci_file);
			@loci_array = <INPUT>;
			close(INPUT);
			##print @loci_array;
			
			if($hap_num >1){
				foreach $line (@loci_array){
					if($line =~ />/){
						$line =~ s/\n//;
						$line = $line.".b".$hap_num."\n";
					}
				}
				
			}
			
			if($loci_hash{$name} != 1){
				##this is a new loci
				open FILE, ">", "$original_dir/$name.fasta";
				print FILE "\n";
				print FILE @loci_array;
				close FILE;
			}else{
				open FILE, ">>", "$original_dir/$name.fasta";
				print FILE "\n";
				print FILE @loci_array;
				close FILE;
			}
			
			$loci_hash{$name} = 1;
		}
	}
}

foreach $key (keys %loci_hash){
	@new_loci=();
	open(INPUT, "$original_dir/$key.fasta");
	@loci_array = <INPUT>;
	close(INPUT);
	foreach $line (@loci_array){
		if ($line =~ /^$/){
		}else{
			push(@new_loci, $line);
		}
	}
	open FILE, ">", "$original_dir/$key.fasta";
	print FILE @new_loci;
	close FILE;
}
			
	
