#!/usr/local/bin/perl
#
# This script is to split a set of MGF files into multiple splices depending on their sizes.
# The limit is set for 1 million (or till the next "End Ions" tag is found) lines/mgf file.
#
# @author - Ritesh Krishna 

use strict;

my $LIMIT = 1000000;  
my $endTag = "END IONS";

my $numArgs = $#ARGV  + 1;

# print usage information if no arguments supplied  
if($numArgs != 2) {
   print "Usage:\n perl mgfSplitter.pl inputMgfDirectory outputMgfDirectory\n";
   exit 1;
}

# get arguments and perform basic error checking
my($inputMgfDirectory);
  
my $inputMgfDirectory = $ARGV[0];
my $outputMgfDirectory = $ARGV[1];

mkdir($outputMgfDirectory, 0777) || die "can't create output dir: $!";

opendir(DIR, $inputMgfDirectory) || die "can't opendir $inputMgfDirectory: $!"; 

my $fileCounter = 0;

my @files = grep(/\.mgf$/, readdir(DIR));

foreach my $file (@files) { 
	
	open (INFILE,  "<$inputMgfDirectory/$file")  || die "Error: cannot open $file";
	
	my $chunk = '';
	my $lineCounter = 0;

	foreach my $line (<INFILE>) {
		#chomp($line);
		$lineCounter++;
		
		if(index($line,$endTag) >= 0){ 
			print "$lineCounter \t $line \t $endTag\n";
		}
		
		$chunk = $chunk . $line;
		
		if( $lineCounter >= $LIMIT){ 
			if (index($line,$endTag) >= 0){
				$lineCounter = 0;
				$fileCounter++;
				my $outputFile = $fileCounter . "_$file"; 
				open (OUTFILE, ">$outputMgfDirectory/$outputFile")  || die "Error: cannot open output file";
				print OUTFILE $chunk;
				close($outputFile);
				$chunk = '';
			}
		} 
	}
	
	# The last remaining bit which might be less than the LIMIT...
	if($lineCounter > 0){
		$fileCounter++;
		my $outputFile = $fileCounter . "_$file"; 
		open (OUTFILE, ">$outputMgfDirectory/$outputFile")  || die "Error: cannot open output file";
		print OUTFILE $chunk;
		close($outputFile);
	}	
	
	close(INFILE);
 	 
} 

closedir(DIR);
