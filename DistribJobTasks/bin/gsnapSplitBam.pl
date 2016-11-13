#!/usr/bin/perl

use strict;

use Getopt::Long;
use lib "$ENV{GUS_HOME}/lib/perl";

use CBIL::TranscriptExpression::SplitBamUniqueNonUnique qw(splitBamUniqueNonUnique);

my ($mainResultDir, $strandSpecific, $isPairedEnd, $bamFile);

&GetOptions("mainResultDir=s"=> \$mainResultDir,
            "strandSpecific=s" => \$strandSpecific,
            "isPairedEnd=s" => \$isPairedEnd,
            "bamFile=s" => \$bamFile,
    );


my $splitExpDir = splitBamUniqueNonUnique($mainResultDir, $strandSpecific, $isPairedEnd, $bamFile);

1;
