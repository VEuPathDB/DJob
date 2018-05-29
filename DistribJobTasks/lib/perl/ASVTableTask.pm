package DJob::DistribJobTasks::ASVTableTask;

use DJob::DistribJob::Task;
use File::Basename;
use Cwd;
use CBIL::Util::Utils;
use Data::Dumper;

@ISA = (DJob::DistribJob::Task);
use strict;

#samplesInfoFile looks for columns 'NAMES' 'BARCODES' 'GROUPS' and anything else will be ignored
#if no samplesInfoFile or if 'GROUPS' column doesnt exist within it then we'll assume to run as one whole set
#should note we're assuming files in dataDir are all from the same run. if they aren't, define groups by run

# format: [name, default (or null if reqd), comment]
my @properties = 
    (
     ["dataDir",   "none",     "full path to directory of demultiplexed reads files"],
     ["multiplexed",   "false",     "default is 'false'"],
     ["barcodesType", "independant", "either 'independant' if barcodes are in a seperate file or 'within' if they are within the reads. default is 'independant'"],
     ["samplesInfoFile", "none", "path to file mapping barcodes to sample names. accepted columns: 'NAMES' 'BARCODES' 'GROUPS'. if file is not specified or GROUPS columns does not exist, will run all samples as a single group. all samples within a group must belong to the same run."],
     ["isPaired", "true", "default is 'true'"],
     ["trimLeft", "none", "number of bases to remove from beginning of forward reads"],
     ["trimLeftR", "none", "number of bases to remove from beginning of reverse reads"],
     ["truncLen", "none", "length in bases to truncate sequence to for forward reads"],
     ["truncLenR", "none", "length in bases to truncate sequence to for reverse reads"],
     ["readLen", 250, "default is 250"],
     ["taxonRefFile", "", "full path to fasta file to be used to assign taxonomy to variants"]
    );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
  my ($self, $inputDir) = @_;

  my $samplesDir = $self->getProperty('dataDir');
  if (!$samplesDir || $samplesDir eq 'none') {
    die "must provide a 'dataDir' argument";
  } elsif (! -e $samplesDir) {
    die "argument to 'dataDir' must be a valid path";
  }

  my $multiplexed = $self->getProperty('multiplexed');
  if ($multiplexed ne 'true' && $multiplexed ne 'false') {
    die "argument to 'multiplexed' must be 'true' or 'false'";
  }

  my $isPaired = $self->getProperty('isPaired');
  if ($isPaired ne 'true' && $isPaired ne 'false') {
    die "argument to 'isPaired' must be 'true' or 'false'";
  }

  my $barcodesType = $self->getProperty('barcodesType');  

  #here set forward and reverse reads and barcodes variables based on dataDir
  #then run checks for all inputs
  my @zipFiles = glob("$samplesDir/*.gz");
  foreach my $zipFile (@zipFiles) {
    print "unzipping $zipFile";
    `gunzip $zipFile`;
  }

  if (glob("$samplesDir/*.fq")) {
    print "changing samples file ext from .fq to .fastq";
    `rename *.fq *.fastq`;
  }

  my $forwardReads = 'none';
  my $forwardBarcodes = 'none';
  my $reverseReads = 'none';
  my $reverseBarcodes = 'none';

  if ($multiplexed eq 'true') {
    if ($isPaired eq 'true') {
      if ($barcodesType eq 'independant') {
        my @seqFiles = glob("$samplesDir/*sequences*");
        if (scalar @seqFiles == 2) {
          #double check perl syntax here
          $forwardReads = grep(/R1/, @seqFiles);
          $reverseReads = grep(/R2/, @seqFiles);
        } else {
          #error out
        }  
        my @barcodesFiles = glob("$samplesDir/*barcodes*");
        if (scalar @seqFiles == 2) {
          #double check perl syntax here
          $forwardBarcodes = grep(/R1/, @barcodesFiles);
          $reverseBarcodes = grep(/R2/, @barcodesFiles);
        } else {
          #error out
        }
      } else {
        my @seqFiles = glob("$samplesDir/*sequences*");
        if (scalar @seqFiles == 2) {
          $forwardReads = grep(/R1/, @seqFiles);
          $reverseReads = grep(/R2/, @seqFiles);
        } else {
          #error out
        }
      }
    } else {
      if ($barcodesType eq 'independant') {
        my @seqFiles = glob("$samplesDir/*sequences*");
        if (scalar @seqFiles == 1) {
          $forwardReads = $seqFiles[0];
        } else {
          #error out
        }
        my @barcodesFiles = glob("$samplesDir/*barcodes*");
        if (scalar @barcodesFiles == 1) {
          $forwardBarcodes = $barcodesFiles[0];
        } else {
          #error out
        }
      } else {
        my @seqFiles = glob("$samplesDir/*sequences*");
        if (scalar @seqFiles == 1) {
          $forwardReads = $seqFiles[0];
        } else {
          #error out
        }
      }
    }
  }

  my $samplesInfo = $self->getProperty('samplesInfoFile');
  if (!$samplesInfo || $samplesInfo eq 'none') {
    if ($multiplexed eq 'true') {
      die "must provide 'samplesInfoFile' argument if multiplexed is true";
    } else {
      print "WARNING: no samplesInfoFile path provided.. will be run as a single group.";
    }
  } elsif (! -e $samplesInfo) {
    die "argument to 'samplesInfoFile' must be a valid path";
  }

  my $trimLeft = $self->getProperty('trimLeft');
  my $trimLeftR = $self->getProperty('trimLeftR');
  my $truncLen = $self->getProperty('truncLen');
  my $truncLenR = $self->getProperty('truncLenR');
  my $readLen = $self->getProperty('readLen');

  if ($multiplexed eq 'true') {
    if (!$forwardReads || $forwardReads eq 'none') {
      die "must provide the 'forwardReads' argument for multiplexed reads.. did you mean to set 'multiplexed' to 'false'?";
    }
    
    if (!$barcodesType || $barcodesType eq '') {
      die "must provide the 'barcodesType' argument for multiplexed reads.. did you mean to set 'multiplexed' to 'false'?";
    } elsif ($barcodesType ne 'within' && $barcodesType ne 'independant') {
      die "argument to 'barcodesType' must be 'within' or 'independant'";
    }

    if ($barcodesType eq 'independant') {
      if (!$forwardBarcodes || $forwardBarcodes eq 'none') {
        die "must provide the 'forwardBarcodes' argument for multiplexed reads with independant barcodes.. did you mean to set 'barcodesType' to 'within'?";
      } elsif (! -e $forwardBarcodes) {
        if (-e "$forwardBarcodes.gz") {
          print "unzipping $forwardBarcodes.gz\n";
          `gunzip $forwardBarcodes.gz`;
        } else {
          die "argument to 'forwardBarcodes' must be a valid path";
        }
      }
    }

    if ($isPaired eq 'true') {
      if (!$reverseReads || $reverseReads eq 'none') {
        die "must provide the 'reverseReads' argument for multiplexed paired reads.. did you mean to set 'isPaired' to 'false'?";
      } elsif (! -e $reverseReads) {
        if (-e "$reverseReads.gz") {
          print "unzipping $reverseReads.gz\n";
          `gunzip $reverseReads.gz`;
        } else {
          die "argument to 'reverseReads' must be a valid path";
        }
      }
     
      if ($barcodesType eq 'independant') {
        if (!$reverseBarcodes || $reverseBarcodes eq 'none') {
          die "must provide the 'reverseBarcodes' argument for multiplexed paired reads with independant barcodes.. did you mean to set 'barcodesType' to 'within'?";
        } elsif (! -e $reverseBarcodes) {
          if (-e "$reverseBarcodes.gz") {
            print "unzipping $reverseBarcodes.gz\n";
            `gunzip $reverseBarcodes.gz`;
          } else {
            die "argument to 'reverseBarcodes' must be a valid path";
          }
        }
      }

      if ($truncLen eq 'none') {
        if ($truncLenR ne 'none') {
          die "truncation length was specified for reverse reads only.";
        } elsif ($readLen eq 'none') {
          die "must provide either the truncation lengths or the read length for paired reads";
        }
      }
    } else {
      if ($trimLeftR ne 'none' || $truncLenR ne 'none') {
        die "invalid arguments found. 'trimLeftR' or 'truncLenR' was specified but 'isPaired' is false.. did you mean to set 'isPaired' to 'true'?";
      }
    }
  } else {
    if (!$samplesDir || $samplesDir eq 'none') {
      die "must provide the 'dataDir' argument for demultiplexed reads.. did you mean to set 'multiplexed' to 'true'?";
    }
  }

  #need a place to demux files to.. 
  if (!$samplesDir || $samplesDir eq 'none') {
    $samplesDir = $inputDir
  }

  my $cmd = "Rscript $ENV{GUS_HOME}/bin/demuxAndBuildErrorModels.R $samplesDir $inputDir $forwardReads $reverseReads $forwardBarcodes $reverseBarcodes $multiplexed $barcodesType $samplesInfo $isPaired $trimLeft $trimLeftR $truncLen $truncLenR $readLen";



  &runCmd($cmd);

  if ($multiplexed eq 'true') {
    $samplesDir = $samplesDir . "/demux/filtered";
  } else {
    $samplesDir = $samplesDir . "/filtered";
  }
  
  $self->setProperty('dataDir', $samplesDir);

  my @fileArr;

  opendir(DIR, $samplesDir) || die "Can't open directory $samplesDir";

  while (defined (my $file = readdir (DIR))) {
    my ($name, $dir, $ext) = fileparse($file, qr/\.[^.]*/);
    next if ($file eq "." || $file eq "..");
    push(@fileArr,$file);
  }

  my $count = scalar(@fileArr);
  if ($isPaired eq 'true') {
    $count = $count / 2;
  }
  print STDERR Dumper($count);
  $self->{fileArray} = \@fileArr;
  $self->{size} = $count;

}

sub initNode {
  my ($self, $node, $inputDir) = @_;
}

sub getInputSetSize {
  my ($self, $inputDir) = @_;
}

sub initSubTask {
  my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir,$subTask) = @_;
 
  my $samplesDir = $self->getProperty("dataDir"); 
  my $samplesInfoFile = $self->getProperty("samplesInfoFile");
  my $isPaired = $self->getProperty("isPaired");
  my $isFeatureTable = 0;


  if(!$subTask->getRedoSubtask()){
    #consider possibility using $start to find my file doesnt work with rds files in this dir also
    my @files = @{$self->{fileArray}}[$start];
    foreach my $file (@files) {
      $self->runCmdOnNode($node, "cp $samplesDir/$file $serverSubTaskDir");
      if ($file eq 'featureTable.tab') {
        $isFeatureTable = 1;
      }
    }
  
    if (!$isFeatureTable) {
      if (!$samplesInfoFile || $samplesInfoFile eq 'none') {
        $self->runCmdOnNode($node, "cp $inputDir/err.rds $serverSubTaskDir"); 
      } else {
        print STDERR Dumper($files[0]);
        my ($sampleName, $dir, $ext) = fileparse($files[0]);
        open my $tempHandle, '<', $samplesInfoFile;
        my $firstLine = <$tempHandle>;
        close $tempHandle;
        my @header = split, $firstLine;
        if (grep(/'GROUPS'/, @header)) {
          if ($isPaired eq 'true') {
            $sampleName = substr($sampleName, 0, index($sampleName, '_'));
          }
          open my $fh, '<', $samplesInfoFile;
          chomp(my @lines = <$fh>);
          close $fh;
          my @match = grep(/$sampleName/, @lines);
          @match = split, $match[0];
          my ( $index ) = grep { $header[$_] eq 'GROUPS' } 0 .. $#header;
          my $group = $match[$index];
          my $errFile = $inputDir . "/" . $group . "_err.rds";
          $self->runCmdOnNode($node, "cp $errFile $serverSubTaskDir");
        } else {
          $self->runCmdOnNode($node, "cp $inputDir/err.rds $serverSubTaskDir");
        } 
      }
    }
  }
  $self->runCmdOnNode($node, "cp -r $serverSubTaskDir/* $nodeExecDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir,$subtaskNumber,$mainResultDir) = @_;

    my $isPaired = $self->getProperty("isPaired");

    #run dada on remaining sample(s). 
    my $cmd = "Rscript $ENV{GUS_HOME}/bin/runDada.R $nodeExecDir $isPaired";

    return($cmd)
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    #just trying to make a file with generic name on node to be specific to our node/task on server
    $self->runCmdOnNode($node, "cp $nodeExecDir/featureTable.tab $mainResultDir/${subTaskNum}_featureTable.tab");

    return $node->getErr();
}

##cleanup masterDir here and remove extra files that don't want to transfer back to compute node
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;

    my $samplesDir = $self->getProperty("dataDir");
    my $taxonRefFile = $self->getProperty("taxonRefFile");
 
    #will combine all output feature tables and assign taxonomy to the resulting final file
    my $cmd = "Rscript $ENV{GUS_HOME}/bin/merge.R $mainResultDir $taxonRefFile";

    $self->runCmdOnNode($node, $cmd); 
}

1;
