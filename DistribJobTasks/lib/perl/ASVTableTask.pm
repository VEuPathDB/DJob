package DJob::DistribJobTasks::ASVTableTask;

use DJob::DistribJob::Task;
use File::Basename;
use File::Path qw/make_path/;
use File::Copy qw/move/;
use Cwd;
use CBIL::Util::Utils;
use Data::Dumper;

@ISA = (DJob::DistribJob::Task);
use strict;

#samplesInfoFile looks for columns 'NAMES',  'GROUPS' and anything else will be ignored
#if no samplesInfoFile or if 'GROUPS' column doesnt exist within it then we'll assume to run as one whole set
#should note we're assuming files in dataDir are all from the same run. if they aren't, define groups by run

# format: [name, default (or null if reqd), comment]
# "none" because a blank default means required
my @properties = 
    (
     ["sraStudyId", "none", "an SRA studyId to use to retrieve all related samples"],
     ["sraSampleAndRunIdsPath", "none", "a file with sample then run IDs to retrieve from SRA, and then name as sample.run"],
     ["dataDir",   "none",     "full path to directory of demultiplexed reads files - alternative to providing a SRA study ID"],
     ["workingDir",   "",     "dir for temporary files"],
     ["platform",  "illumina",  "valid options: 'illumina', '454'"],
     ["samplesInfoFile", "none", "path to TSV file with groups and sample names. Accepted columns: 'NAMES', 'GROUPS'. If present, error models will be constructed and used within groups."],
     ["isPaired", "", "paired library layout as opposed to single"],
     ["trimLeft", "none", "number of bases to remove from beginning of forward reads"],
     ["trimLeftR", "none", "number of bases to remove from beginning of reverse reads"],
     ["truncLen", "none", "length in bases to truncate sequence to for forward reads"],
     ["truncLenR", "none", "length in bases to truncate sequence to for reverse reads"],
     ["maxLen", "none", "max len for 454 reads"],
     ["mergeTechReps", "false" , "default is 'false'"],
     ["trainingSetFile", "", "training set for taxonomic assignment"],
     ["speciesAssignmentFile", "", "reference sequences for species, matched by identity"],
     ["resultFile", "output.tsv", "where to save the final results"],
    );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
  my ($self, $inputDir) = @_;

  my $platform = $self->getProperty('platform');

  (my $samplesDir = $self->getProperty('dataDir')) =~s/^none$//;
  (my $sraStudyId = $self->getProperty('sraStudyId')) =~s/^none$//;
  (my $sraSampleAndRunIdsPath = $self->getProperty('sraSampleAndRunIdsPath')) =~s/^none$//;
  
  die "Must provide either --sraStudyId or --dataDir or --sraSampleAndRunIdsPath"
   unless 1 == grep {$_} ($samplesDir, $sraStudyId, $sraSampleAndRunIdsPath);
  

  my $workingDir = $self->getProperty('workingDir');
  make_path $workingDir; 


  my $isPaired = $self->getProperty('isPaired');
  if (lc($platform) eq 'illumina') {
    if ($isPaired ne 'true' && $isPaired ne 'false') {
      die "argument to 'isPaired' must be 'true' or 'false'";
    }
  } elsif ($platform eq '454') {
    $isPaired = 'false'; 
  } else {
    die "currently only 'illumina' and '454' are supported platforms.";
  } 

  (my $samplesInfo = $self->getProperty('samplesInfoFile'))=~s/^none$//;

  die "File not found: $samplesInfo"
    if $samplesInfo and not -f $samplesInfo;

  (my $trimLeft = $self->getProperty('trimLeft')) =~s/^none$//;
  (my $trimLeftR = $self->getProperty('trimLeftR')) =~s/^none$//;
  (my $truncLen = $self->getProperty('truncLen')) =~s/^none$//;
  (my $truncLenR = $self->getProperty('truncLenR')) =~s/^none$//;
  (my $maxLen = $self->getProperty('maxLen')) =~s/^none$//;

  if (not $maxLen and $platform eq '454') {
    die "must provide maximum read length as parameter 'maxLen' for 454 data.";
  } 


  my $mergeTechReps = $self->getProperty('mergeTechReps');

  my $fastqsInDir; 
  if ($samplesDir){
    my @zipFiles = glob("$samplesDir/*.gz");
    foreach my $zipFile (@zipFiles) {
      print STDERR "unzipping $zipFile\n";
      `gunzip $zipFile`;
    }
  
    if (glob("$samplesDir/*.fq")) {
      print STDERR "changing samples file ext from .fq to .fastq\n";
      `rename .fq .fastq $samplesDir/*.fq`;
    }
    die "$samplesDir contains no fastqs"
      unless glob("$samplesDir/*fastq");
     
    if (glob("$samplesDir/*_R1_001.fastq")){
      `rename _R1_001.fastq _1.fastq $samplesDir/*_R1_001.fastq`;
    }
    if (glob("$samplesDir/*_R2_001.fastq")){
      `rename _R2_001.fastq _2.fastq $samplesDir/*_R2_001.fastq`;
    }

    if($isPaired eq 'true'){
       die "No forward fastqs" unless glob("$samplesDir/*_1.fastq");
       die "No reverse fastqs" unless glob("$samplesDir/*_2.fastq");
    }
    $fastqsInDir = $samplesDir;
  } elsif ($sraStudyId) {
    make_path "$workingDir/fastqsFromSra";
    &runCmd("getFastqFromSra.pl --workingDir $workingDir/fastqsFromSra --studyId '$sraStudyId' --pairs $isPaired"); 
    $fastqsInDir = "$workingDir/fastqsFromSra";
  } elsif($sraSampleAndRunIdsPath){
    make_path "$workingDir/fastqsFromSra";
    &runCmd("getFastqFromSra.pl --workingDir $workingDir/fastqsFromSra --sampleAndRunIdsPath '$sraSampleAndRunIdsPath' --pairs $isPaired");
    $fastqsInDir = "$workingDir/fastqsFromSra";
  }

  if (! -e "$workingDir/filtered"){
  # Takes a while - a thousand samples per day
    my $cmd = <<"EOF";
Rscript $ENV{GUS_HOME}/bin/dada2/filterFastqs.R 
--fastqsInDir $fastqsInDir
--fastqsOutDir $workingDir/filtered.tmp
--isPaired $isPaired 
--maxLen $maxLen 
--platform $platform 
EOF
    $cmd =~ s/\n/ /g;
    $cmd .= " --samplesInfo $samplesInfo" if $samplesInfo;
    $cmd .= " --trimLeft $trimLeft" if $trimLeft;
    $cmd .= " --trimLeftR $trimLeftR" if $trimLeftR;
    $cmd .= " --truncLen $truncLen" if $truncLen;
    $cmd .= " --truncLenR $truncLenR" if $truncLenR;

    &runCmd($cmd);
    die "Failed: no fastqs in output filtering dir $workingDir/filtered.tmp"
       unless glob("$workingDir/filtered.tmp/*.fastq");

    move("$workingDir/filtered.tmp", "$workingDir/filtered");
  }
  if (! glob("$workingDir/*err.rds")){
  # This has ran out of memory in the past. Submit with more if needed.
    my $cmd = <<"EOF";
Rscript $ENV{GUS_HOME}/bin/dada2/buildErrorModels.R 
--fastqsInDir $workingDir/filtered
--errorsOutDir $workingDir 
--errorsFileNameSuffix err.rds
--isPaired $isPaired 
--platform $platform 
EOF
    $cmd =~ s/\n/ /g;
    $cmd .= " --samplesInfo $samplesInfo" if $samplesInfo;
    &runCmd($cmd);
  }

  # These naughty files triggered this bug, fixed in Jan 2020: https://github.com/Bioconductor/ShortRead/pull/2
  # If you're reading this in year 2022 verify Bioconductor has released sometime in 2020 or 2021, and that you have the good ShortRead locally
  # Then delete this code and the message
  if ($sraStudyId eq 'ERP005806'){
    unlink("$workingDir/filtered/ERS454101.ERR503975_filt.fastq");
    unlink("$workingDir/filtered/ERS453453.ERR503327_filt.fastq");
  }

  my @fileArr;

  opendir(DIR, "$workingDir/filtered") || die "$!: Can't open directory $workingDir/filtered";

  while (defined (my $file = readdir (DIR))) {
    my ($name, $dir, $ext) = fileparse($file, qr/\.[^.]*/);
    next if ($file eq "." || $file eq "..");
    next if $file =~ /^.*_2.fastq$/; # don't create jobs twice
    next if ($ext ne '.fastq'); 
    push(@fileArr,$file);
  }

  my $count = scalar(@fileArr);
  print STDERR "Created $count subtasks\n";
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
 
  my $workingDir = $self->getProperty("workingDir");
  my $samplesInfoFile = $self->getProperty("samplesInfoFile");
  my $isPaired = $self->getProperty("isPaired");
  my $isFeatureTable = 0;


  if(!$subTask->getRedoSubtask()){
    $self->runCmdOnNode($node, "mkdir -p $serverSubTaskDir");
    my $file = @{$self->{fileArray}}[$start];
    if ($isPaired eq 'true'){
      # copy both forward and reverse reads
      my ($pfx, @__) = split /_/, $file;
      die "$file doesn't end with fastq, maybe go back to this code"
        unless $file =~ /.fastq$/;
      $self->runCmdOnNode($node, "cp $workingDir/filtered/${pfx}*.fastq $serverSubTaskDir");
    } else {
      $self->runCmdOnNode($node, "cp $workingDir/filtered/${file} $serverSubTaskDir");
    }

    my $errorsFromPath;
    if (!$samplesInfoFile) {
      $errorsFromPath = "$workingDir/err.rds"; 
    } else {
      my ($sampleName, $dir, $ext) = fileparse($file);
      open my $tempHandle, '<', $samplesInfoFile;
      my $firstLine = <$tempHandle>;
      close $tempHandle;
      chomp($firstLine);
      my @header = split /\t/, $firstLine;
      my %headerHash = map { $_ => 1 } @header;
      if (exists($headerHash{GROUPS})) { 
        if ($isPaired eq 'true') {
          # In general not true for single end manual delivery files
          my ($pfx, @__) = split /_/, $file;
          $sampleName = $pfx;
        }
        open my $fh, '<', $samplesInfoFile;
        chomp(my @lines = <$fh>);
        close $fh;
        my @match = grep(/$sampleName/, @lines);
        @match = split /\t/, $match[0];
        my ( $index ) = grep { $header[$_] eq 'GROUPS' } 0 .. $#header;
        my $group = $match[$index];
        $errorsFromPath = "$workingDir/${group}_err.rds";
      } else {
        $errorsFromPath = "$workingDir/err.rds"; 
      } 
    }
    $self->runCmdOnNode($node, "cp $errorsFromPath $serverSubTaskDir/error-models.rds");
  }
  $self->runCmdOnNode($node, "cp -r $serverSubTaskDir/* $nodeExecDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir,$subtaskNumber,$mainResultDir) = @_;

    my $isPaired = $self->getProperty("isPaired");
    my $platform = $self->getProperty("platform");
    my $mergeTechReps = $self->getProperty('mergeTechReps');

    my $cmd = <<EOF;
Rscript $ENV{GUS_HOME}/bin/dada2/fastqToAsv.R
--fastqsInDir $nodeExecDir
--errorsRdsPath $nodeExecDir/error-models.rds
--outRdsPath $nodeExecDir/featureTable.rds
--isPaired $isPaired
--platform $platform
--mergeTechReps $mergeTechReps
EOF
  $cmd =~ s/\n/ /g;
    return($cmd)
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    #just trying to make a file with generic name on node to be specific to our node/task on server
    $self->runCmdOnNode($node, "cp $nodeExecDir/featureTable.rds $mainResultDir/${subTaskNum}_featureTable.rds");

    return $node->getErr();
}

sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;

    my $trainingSetFile = $self->getProperty('trainingSetFile');
    my $speciesAssignmentFile = $self->getProperty('speciesAssignmentFile');
    my $resultFile = $self->getProperty('resultFile');

    my $cmd = <<EOF; 
Rscript $ENV{GUS_HOME}/bin/dada2/mergeAsvsAndAssignToOtus.R
--asvRdsInDir $mainResultDir
--assignTaxonomyRefPath $trainingSetFile
--addSpeciesRefPath $speciesAssignmentFile
--outPath $resultFile
EOF
    $cmd =~ s/\n/ /g;
    unless (-s $resultFile){
       $self->runCmdOnNode($node, $cmd); 
    }
}

1;
