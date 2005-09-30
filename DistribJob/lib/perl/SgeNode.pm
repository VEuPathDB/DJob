package DJob::DistribJob::SgeNode;

use DJob::DistribJob::Node ":states";
use Cwd;
use strict;

our @ISA = qw(DJob::DistribJob::Node);

sub queueNode {
  my $self = shift;
  if (!$self->getJobid()) {     ##need to run qsub 
    ##first create the script...
    my $runFile = $self->{fileName};
    if(!$runFile){
      $runFile = "nodeScript.$$";
    }elsif($self->{filename} =~ /cancel/){
      $runFile =~ s/cancel/run/;
    }else{
      $runFile = "$runFile.run";
    }
    $self->{script} = $runFile;
    if(!-e "$runFile"){
      open(R,">$runFile") || die "unable to create script file '$runFile'\n";
      print R "#!/bin/sh\n\n$ENV{GUS_HOME}/bin/nodeSocketServer.pl $self->{serverHost} $self->{serverPort}\n";
      close R;
      system("chmod +x $runFile");
    }
    my $qsubcmd = "qsub -V -cwd -pe smp 2 $runFile";
    my $tjid = `$qsubcmd`;
    if($tjid =~ /job\s(\d+)/){
      my $jid = $1;
      $self->setJobid($jid);
      open(C,">>$self->{fileName}");
      print C "$self->{jobid} ";
      close C;
    }else{
      die "unable to determine jobid from $tjid\n";
    }
  } 
  $self->setState($QUEUED);
}

sub getNodeAddress {
  my $self = shift;
  if (!defined $self->{nodeNum}) {
    my $getCmd = "qstat";
    my @stat = `$getCmd`;
    return undef if $?;         ##command failed
    foreach my $line (@stat){
      if($line =~ /^\s*$self->{jobid}.*(node\d+).q/){
        $self->{nodeNum} = $1;
      }
    }
  }
  return $self->{nodeNum};
}

1;
