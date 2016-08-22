package DJob::DistribJob::PbsNode;

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
    }elsif($runFile =~ /cancel/){
      $runFile =~ s/cancel/run/;
    }else{
      $runFile = "$runFile.run";
    }
    if(!-e "$runFile"){
      open(R,">$runFile") || die "unable to create script file '$runFile'\n";
      print R "#!/bin/sh\n\n$ENV{GUS_HOME}/bin/nodeSocketServer.pl $self->{serverHost} $self->{serverPort}\n";
      close R;
      system("chmod +x $runFile");
    }
    my $qsubcmd = "qsub -N DistribJob -V -j oe -l nodes=1:ppn=$self->{procsPerNode}".($self->{runTime} ? ",walltime=00:$self->{runTime}:00" : "").($self->{queue} ? " -q $self->{queue}" : "")." $runFile";
##    print STDERR "\n$qsubcmd\n\n";
    my $jid = `$qsubcmd`;
    chomp $jid;
    $self->{workingDir} = "/$self->{nodeWorkingDirsHome}/$jid";
    $self->setJobid($jid);
    if($self->{fileName}){
      open(C,">>$self->{fileName}");
      print C "$self->{jobid} ";
      close C;
    }
  } 
  $self->setState($QUEUED);
}

sub getNodeAddress {
  my $self = shift;
  if (!defined $self->{nodeNum}) {
    my $getCmd = "qstat -n $self->{jobid}";
    my @stat = `$getCmd`;
    return undef if $?;         ##command failed
    if ($stat[-1] =~ /^\s*(\w*\d+)/) {
      $self->{nodeNum} = $1;
    }
  }
  return $self->{nodeNum};
}

# remove this node's job from the queue
sub removeFromQueue {
  my ($self) = @_;
  my $cmd = "qdel $self->{jobid} > /dev/null 2>&1";
  system($cmd) && print STDERR "Failed running command to delete job from queue.  '$cmd'";
}

# the code in this method is suspect.  copied some old code from the old cleanUp method, before it was lost
sub deleteLogFilesAndTmpDir {
  my $self = shift;
  ##delete those pesky OU files that don't do anything
  my $jid = $self->{jobid};
  if($self->{jobid} =~ /^(\d+)/){
    $jid = $1;
  }
  my $delCmd = "/bin/rm $ENV{HOME}/DistribJob.o$jid > /dev/null 2>&1";
  system($delCmd); 
}



1;
