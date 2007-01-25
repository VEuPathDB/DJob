#!@perl@

use strict;
use IO::Socket;
use IPC::Open2;

$| = 1;

my $serverHost = shift;
my $serverPort = shift;

my $endMatchString = 'FooCmdEnd';

my $host = `hostname -s`;
chomp $host;

##get the jobid
my $jobid;
if($ENV{PBS_JOBID}){  ##pbs 
  $jobid = $ENV{PBS_JOBID};
}elsif($ENV{JOB_ID}){  ##sge
  $jobid = $ENV{JOB_ID};
}elsif($ENV{LSB_JOBID}){ ##LSF
  $jobid = $ENV{LSB_JOBID};
}else{
  $jobid = shift;  ##local node
}
die "ERROR:  nodeSocketServer.pl is unable to determine the JOBID on $host\n" unless $jobid;
print STDERR "jobid: $jobid\n";

##now inform the server that this node is ready to run...
my $sock;
my $localPort;
my $numTries = 0;
do {
  die "ERROR: nodeSocketServer.pl serverHost=$serverHost, serverPort=$serverPort\nCould not create socket on node: $!\n"  if $numTries++ > 5;
  $localPort = int(rand(3000)) + 5000;
  $sock = new IO::Socket::INET (
                                LocalHost => $host,
                                LocalPort => $localPort,
                                Proto => 'tcp',
                                Listen => 5,
                                Reuse => 1,
                               );
}until ( $sock );

my $hostSock;
my $ct = 0;
until($hostSock){
  $hostSock = new IO::Socket::INET (
                                PeerAddr => $serverHost,
                                PeerPort => $serverPort,
                                Proto => 'tcp',
                               );
  unless($hostSock){
    die "ERROR: nodeSocketServer.pl serverHost=$serverHost, serverPort=$serverPort\nCould not connect to server socket: $!\n" if $ct++ > 3;
    sleep 3;
  }
}
print $hostSock "$jobid $host $localPort\n";
close($hostSock);

##open the bash shell handles
my ($read,$write);
my $opid = open2($read,$write,"bash -s");
die "ERROR: $@\n" if $opid =~ /^open2/;

##want this thing to die if distribjob doesn't contact it within 20 minutes.
$SIG{'ALRM'} = sub { die "alarm: no connections so timing out\n"; };
alarm(1200);

##want to exit after two connections as no more are required by distribjob and
##need to exit if distribjob exits prematurely
$ct = 0;
while(my $ns = $sock->accept()){
  $ct++;
  alarm(0);
  while( my $input = <$ns>){
#    print STDERR "$ct: Received: $input";
    if($input =~ /^closeAndExit/){
      print $ns "Goodbye\n0.$endMatchString\n";
      close $ns;
      &cleanUp;
      exit(0);
    }
    if($input =~ /subtaskInvoker\s+(\S+)/){
      chomp $input;
      $input .= " > $1/subtask.output 2> $1/subtask.stderr\n";
      print $write "$input"."echo \$?.$endMatchString\n";
    }else{
      print $write "$input"."echo \$?.$endMatchString\n";
    }
    while(<$read>){
      print $ns $_;
      last if /^(\d+)\.$endMatchString/;
    }
  }
  close $ns;
  last if $ct >= 2;
}
&cleanUp;

sub cleanUp {
  close $read;
  close $write;
  close($sock);
}
