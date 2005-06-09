#!@perl@

use strict;
use IO::Socket;
use IPC::Open2;

$| = 1;

my $endMatchString = 'FooCmdEnd';

my $host = `hostname`;
chomp $host;

print "$host\n";

my $sock = new IO::Socket::INET (
                                 LocalHost => $host,
                                 LocalPort => '7070',
                                 Proto => 'tcp',
                                 Listen => 5,
                                 Reuse => 1,
                                );
die "Could not create socket: $!\n" unless $sock;

##open the bash shell handles
my ($read,$write);
my $opid = open2($read,$write,"bash -s");
die "ERROR: $@\n" if $opid =~ /^open2/;

##want this thing to die if distribjob doesn't contact it within 3 minutes.
$SIG{'ALRM'} = sub { die "alarm: no connections so timing out\n"; };
alarm(180);

##want to exit after three connections as no more are required by distribjob and
##need to exit if distribjob exits prematurely
my $ct = 0;
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
    print $write "$input"."echo \$?.$endMatchString\n";
    while(<$read>){
      print $ns $_;
      last if /^(\d+)\.$endMatchString/;
    }
  }
  close $ns;
  last if $ct >= 3;
}
&cleanUp;

sub cleanUp {
  close $read;
  close $write;
  close($sock);
}
