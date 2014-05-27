#!/usr/bin/perl

use strict;
use IO::Socket;
use File::Basename;

my $localDir = shift;
my $jobid = shift;
my $serverHost = shift;
my $serverPort = shift;
my $status = shift;

# print STDERR "informing controller ($localDir, $jobid, $serverHost, $serverPort, $status)\n";

my $hostSock;
my $ct = 0;
until($hostSock){
  $hostSock = new IO::Socket::INET (
                                    PeerAddr => $serverHost,
                                    PeerPort => $serverPort,
                                    Proto => 'tcp',
                                   );
  unless($hostSock){
    if($ct++ > 5){
      die "$jobid: Could not create socket:\n$!\n";
    }
    sleep 3;
  }
}

my $slot = basename($localDir);

print $hostSock "$jobid $slot $status\n" if $hostSock;

close($hostSock) if $hostSock;
