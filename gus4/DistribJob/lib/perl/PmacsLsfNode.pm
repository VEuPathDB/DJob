package DJob::DistribJob::PmacsLsfNode;

use DJob::DistribJob::LsfNode;
use strict;

our @ISA = qw(DJob::DistribJob::LsfNode);

sub getQueue {
  my ($self) = @_;
  my $q = 'normal';
  $q = 'max_mem30' if ($self->{memPerNode} > 3);
  $q = 'max_mem64' if ($self->{memPerNode} > 30);
  return $q;
}

1;
