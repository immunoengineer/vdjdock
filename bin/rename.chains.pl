#!/usr/bin/perl -w

my $file   = $ARGV[0];
my $chain1 = $ARGV[1];
my $chain2 = $ARGV[2];

open(FILE,$file);
my @lines=<FILE>;
close(FILE);
chomp(@lines);

my %replacement = ($chain1,"X",$chain2,"Y");

#ATOM   5096  N   GLY H   2      12.740 -17.748  42.513 19     1 1.63         -0.15
#ATOM   5097  CA  GLY H   2      12.212 -16.608  43.240 19     1 1.99          0.10
#ATOM   5098  C   GLY H   2      12.688 -15.299  42.640 19     1 1.67          0.60
#ATOM   5099  O   GLY H   2      11.873 -14.464  42.237 19     1 1.38         -0.55
#ATOM   5100  HA2 GLY H   2      11.243 -16.627  43.216 19     0 1.90          0.00
#ATOM   5101  HA3 GLY H   2      12.501 -16.649  44.165 19     0 1.90          0.00

for(my $x=0;$x<scalar(@lines);$x++){
  my $first=substr($lines[$x],0,21);
  my $second=substr($lines[$x],21,1);
  my $third=substr($lines[$x],22);
  #print $lines[$x] . "\n";
  print $first . $replacement{$second} . $third . "\n";
  #print "\n";
}
