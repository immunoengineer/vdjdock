#!/usr/bin/perl -w

unless(-f $ARGV[0]){
  print "correct.pl file.pdb 		removes offending non-PDB formatting from ZDOCK\n";
  exit;
}

open(FILE,$ARGV[0]);
my @lines=<FILE>;
chomp(@lines);
close(FILE);

for(my $x=0;$x<scalar(@lines);$x++){
  if($lines[$x]=~m/^ATOM/){
    my $substring=substr($lines[$x],0,54);
    print $substring . "\n";
  }else{
    print $lines[$x] . "\n";
  }
}
