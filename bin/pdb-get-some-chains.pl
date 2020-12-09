#!/usr/bin/perl -w

use strict;

my $file=$ARGV[0];
my $chainA=$ARGV[1];
my $chainB=$ARGV[2];

open(FILE,$file);
my @line=<FILE>;
chomp(@line);

for(my $x=0;$x<scalar(@line);$x++){
  #ATOM    283  CE  LYS A  39      17.796 -27.967 -43.320  1.00 90.52           C
  #ATOM    284  NZ  LYS A  39      16.514 -28.199 -42.571  1.00 90.12           N
  if($line[$x]=~m/^ATOM/){
    #print "found an ATOM line:" . $line[$x] . "\n";
    
    my($atom,$atomnumber,$chars,$aatype,$chain,$aanumber,$altx,$y,$z,$bla,$bla2,$atomtype)=split(/  */,$line[$x]);

    #print "chain=|$chain| and aanumber=|$aanumber|\n";
    if(($chain =~ m/^$chainA/) or ($chain =~ m/^$chainB/)){ # or ($chain =~ m/^$chainC/)){
      #if( ($aanumber < 40) or ($aanumber > 307)){
        print $line[$x] . "\n";
      #}
    }#else{
    #  print $line[$x] . "\n";
    #}
  }else{
    # dont print heteroatoms
    unless($line[$x]=~m/^HETATM/){
      print $line[$x] . "\n";
    }
  }
}
