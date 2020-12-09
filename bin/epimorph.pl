#!/usr/bin/env perl
use Getopt::Long;
use vars qw($rootdir);
BEGIN{
  use File::Basename;
  use Cwd 'abs_path';
  $rootdir = dirname(dirname(abs_path($0)));
};
use lib "$rootdir/lib";
use MSA;
use strict;
use warnings;

# Author:  Jacob Glanville 
# Contact: jake@distributedbio.com

############################### Arguments ###############################

my($epitopedb,$alignment,$numbering,$matchrules,$sanity_check,$centivax_min_conc,$forbidden_contact_list)=GatherOptions();
 
############################### Inputs    ###############################

# numbering
my %msaID2structID=msaID2structID($numbering);		# import the map of of the structural positions to the
							# multiple sequence alignment of all reference sequences
							# msaID    = column in alignment to structural ID.
							# structID = A:15 chain and position in structure
							# Note that each msa column can have multiple
							# structIDs (chain symmetry)
my %structID2msaID=processStructID2msaID($numbering);   # structIDs  have a one-to-one relationship to msaIDs
 
my %msaID2RefAA=processAARefNumberingFile($numbering);  # this stores the {refseq}{MSAID}=AA relationship
							# found in numbering
 
# include an optional forbidden contact list - for structurally forbidden positions
# ignore any epitope that contains these
my %forbidden_structIDs=processForbiddenContactList($forbidden_contact_list);

# upload the alignment of structure and reference sequences
my $msa=MSA->new();
   $msa->loadMSA($alignment);
my $seqs   = $msa->getSeqCount();
my $alnlen = $msa->getMSALength();

# matchrules - find the immunogen matching sequence or sequences in the alignment
# note that the final matching item will become the immunogen
# make very sure that this doesn't cause a problem for centivax tracing  
my ($match_id,$matchrules_nmatches)=getMatchSeq(\$msa,$matchrules);


# a lookup table that will let us go from a sequence name in the MSA and a structural position
# on the structural reference to a specific amino acid   
my %structID2AA=getStructID2AAHash(\$msa,\%msaID2structID); # msaHash{seqName}{structID} = AA 

# upload the epitope database
open(EPITOPE,$epitopedb);
my @elines=<EPITOPE>;
chomp(@elines);

# ok now check the vaccination against all sequences
my %epitope_breadth_counts=(); 

# iterating through every single epitope in the database
my %total_valid_centivax_epitopes=();
for(my $e=0;$e<scalar(@elines);$e++){
  #print "epitope $e\n";
  my $epitope_spreading="";

  # make an array of positions in the epitope
  my @positions=split(/ /,$elines[$e]);

  # get the epitope for the match seq (the immunogen in the case of a single immunogen)
  # it's just one of the sequences otherwise 
  my $match_epitope=getEpitope(\$msa,\@positions,\%structID2AA,$match_id);   

  # run an integrity check to make sure that the reference sequence is indeed
  # what we expected it to be 
  if($sanity_check){
    my $ref_epitope="";
    for(my $p=0;$p<scalar(@positions);$p++){
      $ref_epitope.=$msaID2RefAA{$structID2msaID{$positions[$p]}}
    }
    unless($match_epitope eq $ref_epitope){
      print "e=$e SANITY CHECK FAIL $match_epitope $ref_epitope\n";
      exit;
    }
  }

  # now get the centivax epitopes
  # this has will contain all the unique amino acid sequences at this epitope
  # and the number of times that particular sequence was observed, an important
  # measurement for determining effective concentration of each sequence
  my %centivax_epitopes=();
  if($centivax_min_conc>1){                                           # if this is a centivax run
    for(my $s=0;$s<$seqs;$s++){                                       # go through each sequence
      if($msa->getHeader($s) =~ m/centivax/i){                        # if this is a centivax sequence
        if($msa->getHeader($s) =~ m/$matchrules/){                    # if the header matches the matchrule (H1N1, etc)
          my $epitope=getEpitope(\$msa,\@positions,\%structID2AA,$s); # get the epitope sequence for this sequence
          add2Hash(\%centivax_epitopes,$epitope);                     # and add it to the epitopes for this position
        }
      }
    }
  }

  # now determine if this contains at least one valid centivax epitope
  # it is valid if at least n members of centivax share an epitope here
  # and thus there is sufficient concentration for this site to contribute
  # to the pool
  # note that only the pairs that are above 1 are valid centivax epitopes
  my $valid_centivax_epitope=1;
  if($centivax_min_conc>1){
    $valid_centivax_epitope=0; 
    my @centivax_epitope_list=keys %centivax_epitopes;
    for(my $c=0;$c<scalar(@centivax_epitope_list);$c++){
      if($centivax_epitopes{$centivax_epitope_list[$c]}>=$centivax_min_conc){
        $valid_centivax_epitope=1;
        $total_valid_centivax_epitopes{$centivax_epitope_list[$c]}=1;
        #$c=100000000;
      }     
    }
  }


  # now determine if this is a valid structurally permitted epitope
  # it is invalid if it contains a residue from the $forbidden_contact_list
  my $valid_contact_epitope=1;
  if(-f $forbidden_contact_list){
    for(my $p=0;$p<scalar(@positions);$p++){
      if(defined($forbidden_structIDs{$positions[$p]})){
        $valid_contact_epitope=0;
      }
    }
  }

  # now get the epitope for each member of the alignment

  # handle this differently if you are processing a centivax concentration ensembl
  if($centivax_min_conc>1){					# this is the centivax case where the position is only processed if
    if($valid_centivax_epitope){				# the $valid_centivax_epitope==1
      for(my $s=0;$s<$seqs;$s++){				# and in that case an epitope is a match if it hits any centivax epitope
        my $epitope=getEpitope(\$msa,\@positions,\%structID2AA,$s);							
	if(defined($centivax_epitopes{$epitope})){		# that is >=$centivax_min_conc
          if($centivax_epitopes{$epitope}>=$centivax_min_conc){
            add2Hash(\%epitope_breadth_counts,$msa->getHeader($s));
          }
        }
      }
    }
  }else{ 							# the normal case where there isn't a centivax ensemble to deal with
    if($valid_contact_epitope){       				# in this case just mark all the epitopes that match
      my %these_epitopes=();
      for(my $s=0;$s<$seqs;$s++){   
        my $epitope=getEpitope(\$msa,\@positions,\%structID2AA,$s);
        if($epitope eq $match_epitope){
          $epitope_spreading .= " " . $msa->getHeader($s);
          add2Hash(\%epitope_breadth_counts,$msa->getHeader($s));
          #print $epitope . "\t" . $msa->getHeader($s) . "\n";
        }
      }
      #print $elines[$e] . $epitope_spreading . "\n";
    }
  }
}

# mark any never found as 0 and get max
#my $max=0;
for(my $s=0;$s<$seqs;$s++){
  if(!defined($epitope_breadth_counts{$msa->getHeader($s)})){
    $epitope_breadth_counts{$msa->getHeader($s)}=0;
  }
  #else{
  #  if($epitope_breadth_counts{$msa->getHeader($s)}>$max){
  #    $max=$epitope_breadth_counts{$msa->getHeader($s)};
  #  }
  #}
}

# ok, now print the results

my @sorted_headers = sort keys %epitope_breadth_counts;

for(my $s=0;$s<scalar(@sorted_headers);$s++){
  print   $epitope_breadth_counts{$sorted_headers[$s]} . "\t"
        . ($epitope_breadth_counts{$sorted_headers[$s]}/scalar(keys %total_valid_centivax_epitopes)) . "\t"
        . $sorted_headers[$s] . "\n";
}
exit;

############################### Subroutines #############################

sub processForbiddenContactList {
  my($forbidden_contact_list)=@_;
  my %forbidden_structIDs=();
 
  if(-f $forbidden_contact_list){
    # final-forbidden-contacts-for-stem-binders.txt
    open(LIST,$forbidden_contact_list);
    my @forbidden=<LIST>;
    chomp(@forbidden);
    close(LIST);
    for(my $x=0;$x<scalar(@forbidden);$x++){
      $forbidden_structIDs{$forbidden[$x]}=1;
    }
  }
  return %forbidden_structIDs;
}

sub add2Hash {
  my($this_hash,$this_item)=@_;
  if(defined($$this_hash{$this_item})){
    $$this_hash{$this_item}+=1;
  }else{
    $$this_hash{$this_item}=1;
  }
}

sub getEpitope {
  my($msa,$epitope_positions,$structID2AA,$s)=@_;
  
  my $epitope="";

  my $header=$$msa->getHeader($s);
  for(my $p=0;$p<scalar(@$epitope_positions);$p++){
      $epitope.=$structID2AA{$header}{$$epitope_positions[$p]};
  }
  return $epitope;
}

sub getStructID2AAHash {
  my($msa,$msaID2structID)=@_;

  my %structID2AA=();

  my $seqs   = $$msa->getSeqCount();
  my $alnlen = $$msa->getMSALength();

  for(my $s=0;$s<$seqs;$s++){
    my $header=$$msa->getHeader($s);
    for(my $p=0;$p<$alnlen;$p++){
      my @structIDs=split(/,/,$$msaID2structID{$p});
      my $residue=$$msa->getChar($s,$p);
      #print "$s $p $residue $header " . $structIDs[0] . "\n";
      for(my $i=0;$i<scalar(@structIDs);$i++){
        $structID2AA{$header}{$structIDs[$i]}=$residue;
      }
    }
  }  
  return %structID2AA;
}

# ($match_id,$matchrules_nmatches)=getMatchSeq
 
sub getMatchSeq {
  my($msa,$matchrules)=@_;
  my $match_id="";
  my $matchrules_nmatches=0;
  my $seqs   = $$msa->getSeqCount();

  for(my $s=0;$s<$seqs;$s++){
    if($$msa->getHeader($s)=~m/$matchrules/){
      $match_id=$s;
      $matchrules_nmatches+=1;
      print $$msa->getHeader($s) . "\n";
    }
  }
  return ($match_id,$matchrules_nmatches);
}

# numbering

sub processStructID2msaID {
  my($numbering)=@_;

  my %structID2msaID=();
  open(NUMBER,$numbering);
  my @lines=<NUMBER>;
  chomp(@lines);
  for(my $n=0;$n<scalar(@lines);$n++){
    my($msa_number,$residue_aa,$structure_number)=split(/\t/,$lines[$n]);
    $structID2msaID{$structure_number}=$msa_number;
  }
  return %structID2msaID; 
}

# msaID    = column in alignment to structural ID. Note that each msa column can have multiple structIDs (chain symmetry)
# structID = A:15 chain and position in structure  (this is a one-to-one relationship to msaIDs)

sub msaID2structID {
  my($numbering)=@_;

  my %msaID2structID=();
  open(NUMBER,$numbering);
  my @lines=<NUMBER>;
  chomp(@lines);
  for(my $n=0;$n<scalar(@lines);$n++){
    my($msa_number,$residue_aa,$structure_number)=split(/\t/,$lines[$n]);
    if(defined($msaID2structID{$msa_number})){
      $msaID2structID{$msa_number}.="," . $structure_number;
    }else{
      $msaID2structID{$msa_number} =$structure_number;
    }
  }
  return %msaID2structID;
}

sub processAARefNumberingFile {
  my($numbering)=@_;

  my %msaID2AA=();
  open(NUMBER,$numbering);
  my @lines=<NUMBER>;
  chomp(@lines);
  for(my $n=0;$n<scalar(@lines);$n++){
    my($msa_number,$residue_aa,$structure_number)=split(/\t/,$lines[$n]);
    #if(defined($msaID2AA{$msa_number})){
    #  $msaID2AA{$msa_number}.="," . $residue_aa;
    #}else{
      $msaID2AA{$msa_number} =$residue_aa;
    #}
  }
  return %msaID2AA;
}

sub GatherOptions {
  my $epitopedb     ="";
  my $alignment     ="";
  my $numbering     ="";
  my $antigen       ="";
  my $sanity_check  = 0;
  my $centivax_min_conc = 1;
  my $forbidden_contact_list="";


  GetOptions(
     "--epitopedb=s" => \$epitopedb,
     "--alignment=s" => \$alignment,
     "--numbering=s" => \$numbering,
     "--antigen=s"   => \$antigen,
     "--sanity_check=s" => \$sanity_check,
     "--centivax_min_conc=s" => \$centivax_min_conc,
     "--forbidden_contact_list=s" => \$forbidden_contact_list
  );

  unless( (-f $epitopedb) and (-f $alignment) and (-f $numbering)){
    print "\nUsage: $0\n";
    print "  Standard options\n";
    print "  --epitopedb=epitopes.txt  database of epitopes, with each line containing\n";
    print "                            space separated list of positions that make up \n";
    print "                            the epitope.\n";
    print "                            example: A:119 A:122 A:124 A:125 A:127 A:133 A:146 A:149 A:152 A:67\n";
    print "                                     B:125 B:147 B:150 B:151 B:160 B:161 C:4 C:5 D:11 D:9\n";
    print "                                     B:14 B:150 B:20 B:25 B:28 B:30 B:31 B:32 B:33 B:34\n";
    print "  --alignment=MSA.fa        a set of sequences aligned to the structural   \n";
    print "                            reference antigen. It is an MSA with every residue of\n";
    print "                            the reference structure present and no gaps relative to\n";
    print "                            the reference structure, such that MSA length = ref sequence length\n";
    print "                            example: >structure\n";
    print "                                     EDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKKHNGKLCDLDGVKPLILRDCSVAG\n";
    print "                                     >H1N1-2009\n";
    print "                                     EDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKKHNGKLCDLDGVKPLILRDCSVAG\n";
    print "  --numbering=number.txt    a mapping file assigning a structural number ID\n";
    print "                            to each position in the sequence of the structural\n";
    print "                            reference found in the MSA.fa. Format is three tabs\n";
    print "                            separated seqposition|AA|structuralnumber\n";
    print " 		               example: Residue_Number	Residue_AA	Structure_Number\n";
    print " 					0	E	A:4\n";
    print " 					1	D	A:5\n";
    print " 					2	Q	A:6\n";
    print " 					3	I	A:7\n";
    print "  --antigen=[state]         one of the sequence names as the antigen. then for each epitope\n";
    print "                            code will report all MSA variants that share that epitope with\n";
    print "                            antigen\n";
    print "  --sanity_check=1          [OPTIONAL] if antigen is set to reference structure sequence, it makes sure that\n";
    print "                            the numbering in numbering.txt always matches the alignment or returns\n";
    print "                            an error message. Note --antigen needs to be the structural ref for this!\n";
    print "  --centivax_min_conc=1     [OPTIONAL] minimum number of centivax members to share an epitope here for this epitope\n";
    print "                            to count as above threshold dose. If set to 3, then at least one epitope must be shared\n";
    print "                            between three members of centivax for this site to be counted as present in the mixture\n";
    print "  --forbidden_contact_list  a file that contains forbidden structIDs, one per line. For disallowing certain epitopes\n";
    print "                            (stem binding, etc)\n";
    exit;
  }

  return($epitopedb,$alignment,$numbering,$antigen,$sanity_check,$centivax_min_conc,$forbidden_contact_list);
}
