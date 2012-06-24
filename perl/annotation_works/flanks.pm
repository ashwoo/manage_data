#!/usr/bin/perl
use strict;
use warnings;

=head1 NAME

flanks: A module to identify nearest neighboring features of a gene

=head1 SYNOPSIS

$flanks = flanks->fetch(Organism,
                        Coord_system_name,
                        Sequence_region_name,
                        start,
                        end,
                        FlankBiotype); 
                      

=head1 DESCRIPTION

Organism = Organism to which gene belongs;

Coord_system_name = Chromosome, Scaffold, Contig

Sequence_region_name = 1,2,14 etc

start = Top level coordinate of the feature start

end = Top level coordinate of the feature end


=head1 RETURN VALUE
Returns an array of hashes with each hash containing the values for the following keys

right_flank = "Gives the EnsEMBL stable id of the right flanking gene"
left_flank = "Gives the EnsEMBL stable id of the left flanking gene"

=head1 AUTHORS

Remo Sanges & Swaraj Basu

=cut



package flanks;
my $registry = 'Bio::EnsEMBL::Registry';

#NEAREST NEIGHBORS OF GENE
sub flanking {
  my $self = shift;
  my $org = shift;
  my $csn = shift;
  my $srn = shift;
  my $gene_start = shift;
  my $gene_end = shift;
  my $biotype = shift;
  #print "$self\t$org\t$csn\t$srn\t$gene_start\t$gene_end\t$biotype\n";
  
  #PARAMETERS TO BE USED FOR BOTH FLANKS
  
  my $ga = $registry->get_adaptor ("$org", "Core", "Gene");
  my $sa = $registry->get_adaptor ("$org", "Core", "Slice");
  my $interval = 10000;
  my $lid;
  my $rid;
    
  
  #LOOP TO FIND THE RIGHT FLANKING GENE
  my $rcounter = 0;
  my $rstart = $gene_end + 1;
  my $rend;
  my $rlast = 0;
  my $seq_region_length = $sa->fetch_by_region($csn, $srn)->length;
  
  until ($rcounter > 0 || $rlast > 0) {
    $rid = '\N';
    #MOVE IN FORWARD DIRECTION USING THE WINDOW SIZE SPECIFIED BY $interval
    $rend = $rstart + $interval;
    
    #CHECK THAT THE WINDOW DOESNT EXCEED THE SEQUENCE REGION LENGTH 
    if ($rend > $seq_region_length) {
      $rend = $seq_region_length;
      $rlast = 1;
    }
        
    my $slice_right = $sa->fetch_by_region ($csn, $srn, $rstart,  
                                            $rend);
    my $nearest_genes = $slice_right->get_all_Genes_by_type($biotype);
    $rstart = $rend + 1;
    
    #MOVE AHEAD IF THE REGION DOES NOT HAVE A PROTEIN CODING GENE
    next unless exists $nearest_genes->[0];
        
    #COUNT THE NUMBER OF GENES FALLING WITHIN WINDOW; THE 1st POSITION IS NEAREST GENE
    $rcounter = @{$nearest_genes};
    
    #DONT WANT AN OVERLAPPING GENE
    if ($nearest_genes->[0]->seq_region_start < $gene_end) {
      $rstart = $nearest_genes->[0]->seq_region_end + 1;
      $rcounter = 0;
    }
    next if $rcounter == 0;
    
    #GET THE OVERLAPPING GENE ID
    $rid = $nearest_genes->[0]->stable_id;
    
  }
  
  
  #LOOP TO FIND THE LEFT FLANKING GENE
  my $lcounter = 0;
  my $lstart;
  my $lend = $gene_start - 1;
  my $llast = 0;
  
  until ($lcounter > 0 || $llast > 0) {
    $lid = '\N';
    #MOVE IN REVERSE DIRECTION USING THE WINDOW SIZE SPECIFIED BY $interval
    $lstart = $lend - $interval;
    
    #CHECK THAT THE WINDOW DOESNT EXCEED THE SEQUENCE REGION LENGTH
    ($lstart = 1 && $llast = 1 ) if $lstart <= 0;

    my $slice_left = $sa->fetch_by_region ($csn, $srn, $lstart,  
                                           $lend);
    my $nearest_genes = $slice_left->get_all_Genes_by_type($biotype);
    $lend = $lstart - 1;
    
    #MOVE AHEAD IF THE REGION DOES NOT HAVE A PROTEIN CODING GENE
    next unless exists $nearest_genes->[0];
    
    #COUNT THE NUMBER OF GENES FALLING WITHIN WINDOW AND POSITION OF LAST GENE
    $lcounter = @{$nearest_genes};
    my $array_position = $lcounter - 1;
    
   
    #DONT WANT AN OVERLAPPING GENE
    if ($nearest_genes->[$array_position]->seq_region_end > $gene_start) {
      $lend = $nearest_genes->[$array_position]->seq_region_start - 1;
      $lcounter = 0;
    }
    next if $lcounter == 0;
    
    #GET THE OVERLAPPING GENE ID
    $lid = $nearest_genes->[$array_position]->stable_id;
    
  }
  
  #RETURN THE FLANKING GENES AS Bio::EnsEMBL::Gene OBJECTS
  return { right_flank => $rid,
           left_flank  => $lid };
}


