#!/usr/bin/perl
use strict;
use warnings;

=head1 NAME

neighbors: A module to identify overlapping and nearest neighboring features of a gene

=head1 SYNOPSIS

$Neighbors = $PATH::neighbors->fetch(Organism,
                                     Coord_system_name,
                                     Sequence_region_name,
                                     start,
                                     end,
                                     Overlap_biotype); 
                                     


=head1 DESCRIPTION

Organism = Organism to which gene belongs;

Coord_system_name = Chromosome, Scaffold, Contig

Sequence_region_name = 1,2,14 etc

start = Top level coordinate of the feature start

end = Top level coordinate of the feature end

Overlap_biotype = Which type of genes flanking the query would you like to see eg. protein_coding, miRNA, lincRNA

=head1 RETURN VALUE
Returns an array of hashes with each hash containing the values for the following keys
overlapping_gene = "Gives the EnsEMBL stable id of the overlapping gene"
overlapping_gene type = "Gives the type of the overlapping gene"
overlapping_transcript = "Gives the EnsEMBL stable id of the overlapping transcripts"
overlapping_transcript type = "Gives the type of the overlapping transcripts"
overlapping_exon = "Gives the EnsEMBL stable id of the overlapping exon of each 
                    transcript"


=head1 AUTHORS

Remo Sanges & Swaraj Basu

=cut



package overlaps;
my $registry = 'Bio::EnsEMBL::Registry';


sub overlapping {

#STORE AS VARIABLES THE FEATURES OF THE GENE 
my $self = shift;
my $organism = shift;
my $csn = shift;
my $srn = shift;
my $start = shift;
my $end = shift;

#LOAD THE GENE AND SLICE ADAPTORS
my $gene_adaptor = $registry->get_adaptor($organism,"Core","Gene");
my $slice_adaptor = $registry->get_adaptor($organism,"Core","Slice");

#FETCH THE GENE AND ITS EXON FEATURES WHICH OVERLAP THE REGION WHERE THE miRNA GENE
#FALLS   
my $slice = $slice_adaptor->fetch_by_region($csn,$srn,$start,$end);
my $gene_overlaps = $slice->get_all_Genes;
my $count = 1;

#CREATE AN EMPTY ARRAY TO RETURN THE RESULTS AS AN ARRAY OF HASHES
my @return;
  
  #ITERATE FOR EACH GENE OVERLAPPING THE SLICE               
  while (my $goverlap = shift(@{$gene_overlaps})) {
    $count = 2;
    my $gtype = $goverlap->biotype;
    my $transcripts = $goverlap->get_all_Transcripts;
    
    while (my $transcript = shift (@{$transcripts})) {
      next unless overlaps ($transcript, $start, $end);
      my $ttype = $transcript->biotype;
      #print $gene->stable_id."\t".$transcript->stable_id."\n";
      
      my $eoverlap = 'intron';
      my $exons = $transcript->get_all_Exons;
      
      while (my $exon = shift(@{$exons})) {
        $eoverlap = $exon->stable_id if overlaps($exon, $start, $end);
      }
    
      #print  $gene_id."\t".$srn."\t".$csn."\t".$start."\t".$end."\t".
      #       $goverlap->stable_id."\t".$goverlap->biotype."\t".$overlap."\t";
    
      #print  $flanking_right."\t".$flanking_left."\n";
      push @return, {gene => $goverlap->stable_id,
                     gene_type => $gtype,
                     transcript => $transcript->stable_id,
                     transcript_type => $ttype,
                     exon => $eoverlap};
    }
  }
  
  if($count < 2) {
   
    #print  $gene_id."\t".$srn."\t".$csn."\t".$start."\t".$end."\t".
    #      "NA\tNA\tIntergenic\t";
    
    #print  $flanking_right."\t".$flanking_left."\n";
    push @return, {gene => 'Intergenic',
                   gene_type => '\N',
                   transcript => '\N',
                   transcript_type => '\N',
                   exon => '\N'};
    return \@return;   
  }
  else{
    
    return \@return;
  }
}
###########################################################################################
###########################################################################################
###########################################################################################

#FIND WHETHER TWO ENTITIES ARE OVERLAPPING OR NOT
sub overlaps {
  my $f1 = shift;
  my $f2start = shift;
  my $f2end = shift;

  return($f1->seq_region_end >= $f2start and 
         $f1->seq_region_start <= $f2end);
}


