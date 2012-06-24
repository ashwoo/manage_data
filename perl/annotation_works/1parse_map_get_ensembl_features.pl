##CONSIDERS MAPPING SEQUENCE TO GENOME TAKING CARE OF POSSIBLE INTRONIC REGIONS. HENCE SUITABLE TO BACKMAP FULL LENGTH TRANSCRIPTS. A FILE CONTAINING BLAST RESULTS OF SEQUENCES AGAINST GENOME IS REQUIRED. THE FILE IS CONVERTED INTO BED FORMAT WITH 12 COLUMNS. FOR EACH MAPPED FEATURE PROXIMAL PROTEIN CODING GENES ARE EXTRACTED FROM THE ENSEMBL VIA API MODULES. A FEATURE FILE IS GENERATED WHICH CAN BE UPLOADED INTO SQL.

#BLAST COMMANDS SUITABLE FOR MAPPING SEQUENCES TO GENOME WITHIN SAME SPECIES AND BETWEEN DIFFERENT SPECIES
#blastn -query query_filtered.fa -db /home/swaraj/db/genomes/mm9/ensembl63/mm9 -out query_filtered.backmapped.mm9 -num_threads 22 -culling_limit 1 : MAPPING PARAMETERS FOR BLAST
#blastn -query zebrafish.fasta -db /home/swaraj/db/genomes/mm9/ensembl63/mm9 -out subhit_filtered.backmapped.mm9 -num_threads 22 -culling_limit 1 -task blastn : FOR INTER SPECIES

#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::Search::Tiling::MapTiling;
use Bio::SearchIO;
use Bio::EnsEMBL::Registry;

use FindBin;
use lib $FindBin::Bin;
use flanks;
use overlaps;
#GLOBAL PARAMETERS 
my $org = "Danio rerio";
my $ftype = "protein_coding";
my $filename = 'filtered.backmapped.zv9'||die "No such file or directory\n";
my $percent = 90;#KEEP 90 FOR SAME SPECIES AND 70 FOR DIFFFERENT SPECIES
open(FH, ">>backmapped.features");
open(FV, ">>backmapped.bed");

print FH "ID\tRegion\tStart\tEnd\tStrand\tRight Flank\tDesc\tLeft Flank\tDesc\tGeneOverlap\tDesc\tGeneType\tTranscriptOverlap\tTranscriptType\tExonOverlap\n";

#CONNECT TO THE ENSEMBL DATABASE                               
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db (-host => 'ensembldb.ensembl.org',-user => 'anonymous');

#GET ALL THE CORE ADAPTORS
my $gene_adaptor = $registry->get_adaptor ("$org", "Core", "Gene");
my $transcript_adaptor = $registry->get_adaptor ("$org", "Core", "Transcript");
my $exon_adaptor = $registry->get_adaptor ("$org", "Core", "Exon");
my $slice_adaptor = $registry->get_adaptor ("$org", "Core", "Slice");
#GET ALL THE OTHER FEATURE ADAPTOR
my $dafa = $registry->get_adaptor( "$org", "OtherFeatures", "DnaAlignFeature" );
my $align_slice_adaptor = $registry->get_adaptor( "$org", "OtherFeatures", "Slice" );  

#PARSE THE BLASTOUTPUT
my $searchio = Bio::SearchIO->new( -format => 'blast', -file   => $filename );
                                   
while ( my $result = $searchio->next_result() ) {
  my $query = $result->query_name;
  my $description = $result->query_description;

  while( my $hit = $result->next_hit()) {
    #EXTRACT THE COORDINATE NAME AND SEQ REGION NAME
    my @DATA = split("\:| ", $hit->description());
    my $csn = $DATA[1];
    my $csn1 = $csn;
    $csn1 =~s/chromosome/chr/;
    my $srn = $hit->name();
    my $region = $csn1.$srn;
    my @blocks;
    my @starts;
    my $tiling = Bio::Search::Tiling::MapTiling->new($hit); 
    my $qprevious;
    my $qprevious_minus;
    my %check;
    #COLLECT STRAND INFORMATION FOR BEST HIT TO BUILD A CONTEXT FOR THE TILING
    my @chkhsp = $hit->hsps();
    my $strand = $chkhsp[0]->strand('hit');
    my $context = $tiling->_context( -type => 'subject', -strand=> $strand);
    #print $tiling->coverage_map_as_text('subject',$context,'LEGEND'); #USE TO VISUALIZE RESULTS
    #TILE HSPS FOR THE BEST HIT
    my @hsps = $tiling->next_tiling('subject', $context);
    my @hstarts; 
    my @hends;

    foreach my $hsp(@hsps){  
      my $percentID = $hsp->percent_identity();

      next if ($percentID < $percent);
      next if exists $check{ $hsp->start('query') };
      #CHECK FOR A HSP TO BE AGAIN REPEATED IN SOME OTHER PART OF THE GENOME
      if($strand == -1){ $check{ $hsp->start('query')} = $hsp->end('hit'); }
      if($strand == +1){ $check{ $hsp->start('query')} = $hsp->start('hit'); }
      
      my $start = $hsp->start('hit')-1; #IN A BED FILE THE START POSTION IS ONE NUCLEOTIDE BEHIND THE ACTUAL START  AND CONSIDERED AS 0   
      my $end = $hsp->end('hit');
      push(@hstarts, $start);
      push(@hends, $end);
      #NULLIFY THE EFFECT OF HSP OVERLAP ON THE BED FILE BLOCK DESCRIPTION
      if($qprevious && $strand == +1){
        if($qprevious > $hsp->start('query')){ 
          my $difference = $qprevious - $hsp->start('query') + 1;
          $start = $start + $difference;
        }
        if($qprevious == $hsp->start('query')){ 
          my $difference = 1;
          $start = $start + $difference;
        }
      }
      #NULLIFY THE EFFECT OF HSP OVERLAP ON THE BED FILE BLOCK DESCRIPTION IF THE ALIGNMENT IS WITH THE MINUS STRAND
      if($qprevious_minus && $strand == -1){
        if($qprevious_minus < $hsp->end('query')){ 
          my $difference = $hsp->end('query') - $qprevious_minus + 1;
          $start = $start + $difference;
        }
        if($qprevious_minus == $hsp->end('query')){ 
          my $difference = 1;
          $start = $start + $difference;
        }
      }

      my $bedstart = $start-$hstarts[0];  
      my $bedlen = $end - $start;
      #print "$start\t$end\t$bedlen\n";
      push(@blocks, $bedlen);
      push(@starts, $bedstart); 
      $qprevious = $hsp->end('query');
      $qprevious_minus = $hsp->start('query');
    }
    my $hstart = $hstarts[0];
    my $hend = $hends[-1];
    my $block = join(',', @blocks);
    my $start = join(',', @starts);
    my $count = @blocks;
    my $symbol;
    $symbol = '-' if $strand == -1;
    $symbol = '+' if $strand == +1;
    print FV "$region\t$hstart\t$hend\t$query:$description\t0\t$symbol\t$hstart\t$hend\t0,0,0\t$count\t$block\t$start\n";

    #GET FLANKING PROTEINS FOR EACH FEATURE
    my $flanks = flanks->flanking($org,$csn,$srn,$hstart,$hend,$ftype);
    my $rflank = $flanks->{right_flank};
    my $lflank = $flanks->{left_flank};
    #BUILD FLANK OBJECTS AND GET DESCRIPTIONS
    my $rflank_desc = '\N';
    if($rflank ne '\N'){ $rflank_desc = $gene_adaptor->fetch_by_stable_id($rflank)->description; }
    my $lflank_desc = '\N';
    if($lflank ne '\N'){ $lflank_desc = $gene_adaptor->fetch_by_stable_id($lflank)->description; }      

    #FIND ALL GENE OVERLAPS FOR EACH FEATURE
    my $list_overlaps = overlaps->overlapping($org,$csn,$srn,$hstart,$hend);
    while (my $overlap = shift(@$list_overlaps)){

      my $ogene = $overlap->{gene};
      my $ogene_desc = '\N';
      
      if($ogene ne 'Intergenic'){ $ogene_desc = $gene_adaptor->fetch_by_stable_id($ogene)->description; }
      my $ogene_type = $overlap->{gene_type};
      my $otranscript = $overlap->{transcript};
      my $otranscript_type = $overlap->{transcript_type};
      my $oexon = $overlap->{exon};
        
      #my $estart = $exon_adaptor->fetch_by_stable_id($oexon)->seq_region_start;
      #my $eend = $exon_adaptor->fetch_by_stable_id($oexon)->seq_region_end;
        
      print FH $query."\t".$csn.$srn."\t".$hstart."\t".$hend."\t".$strand."\t".$rflank."\t".$rflank_desc."\t".$lflank."\t".
               $lflank_desc."\t".$ogene."\t".$ogene_desc."\t".$ogene_type."\t".$otranscript."\t".$otranscript_type."\t".
               $oexon."\n";
    }


    last
  }
}

  
