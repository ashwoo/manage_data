#THIS SCRIPT USES A 1KB WINDOW ON EITHER SIDE OF THE HIT COORDINATE TO FIND ESTS USING LOCAL EST DATABASE
#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

use Bio::SearchIO;
use Bio::EnsEMBL::Registry;

use FindBin;
use lib $FindBin::Bin;
#use flanks;
use overlaps;

#USE THE SQL DATA MANAGEMENT MODULES
use DBI;
use DBD::mysql;

#MySQL CONFIG VARIABLES
my $platform  = "mysql";
my $host      = "localhost";
my $port      = "3306";
my $user      = "mysql_dev";
my $pw        = "dEvEl0pEr";
my $database  = "danRer7";

#DATA SOURCE NAME
my $dsn = "dbi:$platform:$database:$host:$port";
#PERL DBI CONNECT AND CREATE TABLES
my $dbh =  DBI->connect($dsn, $user, $pw, { RaiseError => 1, AutoCommit => 0 } )  || die "Could not connect to database: $DBI::errstr";
#GLOBAL PARAMETERS 
my $org = "Danio rerio";
my $ftype = "protein_coding";

my $window = 1000;#A WINDOW TO EXTEND ON EITHER SIDE OF THE QUERY COORDINATE IN SEARCH OF EST

#print "ID\tname\tregion\tqstart\tqend\tfID\tfstart\tfend\tfstatus\n";

#CONNECT TO THE ENSEMBL DATABASE                               
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db (-host => 'ensembldb.ensembl.org',-user => 'anonymous');
my $dafa = $registry->get_adaptor( "$org", "OtherFeatures", "DnaAlignFeature" );
my $slice_adaptor = $registry->get_adaptor( "$org", "OtherFeatures", "Slice" );  
my $gene_adaptor = $registry->get_adaptor ("$org", "Core", "Gene");

open(FH, "<query.txt");

while(<FH>){  
                                 
  my $line = $_;
  chomp($line);
  my @DATA = split("\t", $line);
  my $id = $DATA[0];
  my $name = $DATA[1];
  my $srn1 = $DATA[2];
  my $start = $DATA[3];
  my $end = $DATA[4];
  
  my $range_start = $start - $window;
  my $range_end = $start + $window;
  my $csn;
  my $srn;
  if($srn1 =~m/chromosome/){
    $srn1 =~s/chromosome/chr/;
    $csn = 'chromosome';
    $srn = $srn1;
    $srn =~s/chr//;
  }
  else{
    $csn = 'scaffold';
    $srn = $srn1;
  }

  #print "$query\t$csn\t$srn\t$range_start\t$range_end\t";
  my $sth = $dbh->selectall_arrayref("SELECT * FROM all_est where tEnd >= '$range_start' AND tStart <= '$range_end' AND tName = '$srn1';", {Slice => {} });
  if(exists $sth->[0]->{qName}){
    #print $id."\n\n";
    foreach my $row(@$sth){
      my $fname = $row->{qName};
      my $chr = $row->{tName};
      my $ori = $row->{tStart};
      my $cul = $row->{tEnd};
      my @blocks = split("\,",$row->{blockSizes});
      my @fstarts = split("\,", $row->{tStarts});
      my $count = @blocks;
      my $fstatus = 'neighboring';
      my $foverlap = 'protein_coding';
      my $fsubfeature = 'exonic';
      foreach (my $i = 1; $i <= $count; $i++){
        my $fstart = $fstarts[$i] + 1;
        my $fend = $fstarts[$i] + $blocks[$i];
        next if $fstatus eq 'overlapping';
        if($end >= $fstart && $start <= $fend){ $fstatus = 'overlapping';}  
        #FIND ALL GENE OVERLAPS FOR EACH FEATURE
        #my $list_overlaps = overlaps->overlapping($org,$csn,$srn,$fstart,$fend);
        #while (my $overlap = shift(@$list_overlaps)){
         # my $ogene = $overlap->{gene};
         # my $ogene_type = $overlap->{gene_type};
         # my $otranscript = $overlap->{transcript};
         # my $otranscript_type = $overlap->{transcript_type};
         # my $oexon = $overlap->{exon};
         # next if $ogene_type eq 'protein_coding';
          
         # if($ogene_type eq 'protein_coding'){
         #   if($oexon eq '\N'){ $fsubfeature = "intronic"; }
         # }
         # elsif($ogene_type eq '\N'){
         #   $foverlap = "intergenic";
          #  $fsubfeature = '\N';
         # }
         # else{
          #  $foverlap = 'non_protein_coding';
          #  if($oexon eq '\N'){ $fsubfeature = "intronic"; }
          #}
        #}
      }
      print $id."\t".$name."\t".$srn1."\t".$start."\t".$end."\t".$fname."\t".$ori."\t".$cul."\t".$fstatus."\n";   
    }
  }
  else{
    print $id."\t".$name."\t".$srn1."\t".$start."\t".$end."\t\\N\t\\N\t\\N\t\\N\n";
  }
}

#PERL DBI DATA DISCONNECT
$dbh->disconnect();
  
