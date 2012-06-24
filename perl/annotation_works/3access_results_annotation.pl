#!/usr/bin/perl
#THE SCRIPT COMPARES BETWEEN TWO DATASETS TO FILTER OUT RESULTS WHICH SHOW HOMOLOGY AND OVERLAP ESTS WHILE NOT OVERLAPPING ANY PROTEIN CODING GENE. THE TWO DATASETS CAN BE CHOSEN BY CONSIDERING A PROPER ANALYSIS ID FROM THE DATABASE. THE THREE VARIABLES TO BE CHECKED BEFORE RUNNING THE SCRIPT ARE $datatype1, $datatype2 AND $analysis_id WHICH ARE SUPPOSED TO BELONG TO THE SAME DATASET
use strict;
use warnings;
use Data::Dumper;

#USE THE SQL DATA MANAGEMENT MODULES
use DBI;
use DBD::mysql;

#MySQL CONFIG VARIABLES
my $platform  = "mysql";
my $host      = "localhost";
my $port      = "3306";
my $user      = "mysql_dev";
my $pw        = "dEvEl0pEr";
my $database  = "storehouse";
my $analysis_id = 11;
my $datatype1 = "ZebrafishCNS";
my $datatype2 = "MouseCNS";
#DATA SOURCE NAME
my $dsn = "dbi:$platform:$database:$host:$port";
#PERL DBI CONNECT AND CREATE TABLES
my $dbh =  DBI->connect($dsn, $user, $pw, { RaiseError => 1, AutoCommit => 0 } )  || die "Could not connect to database: $DBI::errstr";



#ACCESS THE QUERY IDENTIFIER SUITABLE TO THE CUT OFF
my $sthA = $dbh->selectall_arrayref("select distinct a1.name,a1.region,a1.start,a1.end,a2.name,a2.region,a2.start,a2.end,a1.strand, a2.strand from feature_pair inner join annotation as a1 on feature_pair.h_id = a1.ID inner join annotation as a2 on feature_pair.q_ID = a2.ID where feature_pair.analysis_id = $analysis_id and a1.datatype = '$datatype1' and a2.datatype = '$datatype2' AND a1.exon_overlap IS NOT NULL AND a2.exon_overlap IS NULL;"); #ADD (AND a2.strand = '-1') if strand specific results are to extracted
foreach my $row (@$sthA){

  my $qNAME = $row->[0];
  my $qREG = $row->[1];
  my $qSTART = $row->[2];
  my $qEND = $row->[3];

  my $hNAME = $row->[4];
  my $hREG = $row->[5];
  my $hSTART = $row->[6];
  my $hEND = $row->[7];
  my $qstrand = $row->[8];
  my $hstrand = $row->[9];
  
  #print $qNAME."\t".$qREG.":".$qSTART."-".$qEND."\t".$hNAME."\t".$hREG.":".$hSTART."-".$hEND."\n";

  my $sthB = $dbh->selectall_arrayref("select ftype from est where name = '$qNAME' and ftype = 'overlapping';");
  next unless exists $sthB->[0]->[0];
  print $qNAME."\t".$qREG.":".$qSTART."-".$qEND."\t".$qstrand."\t".$hNAME."\t".$hREG.":".$hSTART."-".$hEND."\t".$hstrand."\n";
    
}

#PERL DBI DATA DISCONNECT
$dbh->disconnect();










