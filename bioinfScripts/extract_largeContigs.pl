
use warnings;
use strict;


use Bio::SeqIO::fasta;
use Bio::SearchIO;

my $inputData=$ARGV[0];
my $outputdir=$ARGV[1];
my $kmer=$ARGV[2];
my $cutoffLength=$ARGV[3];
my $numArgs = $#ARGV + 1;

print "$numArgs,  Output directory is $outputdir, input dataset is $inputData, kmer chosen is $kmer \n";

open(IN,  "$inputData") or die "cannot open input file\n";


my $resFile = ${outputdir}."/output".$kmer."_contigs_grt150.txt";
open (OUT2, " > $resFile") or die "cannot open output file\n";
print "Output in $resFile\n";


my $resFile2 = ${outputdir}."/output".$kmer."_contigs_grt150.fa";
open (OUT3, " > $resFile2") or die "cannot open output file\n";
print "FASTA Output in $resFile2\n";


my $resFile3 = ${outputdir}."/output".$kmer."_node_ids.txt";
open (OUT4, " > $resFile3") or die "cannot open output file\n";
print "Contig ids in $resFile3\n";

my $resFile4 = ${outputdir}."/".$kmer.".txt";
open (OUT5, " > $resFile4") or die "cannot open output file\n";
print "Output in $resFile4\n";


my $seqio = Bio::SeqIO->new(-file=> $inputData,-format => 'Fasta' );

my $sum_contig_length=0;
while ( my $seq_obj = $seqio->next_seq()){
    my $read_id=$seq_obj->id;
    my $read_seq=$seq_obj->seq;
     
    
    my $length = (split /\_/, $read_id)[-3];
    my $contigID = (split /\_/, $read_id)[-5];

    my $kvalue = substr($kmer, 1, 2);

    my $realLength = $length+$kvalue-1 ;
    

    if ($realLength>=$cutoffLength){
	print OUT2 $read_id, "\t", $read_seq, "\n";      
	print OUT3 ">",$read_id, "\n", $read_seq, "\n";
	print OUT4 $contigID, "\t", $read_id, "\n";
	$sum_contig_length+=$realLength;
    }
}
print OUT5  $kmer, "\t", $sum_contig_length, "\n";





