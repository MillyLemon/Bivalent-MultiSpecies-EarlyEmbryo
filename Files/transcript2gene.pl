#!/usr/bin/perl -w
use strict;

my $file = shift;
open F,"gzip -cd $file|" or die $!;
while (<F>){
    chomp;next if (/#/); 
    my @a = split /\t/; # 1       ensembl gene    226161299       226217308       .       -       .       gene_id "ENSSSCG00000028996"; gene_version "4"; gene_name "ALDH1A1"; gene_source "ensembl"; gene_biotype "protein_coding";
    #gene_id "ENSSSCG00000005273"; gene_version "5"; transcript_id "ENSSSCT00000042535"; transcript_version "2"; gene_name "OSTF1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "OSTF1-202"; transcript_source "ensembl"; transcript_biotype "protein_coding";
    #gene_id "ENSSSCG00000045514"; gene_version "1"; transcript_id "ENSSSCT00000072157"; transcript_version "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding";
# chr1    HAVANA  transcript      3073253 3074322 .       +       .       gene_id "ENSMUSG00000102693.1"; transcript_id "ENSMUST00000193812.1"; gene_type "TEC"; gene_name "4933401J01Rik"; transcript_type "TEC"; transcript_name "4933401J01Rik-201"; level 2; transcript_support_level "NA"; mgi_id "MGI:1918292"; tag "basic"; havana_gene "OTTMUSG00000049935.1"; havana_transcript "OTTMUST00000127109.1";
#
    if ($a[2] eq "transcript"){
        my @b = split /\s+/,$a[8];
        my $line; 
        if(/gene_name/){
            $line = "$b[1]\t$b[3]\t$b[7]"; #print "$line\n"; 
        }
        my @c = split /"/,$line;
        my $info = "$c[3]\t$c[1]\t$c[5]";
        print "$info\n";
    }
}
close F;
