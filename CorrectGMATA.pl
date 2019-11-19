#!/usr/bin/env perl

use strict;
use warnings;

my ($frgFile, $emapFile, $mkFile, $seqFile, $ssrFile) = @ARGV;

my %ssrMaker;
open my $iFRG, '<', $frgFile or die "can't open file $frgFile : $!";
while (<$iFRG>) {
    my ($name, $judge, $len) = split;
    next if $judge != 1;
    #my ($product, $theory) = $len =~ /\A(\d+)\/(\d+)/;
    #$ssrMaker{$name}{allele} = $product == $theory ? 'mono' : 'poly';
    $ssrMaker{$name} = undef;
}
close $iFRG;

open my $iEMAP, '<', $emapFile or die "can't open file $emapFile : $!";
while (<$iEMAP>) {
    my ($seqId, $name, $strand, $start, $end, undef, undef, undef) = split;
    if (exists $ssrMaker{$name}) {
        $ssrMaker{$name}{id}    = $seqId;
        $ssrMaker{$name}{start} = $start;
        $ssrMaker{$name}{end}   = $end;
    }
}
close $iEMAP;

my %info;
open my $iMK, '<', $mkFile or die "can't open file $mkFile : $!";
<$iMK>;
my @fileInfo = <$iMK>;
close $iMK;
my %mks;
foreach my $each (@fileInfo) {
    my (undef, $name, @noUse) = split /\s+/, $each;
    next if ! defined $name;
    $mks{$name}++;
}
foreach my $each (@fileInfo) {
    my ($ssrInfo, $name, @noUse) = split /\s+/, $each;
    if (defined $name and exists $ssrMaker{$name}) {
        next if $mks{$name} == 0; #Defect: if a maker appears multiple times, only keep the first time
        my ($seqId, undef, $ssrPos) = split '\|', $ssrInfo;
        #$info{"$seqId-$ssrPos"}{allele} = $ssrMaker{$name}{allele};
        $info{"$seqId-$ssrPos"}{start}  = $ssrMaker{$name}{start};
        $info{"$seqId-$ssrPos"}{end}    = $ssrMaker{$name}{end};
        $info{"$seqId-$ssrPos"}{id}     = $ssrMaker{$name}{id};
        $mks{$name} = 0 if $mks{$name} > 1;
    }
}
undef %mks;
undef %ssrMaker;

my %seqs;
open my $iSEQ, '<', $seqFile or die "can't open file $seqFile : $!";
my $id;
while (<$iSEQ>) {
	s/[\r\n]+//;
	next unless $_;
	/\A>(\S+)/ ? ( $id = $1 ) : ( s/[^a-zA-Z]//g, $seqs{$id} .= uc );
}
close $iSEQ;

my ($keepTotal, $keepPoly, $keepMono, $filter)   = (0, 0, 0, 0);

my %threshold = (1 => 10, 2 => 7, 3 => 6, 4 => 5, 5 => 4, 6 => 4);

open my $iSSR, '<', $ssrFile or die "can't open file $ssrFile : $!";
<$iSSR>;
while (<$iSSR>) {
    my ($seqId, undef, $start, $end, $repeat, $motif) = split;
    my $len = length($motif)*$repeat;
    my $key = "$seqId-$start:$end";
    if ($info{$key}) {
        my $seq = substr $seqs{$info{$key}{id}}, $info{$key}{start}-1, $info{$key}{end}-$info{$key}{start}+1;
        my $number = $threshold{length $motif};
        my ($ssr) = $seq =~ /(($motif){$number,})/; #Defect: if have more the same SSR, just keep the first one
        if (defined $ssr) {
            length($ssr) != $len ? $keepPoly++ : $keepMono++;
            $keepTotal++;
        }
        else {
            $filter++;
        }
    }
}
close $iSSR;

print "keep Total = $keepTotal\n";
print "keep Poly  = $keepPoly\n";
print "keep Mono  = $keepMono\n";
print "filter     = $filter\n";

__END__

#Glycine max
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_1/Glycine_max-ZH13.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_1/Glycine_max-ZH13.fa.eMap /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_1/Glycine_max-Williams82.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Glycine_max-ZH13.fa /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_1/Glycine_max-Williams82.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_2/Glycine_max-ZH13.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_2/Glycine_max-ZH13.fa.eMap /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_2/Glycine_max-Williams82.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Glycine_max-ZH13.fa /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_2/Glycine_max-Williams82.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_3/Glycine_max-ZH13.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_3/Glycine_max-ZH13.fa.eMap /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_3/Glycine_max-Williams82.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Glycine_max-ZH13.fa /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_3/Glycine_max-Williams82.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_4/Glycine_max-ZH13.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_4/Glycine_max-ZH13.fa.eMap /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_4/Glycine_max-Williams82.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Glycine_max-ZH13.fa /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_4/Glycine_max-Williams82.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_5/Glycine_max-ZH13.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_5/Glycine_max-ZH13.fa.eMap /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_5/Glycine_max-Williams82.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Glycine_max-ZH13.fa /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_5/Glycine_max-Williams82.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_6/Glycine_max-ZH13.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_6/Glycine_max-ZH13.fa.eMap /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_6/Glycine_max-Williams82.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Glycine_max-ZH13.fa /lthpcfs/home/gouxj/SSR/output/soybean/GMATA/Glycine_max_6/Glycine_max-Williams82.fa.ssr

#Gossypium hirsutum
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum1/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum1/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum1/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum1/Gossypium_hirsutum-TM1.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum2/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum2/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum2/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum2/Gossypium_hirsutum-TM1.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum3/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum3/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum3/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum3/Gossypium_hirsutum-TM1.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum4/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum4/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum4/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum4/Gossypium_hirsutum-TM1.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum5/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum5/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum5/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum5/Gossypium_hirsutum-TM1.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum6/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum6/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum6/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum6/Gossypium_hirsutum-TM1.fa.ssr

#Triticum aestivum
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum1/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum1/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum1/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum1/Triticum_aestivum-CS-7D.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum2/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum2/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum2/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum2/Triticum_aestivum-CS-7D.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum3/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum3/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum3/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum3/Triticum_aestivum-CS-7D.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum4/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum4/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum4/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum4/Triticum_aestivum-CS-7D.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum5/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum5/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum5/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum5/Triticum_aestivum-CS-7D.fa.ssr
perl CorrectGMATA.pl /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum6/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum6/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum6/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum6/Triticum_aestivum-CS-7D.fa.ssr
