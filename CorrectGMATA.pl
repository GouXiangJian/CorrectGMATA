#!/usr/bin/env perl

#date : 2020-03-25

use strict;
use warnings;
use File::Basename;

#get all input
my ($motifType, $frgFile, $emapFile, $mkFile, $seqFile, $ssrFile) = @ARGV;

#set SSR threshold
my %threshold = (1 => 10, 2 => 7, 3 => 6, 4 => 5, 5 => 4, 6 => 4);

#===================================================

my %ssrMaker;
open my $iFRG, '<', $frgFile or die "can't open file $frgFile : $!";
while (<$iFRG>) {
    my ($name, $judge, $len) = split;
    next if $judge != 1;
    my ($product, $theory) = $len =~ /\A(\d+)\/(\d+)/;
    $ssrMaker{$name}{allele} = $product == $theory ? 'mono' : 'poly';
    $ssrMaker{$name}{other}  = undef;
}
close $iFRG;

#===================================================

open my $iEMAP, '<', $emapFile or die "can't open file $emapFile : $!";
while (<$iEMAP>) {
    my ($seqId, $name, $strand, $start, $end, undef, undef, undef) = split;
    if (exists $ssrMaker{$name}) {
        $ssrMaker{$name}{other}{id}     = $seqId;
        $ssrMaker{$name}{other}{start}  = $start;
        $ssrMaker{$name}{other}{end}    = $end;
        $ssrMaker{$name}{other}{strand} = $strand;
    }
}
close $iEMAP;

#===================================================

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
        next if $mks{$name} == 0; #if a maker appears multiple times, only keep the first time
        my ($seqId, undef, $ssrPos) = split '\|', $ssrInfo;
        $info{"$seqId-$ssrPos"}{allele} = $ssrMaker{$name}{allele};
        $info{"$seqId-$ssrPos"}{start}  = $ssrMaker{$name}{other}{start};
        $info{"$seqId-$ssrPos"}{end}    = $ssrMaker{$name}{other}{end};
        $info{"$seqId-$ssrPos"}{id}     = $ssrMaker{$name}{other}{id};
        $info{"$seqId-$ssrPos"}{strand} = $ssrMaker{$name}{other}{strand};
        $mks{$name} = 0 if $mks{$name} > 1;
    }
}
undef %mks;
undef %ssrMaker;

#===================================================

my %seqs;
open my $iSEQ, '<', $seqFile or die "can't open file $seqFile : $!";
my $id;
while (<$iSEQ>) {
	s/[\r\n]+//;
	/\A>(\S+)/ ? ( $id = $1 ) : ( s/[^a-zA-Z]//g, $seqs{$id} .= uc );
}
close $iSEQ;

#===================================================

open my $iSSR, '<', $ssrFile or die "can't open file $ssrFile : $!";
my $outName = basename $ssrFile;
open my $oSSR, '>', "$outName.correct.$motifType";
my @first = split /\s+/, <$iSSR>;
print $oSSR join("\t", @first, 'old', 'new'), "\n";

while (<$iSSR>) {
    my @row = split;
    my ($seqId, $seqLen, $start, $end, $repeat, $motif) = @row;
    my $len = length($motif)*$repeat;
    my $key = "$seqId-$start:$end";
    if ($info{$key}) {
        my $seq = substr $seqs{$info{$key}{id}}, $info{$key}{start}-1, $info{$key}{end}-$info{$key}{start}+1;
        if ($info{$key}{strand} eq '-') {
            $seq = reverse $seq;
            $seq =~ tr/AGTC/TCAG/;
        }
        my $number = $threshold{length $motif};
        my (@ssr) = $seq =~ /((?:$motif){$number,})/g;
        if (@ssr == 0) {
            print $oSSR join("\t", @row, $info{$key}{allele}, 'be-filtered'), "\n";
        }
        else {
            my $mark = 0;
            foreach my $ssr (@ssr) {
                if (length($ssr) != $len) {
                    $mark = 1;
                    last;
                }
            }
            if ($mark == 1) {
                print $oSSR join("\t", @row, $info{$key}{allele}, 'poly'), "\n";
            }
            else {
                print $oSSR join("\t", @row, $info{$key}{allele}, 'mono'), "\n";
            }
        }
    }
}

close $iSSR;
close $oSSR;

__END__

#Gossypium hirsutum
perl CorrectGMATA.pl 1 /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum1/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum1/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum1/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum1/Gossypium_hirsutum-TM1.fa.ssr
perl CorrectGMATA.pl 2 /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum2/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum2/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum2/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum2/Gossypium_hirsutum-TM1.fa.ssr
perl CorrectGMATA.pl 3 /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum3/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum3/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum3/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum3/Gossypium_hirsutum-TM1.fa.ssr
perl CorrectGMATA.pl 4 /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum4/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum4/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum4/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum4/Gossypium_hirsutum-TM1.fa.ssr
perl CorrectGMATA.pl 5 /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum5/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum5/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum5/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum5/Gossypium_hirsutum-TM1.fa.ssr
perl CorrectGMATA.pl 6 /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum6/Gossypium_hirsutum-ZM24.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum6/Gossypium_hirsutum-ZM24.fa.eMap /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum6/Gossypium_hirsutum-TM1.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Gossypium_hirsutum-ZM24.fa /lthpcfs/home/gouxj/SSR/output/cotton/GMATA/Gossypium_hirsutum6/Gossypium_hirsutum-TM1.fa.ssr

#Triticum aestivum
perl CorrectGMATA.pl 1 /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum1/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum1/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum1/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum1/Triticum_aestivum-CS-7D.fa.ssr
perl CorrectGMATA.pl 2 /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum2/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum2/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum2/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum2/Triticum_aestivum-CS-7D.fa.ssr
perl CorrectGMATA.pl 3 /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum3/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum3/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum3/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum3/Triticum_aestivum-CS-7D.fa.ssr
perl CorrectGMATA.pl 4 /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum4/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum4/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum4/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum4/Triticum_aestivum-CS-7D.fa.ssr
perl CorrectGMATA.pl 5 /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum5/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum5/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum5/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum5/Triticum_aestivum-CS-7D.fa.ssr
perl CorrectGMATA.pl 6 /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum6/Triticum_aestivum-AK58.fa.eMap.frg /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum6/Triticum_aestivum-AK58.fa.eMap /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum6/Triticum_aestivum-CS-7D.fa.ssr.mk /lthpcfs/home/gouxj/SSR/sequence/Triticum_aestivum-AK58.fa /lthpcfs/home/gouxj/SSR/output/wheat_poly/GMATA/Triticum_aestivum6/Triticum_aestivum-CS-7D.fa.ssr
