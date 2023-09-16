#!/usr/bin/perl

use warnings;
use strict;  
use File::Basename;
use List::MoreUtils qw(uniq);   

sub program_info {
    print "\n\t--------------------\n\trMATS_delta_IncDiff_inverter.pl will generate a new file in which\n\tall dPSI values are reversed (i.e. positive values are converted to negative\n\tvalues and vice versa)\n\n\t- This obviates the need to re-run rMATS under reversed conditions\n\t  (i.e. test vs control or control vs test).\n\n\t- Useful for SpliceCompare and Translate_NMD modules\n\n\t- File names must be appended with one of the following:\n\n\t\tA3SS.MATS.JCEC.txt\n\t\tA5SS.MATS.JCEC.txt\n\t\tMXE.MATS.JCEC.txt\n\t\tRI.MATS.JCEC.txt\n\t\tSE.MATS.JCEC.txt\n\n\tUsage: perl rMATS_delta_IncDiff_inverter.pl rMATS_file(s)\n\t--------------------\n\n";
    exit;
}

sub qc {
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }    
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] =~ m/A3SS\.MATS\.JCEC\.txt$/ || $ARGV[$i] =~ m/A5SS\.MATS\.JCEC\.txt$/ || $ARGV[$i] =~ m/MXE\.MATS\.JCEC\.txt$/ || $ARGV[$i] =~ m/RI\.MATS\.JCEC\.txt$/ || $ARGV[$i] =~ m/SE\.MATS\.JCEC\.txt$/) {
        }
        else {
            my $file = basename($ARGV[$i]);
            print "\n\n", $file, " doesn't match rMATS file extension!!\tEXITING PROGRAM!\n";
            program_info;
            exit;
        }
    }
}

sub invert {
    foreach my $file(@ARGV) {
        my $file_name = basename($file);
        my $directory_name = dirname($file);
        my @split_file_name = split("\_", $file_name);
        my $splice_type_appendage = pop(@split_file_name);

        my $output_directory_file_name = $directory_name."/".join("\_", @split_file_name)."_inverted_dPSI_".$splice_type_appendage;
        open(INF, "<$file") or die "couldn't open input file";
        open(OUT, ">$output_directory_file_name") or die "couldn't open output file";

        if ($file =~ m/MXE.MATS.JCEC.txt/g) {
            while(my $line = <INF>) {
                chomp($line);
                if ($. == 1) {
                    print OUT $line, "\n";
                }
                else {
                    my @split_line = split("\t", $line);
                    print OUT join("\t", @split_line[0..23]), "\t", -$split_line[24], "\n";
                }
            }
        }
        else {
            while(my $line = <INF>) {
                chomp($line);
                if ($. == 1) {
                    print OUT $line, "\n";
                }
                else {
                    my @split_line = split("\t", $line);
                    print OUT join("\t", @split_line[0..21]), "\t", -$split_line[22], "\n";
                }
            }
        }
        close(INF);
        close(OUT);
    }
}
qc;
invert;