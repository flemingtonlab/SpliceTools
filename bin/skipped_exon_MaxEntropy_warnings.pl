# Notes: 
    # Takes skipped exon file from rMATS and determines MaxEntropy (Burge lab) splice site scores for upstream splice donors, SE acceptor, SE donor and downstream acceptor for negative dpsis and for positive dpsis. Also determines these values for all candidate skipped exon events from every 10th annotated transcripts from annotation file (only every 10th in order to decrease processing time).
# Additional Note:
    # Annotation file must be in strange format (script must be adjusted to accomodate annotation file in normal configuration)
#
#USAGE:
	# perl skipped_exon_MaxEntropy.pl $skipped_exon_file $FDR $genome_fasta_file $annotation_file 
   
# EKF (Flemington lab) 06/12/2021.

use warnings;
use strict;  
#use diagnostics;
use File::Basename;

my ($skipped_exon_file, $FDR, $genome_fasta_file, $annotation_file) = (@ARGV);

my $script_dirname = dirname(__FILE__);
my @genome_array = ();
my $genome_array_length;
my @annotation_array = ();
my $directory;
my $basename;

sub makedirectories {
    $directory = dirname($skipped_exon_file);
    $basename = basename($skipped_exon_file);
    `mkdir $directory"\/1_output_"$basename"_FDR_"$FDR`;
    `mkdir $directory"\/1_output_"$basename"_FDR_"$FDR"\/1_neg_dpsi"`;
    `mkdir $directory"\/1_output_"$basename"_FDR_"$FDR"\/2_pos_dpsi"`;
    `mkdir $directory"\/1_output_"$basename"_FDR_"$FDR"\/3_annotated"`;
}

sub genome_fasta_file_unwrapper {
    print "\nUnwrapping genome fasta file and putting into memory...\n\n";
    open(INF, "<$genome_fasta_file") or die "couldn't open genome fasta file";
    my @line_array = ();
    while(my $line = <INF>) {
        chomp($line);
        if ($. == 1) {
            push(@genome_array, $line);
        }
        elsif ($line =~ m/^\>chr/ || $line =~ m/^\> chr/) {
            push(@genome_array, join("", @line_array));
            push(@genome_array, $line);
            @line_array = ();
        }
        elsif (eof(INF)) {
            push(@line_array, $line);
            push(@genome_array, @line_array);
        }
        else {  
            push(@line_array, $line);
        }
    }
    close(INF);
    $genome_array_length = @genome_array;
}

sub SE_scoring {
    print "Calculating MaxEntropy scores for skipped exons...\n\n";
    my $upstream_exon_donor_neg_dpsi = "\>upstream_exon_donor_neg_dpsi\n";
    my $skipped_exon_acceptor_neg_dpsi = "\>skipped_exon_acceptor_neg_dpsi\n";
    my $skipped_exon_donor_neg_dpsi = "\>skipped_exon_donor_neg_dpsi\n";
    my $downstream_exon_acceptor_neg_dpsi = "\>downstream_exon_acceptor_neg_dpsi\n";

    my $upstream_exon_donor_pos_dpsi = "\>upstream_exon_donor_pos_dpsi\n";
    my $skipped_exon_acceptor_pos_dpsi = "\>skipped_exon_acceptor_pos_dpsi\n";
    my $skipped_exon_donor_pos_dpsi = "\>skipped_exon_donor_pos_dpsi\n";
    my $downstream_exon_acceptor_pos_dpsi = "\>downstream_exon_acceptor_pos_dpsi\n";

    open(INF, "<$skipped_exon_file") or die "couldn't open skipped exon input file";
    while(my $line = <INF>) {
        chomp($line);
        next if ($. == 1);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        if ($split_line[19] <= $FDR)  {
            my $strand = $split_line[4];
            my $amplitude = $split_line[22];
            my $chr = $split_line[3];
            if ($strand eq "+") {
                my $chr_hit = 0;
                for (my $i = 0; $i < $genome_array_length; $i++) {
                    if ($chr_hit == 1) {
                        if ($amplitude < 0) {
                            $upstream_exon_donor_neg_dpsi = $upstream_exon_donor_neg_dpsi.uc(substr $genome_array[$i], $split_line[8]-3, 9)."\n";
                            $skipped_exon_acceptor_neg_dpsi = $skipped_exon_acceptor_neg_dpsi.uc(substr $genome_array[$i], $split_line[5]-20, 23)."\n";
                            $skipped_exon_donor_neg_dpsi = $skipped_exon_donor_neg_dpsi.uc(substr $genome_array[$i], $split_line[6]-3, 9)."\n";
                            $downstream_exon_acceptor_neg_dpsi = $downstream_exon_acceptor_neg_dpsi.uc(substr $genome_array[$i], $split_line[9]-20, 23)."\n";
                        }
                        elsif ($amplitude > 0) {
                            $upstream_exon_donor_pos_dpsi = $upstream_exon_donor_pos_dpsi.uc(substr $genome_array[$i], $split_line[8]-3, 9)."\n";
                            $skipped_exon_acceptor_pos_dpsi = $skipped_exon_acceptor_pos_dpsi.uc(substr $genome_array[$i], $split_line[5]-20, 23)."\n";
                            $skipped_exon_donor_pos_dpsi = $skipped_exon_donor_pos_dpsi.uc(substr $genome_array[$i], $split_line[6]-3, 9)."\n";
                            $downstream_exon_acceptor_pos_dpsi = $downstream_exon_acceptor_pos_dpsi.uc(substr $genome_array[$i], $split_line[9]-20, 23)."\n";
                        }
                        last;
                    }
                    elsif ($genome_array[$i] =~ m/^\>/) {
                        my $genome_array_chromosome = $genome_array[$i];
                        $genome_array_chromosome =~ s/ //g;
                        $genome_array_chromosome =~ s/\>//g;
                        if ($genome_array_chromosome eq $chr) {
                            $chr_hit++;
                        }
                    }
                }
            }
            elsif ($strand eq "-") {
                my $chr_hit = 0;
                for (my $i = 0; $i < $genome_array_length; $i++) {
                    if ($chr_hit == 1) {
                        my $upstream_donor = uc(substr $genome_array[$i], $split_line[9]-6, 9);
                        $upstream_donor =~ tr/AGCT/TCGA/;
                        $upstream_donor = reverse($upstream_donor);
                        my $SE_acceptor = uc(substr $genome_array[$i], $split_line[6]-3, 23);
                        $SE_acceptor =~ tr/AGCT/TCGA/;
                        $SE_acceptor = reverse($SE_acceptor);
                        my $SE_donor = uc(substr $genome_array[$i], $split_line[5]-6, 9);
                        $SE_donor =~ tr/AGCT/TCGA/;
                        $SE_donor = reverse($SE_donor);
                        my $downstream_acceptor = uc(substr $genome_array[$i], $split_line[8]-3, 23);
                        $downstream_acceptor =~ tr/AGCT/TCGA/;
                        $downstream_acceptor = reverse($downstream_acceptor);

                        if ($amplitude < 0) {
                            $upstream_exon_donor_neg_dpsi = $upstream_exon_donor_neg_dpsi.$upstream_donor."\n";
                            $skipped_exon_acceptor_neg_dpsi = $skipped_exon_acceptor_neg_dpsi.$SE_acceptor."\n";
                            $skipped_exon_donor_neg_dpsi = $skipped_exon_donor_neg_dpsi.$SE_donor."\n";
                            $downstream_exon_acceptor_neg_dpsi = $downstream_exon_acceptor_neg_dpsi.$downstream_acceptor."\n";
                        }
                        elsif ($amplitude > 0) {
                            $upstream_exon_donor_pos_dpsi = $upstream_exon_donor_pos_dpsi.$upstream_donor."\n";
                            $skipped_exon_acceptor_pos_dpsi = $skipped_exon_acceptor_pos_dpsi.$SE_acceptor."\n";
                            $skipped_exon_donor_pos_dpsi = $skipped_exon_donor_pos_dpsi.$SE_donor."\n";
                            $downstream_exon_acceptor_pos_dpsi = $downstream_exon_acceptor_pos_dpsi.$downstream_acceptor."\n";
                        }
                        last;
                    }
                    elsif ($genome_array[$i] =~ m/^\>/) {
                        my $genome_array_chromosome = $genome_array[$i];
                        $genome_array_chromosome =~ s/ //g;
                        $genome_array_chromosome =~ s/\>//g;
                        if ($genome_array_chromosome eq $chr) {
                            $chr_hit++;
                        }
                    }
                }
            }
        }
    }
    close(INF);

    my $neg_dpsi_directory = $directory."\/1_output_".$basename."_FDR_".$FDR."\/1_neg_dpsi\/out";
    my $pos_dpsi_directory = $directory."\/1_output_".$basename."_FDR_".$FDR."\/2_pos_dpsi\/out";

    open(OUT1, ">$neg_dpsi_directory.upstream_donor_neg") or die "couldn't open upstream_donor_neg file";
    open(OUT2, ">$neg_dpsi_directory.SE_acceptor_neg") or die "couldn't open SE_acceptor_neg file"; 
    open(OUT3, ">$neg_dpsi_directory.SE_donor_neg") or die "couldn't open SE_donor_neg file"; 
    open(OUT4, ">$neg_dpsi_directory.downstream_acceptor_neg") or die "couldn't open downstream_acceptor_neg file";  
    open(OUT5, ">$pos_dpsi_directory.upstream_donor_pos") or die "couldn't open upstream_donor_pos file";
    open(OUT6, ">$pos_dpsi_directory.SE_acceptor_pos") or die "couldn't open SE_acceptor_pos file"; 
    open(OUT7, ">$pos_dpsi_directory.SE_donor_pos") or die "couldn't open SE_donor_pos file"; 
    open(OUT8, ">$pos_dpsi_directory.downstream_acceptor_pos") or die "couldn't open downstream_acceptor_pos file";  

    print OUT1 $upstream_exon_donor_neg_dpsi;
    print OUT2 $skipped_exon_acceptor_neg_dpsi;
    print OUT3 $skipped_exon_donor_neg_dpsi;
    print OUT4 $downstream_exon_acceptor_neg_dpsi;
    print OUT5 $upstream_exon_donor_pos_dpsi;
    print OUT6 $skipped_exon_acceptor_pos_dpsi;
    print OUT7 $skipped_exon_donor_pos_dpsi;
    print OUT8 $downstream_exon_acceptor_pos_dpsi;

    close(OUT1);
    close(OUT2);
    close(OUT3);
    close(OUT4);
    close(OUT5);
    close(OUT6);
    close(OUT7);
    close(OUT8);

    print OUTsummary "\tupstream exon donor\tSE acceptor\tSE donor\tdownstream exon acceptor\n";
    print OUTsummary "neg dpsi\t", score5($neg_dpsi_directory.".upstream_donor_neg"), "\t", score3($neg_dpsi_directory.".SE_acceptor_neg"), "\t", score5($neg_dpsi_directory.".SE_donor_neg"), "\t", score3($neg_dpsi_directory.".downstream_acceptor_neg"), "\npos dpsi\t", score5($pos_dpsi_directory.".upstream_donor_pos"), "\t", score3($pos_dpsi_directory.".SE_acceptor_pos"), "\t", score5($pos_dpsi_directory.".SE_donor_pos"), "\t", score3($pos_dpsi_directory.".downstream_acceptor_pos"), "\n";
}

sub prepare_annotation_file {
    print "Sorting annotation file and putting into array...\n\n";
    `sort -k1,1 -k2,2n -k3,3n $annotation_file > $annotation_file.sorted.tsv`;

    open(INF, "<$annotation_file.sorted.tsv") or die "couldn't open sorted annotation file";
    my $count = 0;
    while(my $line = <INF>) {
        if ($count == 10)  {
            chomp($line);
            $line =~ s/\"//g;
            my @split_line = split("\t", $line);
            my $fixed_line = join("\t", @split_line[0..3])."_".$split_line[4]."\t1000\t".join("\t", @split_line[5..11]); 
            push(@annotation_array, $fixed_line);
            $count = 0;
        }
        else {
            $count++;
        }
    }
    `rm $annotation_file.sorted.tsv`;
    close(INF);
}

sub annotated_scoring {
    print "Calculating MaxEntropy scores for annotated exons...\n\n";
    my $upstream_donors = "\>upstream_donors\n";
    my $middle_acceptors = "\>middle_acceptors\n";
    my $middle_donors = "\>middle_donors\n";
    my $downstream_acceptors = "\>downstream_acceptors\n";
    
    my $annotation_array_size = @annotation_array;
    for (my $j = 0; $j < $annotation_array_size; $j++) {
        my $annotation_array_line = $annotation_array[$j];
        my @split_annotation_array_line = split("\t", $annotation_array_line);
        my $transcript_start = $split_annotation_array_line[1];
        my @block_sizes = split("\,", $split_annotation_array_line[10]);#
        my @block_starts = split("\,", $split_annotation_array_line[11]);#
        my $block_number = $split_annotation_array_line[9];#
        my $strand = $split_annotation_array_line[5];
        my $chr = $split_annotation_array_line[0];

        my $chr_hit = 0;
        for (my $i =  0; $i < $genome_array_length; $i++) {
            if ($chr_hit == 1) {
                if ($strand eq "+") {
                    if ($block_number >2) {
                        for (my $z = 2; $z < $block_number; $z++) {
                            $upstream_donors = $upstream_donors.uc(substr $genome_array[$i], $transcript_start+$block_starts[$z-2]+$block_sizes[$z-2]-3, 9)."\n";
                            $middle_acceptors = $middle_acceptors.uc(substr $genome_array[$i], $transcript_start+$block_starts[$z-1]-20, 23)."\n";
                            $middle_donors = $middle_donors.uc(substr $genome_array[$i], $transcript_start+$block_starts[$z-1]+$block_sizes[$z-1]-3, 9)."\n";
                            $downstream_acceptors = $downstream_acceptors.uc(substr $genome_array[$i], $transcript_start+$block_starts[$z]-20, 23)."\n";
                        }
                    }
                }
                elsif ($strand eq "-") {
                    if ($block_number >2) {
                        for (my $z = 2; $z < $block_number; $z++) {
                            my $upstream_donor = uc(substr $genome_array[$i], $transcript_start+$block_starts[$z]-6, 9);
                            $upstream_donor =~ tr/AGCT/TCGA/;
                            $upstream_donor = reverse($upstream_donor);
                            my $middle_acceptor = uc(substr $genome_array[$i], $transcript_start+$block_starts[$z-1]+$block_sizes[$z-1]-3, 23);
                            $middle_acceptor =~ tr/AGCT/TCGA/;
                            $middle_acceptor = reverse($middle_acceptor);
                            my $middle_donor = uc(substr $genome_array[$i], $transcript_start+$block_starts[$z-1]-6, 9);
                            $middle_donor =~ tr/AGCT/TCGA/;
                            $middle_donor = reverse($middle_donor);
                            my $downstream_acceptor = uc(substr $genome_array[$i], $transcript_start+$block_starts[$z-2]+$block_sizes[$z-2]-3, 23);
                            $downstream_acceptor =~ tr/AGCT/TCGA/;
                            $downstream_acceptor = reverse($downstream_acceptor);

                            $upstream_donors = $upstream_donors.$upstream_donor."\n";
                            $middle_acceptors = $middle_acceptors.$middle_acceptor."\n";
                            $middle_donors = $middle_donors.$middle_donor."\n";
                            $downstream_acceptors = $downstream_acceptors.$downstream_acceptor."\n";
                        }
                    }
                }
                last;
            }
            elsif ($genome_array[$i] =~ m/^\>/) {
                my $genome_array_chromosome = $genome_array[$i];
                $genome_array_chromosome =~ s/ //g;
                $genome_array_chromosome =~ s/\>//g;
                if ($genome_array_chromosome eq $chr) {
                    $chr_hit++;
                }
            }
        }
    }

    my $annotated_directory = $directory."\/1_output_".$basename."_FDR_".$FDR."\/3_annotated\/out";

    open(OUT1, ">$annotated_directory.upstream_donor") or die "couldn't open upstream_donor file";
    open(OUT2, ">$annotated_directory.middle_acceptor") or die "couldn't open SE_acceptor file"; 
    open(OUT3, ">$annotated_directory.middle_donor") or die "couldn't open SE_donor file"; 
    open(OUT4, ">$annotated_directory.downstream_acceptor") or die "couldn't open downstream_acceptor file";  

    print OUT1 $upstream_donors;
    print OUT2 $middle_acceptors;
    print OUT3 $middle_donors;
    print OUT4 $downstream_acceptors;

    close(OUT1);
    close(OUT2);
    close(OUT3);
    close(OUT4);

    print OUTsummary "annotated\t", score5($annotated_directory.".upstream_donor"), "\t", score3($annotated_directory.".middle_acceptor"), "\t", score5($annotated_directory.".middle_donor"), "\t", score3($annotated_directory.".downstream_acceptor"); 
}

sub score3 {
    my $score3 = $script_dirname."\/score3_EKF.pl";
    foreach my $input(@_) {
        my $output_file = $input.".scores.txt";
        `perl $score3 $input > $output_file`;
        `rm $input`;
        open(INF, "<$output_file") or die "couldn't open input file";
        my @score_array = ();
        while(my $line = <INF>) {
            chomp($line);
            my @split_line = split("\t", $line);
            push(@score_array, $split_line[1]);
        }
        close(INF);
        my $average = average_of_array(@score_array);
        return($average);
    }
}

sub score5 {
    my $score5 = $script_dirname."\/score5_EKF.pl";
    foreach my $input(@_) {
        my $output_file = $input.".scores.txt";
        `perl $score5 $input > $output_file`;
        `rm $input`;
        open(INF, "<$output_file") or die "couldn't open input file";
        my @score_array = ();
        while(my $line = <INF>) {
            chomp($line);
            my @split_line = split("\t", $line);
            push(@score_array, $split_line[1]);
        }
        close(INF);
        my $average = average_of_array(@score_array);
        return($average);
    }
}

sub average_of_array {
    my $array_length = @_;
    my $total = 0;
    for my $element(@_) {
        $total += $element;
    }
    return $total/$array_length;
}

makedirectories;
my $summary_output_file = $directory."\/1_output_".$basename."_FDR_".$FDR."\/1_1_summary.tsv";
open(OUTsummary, ">$summary_output_file") or die "couldn't open upstream_donor_neg file";
genome_fasta_file_unwrapper;
SE_scoring;
if (-e $annotation_file) {
    prepare_annotation_file;
    annotated_scoring;
}
close(OUTsummary);