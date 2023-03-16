use warnings;
use strict;  
use File::Basename;

my $skipped_exon_file;
my $annotation_file;
my $FDR;
my @skipped_exon_array = ();
my @annotation_array = ();
my @merged_SE_ann_array = ();

my $outdir;
my $basename;

sub program_info {
    print "\n\tSEUnannotated.pl will generate statistics on the number of skipping events that are present\n\t(annotated) and not present (unannotated) in a bed12 input annotation file. Also outputs\n\tlists of annotated and unannotated events (pos IncDiff and neg IncDiff).\n\n\tUsage: perl SEUnannotated.pl [OPTIONS] -s <skipped exon file (rMATS)> -a <bed12 annotation file>  -f <FDR>\n\n\tRequired:\n\t\t-s <skipped exon file>\n\t\t-a <bed12 annotation file>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl SEUnannotated.pl -s PATH/SEfile.txt -a PATH/bed12_annotation.bed -f 0.05\n\n";
    exit;
}

sub options {
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }    
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] eq "\-s") {
            $skipped_exon_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-a") {
            $annotation_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-f") {
            $FDR = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-h") {
            program_info;
            exit;
        }
    }
}

sub qc {
    if (not defined($FDR)) {
        print "\nFDR not defined!\n";
        program_info;
        exit;
    }
    elsif ($FDR !~ m/^\d/) {
        print "\nFDR is not numeric!\n";
        program_info;
        exit;
    } 
    elsif (not defined($skipped_exon_file)) {
        print "\nSkipped exon file not defined!\n";
        program_info;
        exit;
    }
    open(INF, "<$skipped_exon_file") or die "couldn't open skipped exon file";
    while(my $line = <INF>) {
        chomp($line);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        if ($. == 2) {
            
            if (scalar @split_line != 23 || $split_line[4] !~ m/[+-]/ || $split_line[5] !~ m/^\d/ || $split_line[6] !~ m/^\d/ || $split_line[7] !~ m/^\d/ || $split_line[8] !~ m/^\d/ || $split_line[9] !~ m/^\d/ || $split_line[10] !~ m/^\d/) {
                print "\nSkipped exon file may not be an rMATS JCEC file format!\n\n";
                program_info;
                exit;
            }
            elsif ($split_line[3] !~ m/^chr/) {
                print "\nSkipped exon file may not be an rMATS JCEC file format!\n     Is fourth column chromosome? Fourth column doesn't begin with \"chr\"\n\n";
                program_info;
                exit;
            }
        }
        elsif ($. == 3) {
            last;
        }
    }
    close(INF);
    if ($annotation_file ne "none") {   
        open(INF, "<$annotation_file") or die "couldn't open annotation file";
        while(my $line = <INF>) {
            chomp($line);
            my @split_line = split("\t", $line);
            if ($. == 2) {
                if (scalar @split_line != 12 || $split_line[0] !~ m/^chr/ || $split_line[5] !~ m/[+-]/ || $split_line[6] !~ m/^\d/ || $split_line[7] !~ m/^\d/ || $split_line[8] !~ m/^\d/ || $split_line[9] !~ m/^\d/ || $split_line[10] !~ m/^\d/ || $split_line[11] !~ m/^\d/) {
                    print "\nAnnotation file may not be in bed12 format!\n\n";
                    program_info;
                    exit;
                }
            }
            elsif ($. == 3) {
                last;
            }
        }
        close(INF);  
    }
}

sub makedirectories {
    $basename = basename($skipped_exon_file);
    $basename =~ s/\.txt$//g;
    $outdir = dirname($skipped_exon_file)."\/SEUnannotated_".$basename."_FDR_".$FDR;
    `mkdir $outdir`;
}

sub print_input_parameters {
    open(OUT, ">$outdir\/1_1_input_parameters.txt") or die "couldn't open run parameters output file";

    print OUT "Skipped exon file: ", $basename, "\nAnnotation file: ", basename($annotation_file), "\nFDR: ", $FDR, "\n\n";
    print "\nINPUT\nSkipped exon file: ", $basename, "\nAnnotation file: ", basename($annotation_file), "\nFDR: ", $FDR, "\n\n";
    close(OUT);
}

sub prepare_skipped_exon_array {
    print "Processinging skipped exon file...\n\n";
    open(INF, "<$skipped_exon_file") or die "couldn't open skipped exon input file";
    open(OUT, ">$skipped_exon_file.temp") or die "couldn't open temporary skipped exon output file";
    while(my $line = <INF>) {
        chomp($line);
        next if ($. == 1);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        if ($split_line[19] <= $FDR) {
            print OUT $split_line[3], "\t", $split_line[8], "\t", $split_line[9], "\t", $split_line[2],  "\t",$split_line[22], "\t", $split_line[4],  "\n";
        }
    }
    close(INF);
    close(OUT);

    `sort -k1,1 -k2,2n -k3,3n $skipped_exon_file.temp > $skipped_exon_file.sorted.tsv`;
    `rm $skipped_exon_file.temp`;

    open(INF, "<$skipped_exon_file.sorted.tsv") or die "couldn't open sorted skipped exon file";
    my @prev_amp = ();
    my $prev_line_min_amp;
    my $previous_strand;
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        my $line_minus_amp = join("\t", @split_line[0..3]);
        my $amp = $split_line[4];
        my $strand = $split_line[5];
        if ($. == 1) {
            push(@prev_amp, $amp);
            $prev_line_min_amp = $line_minus_amp;
            $previous_strand = $strand;
        }
        elsif ($line_minus_amp eq $prev_line_min_amp) {
            push(@prev_amp, $amp);
        }
        else {
            push(@skipped_exon_array, $prev_line_min_amp."\t".(eval join '+', @prev_amp)/scalar (@prev_amp)."\t".$previous_strand);
            $prev_line_min_amp = $line_minus_amp;
            @prev_amp = ();
            push(@prev_amp, $amp);
            $previous_strand = $strand;
        }
    }
    close(INF);
    `rm $skipped_exon_file.sorted.tsv`;
}

sub prepare_annotation_file {
    print "Processing annotation file...\n\n"; #Sort and put into array
    `sort -k1,1 -k2,2n -k3,3n $annotation_file > $annotation_file.sorted.tsv`;

    open(INF, "<$annotation_file.sorted.tsv") or die "couldn't open sorted annotation file";
    while(my $line = <INF>) {
        chomp($line);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        push(@annotation_array, $line);
    }
    `rm $annotation_file.sorted.tsv`;
    close(INF);
}

sub link_annotation_to_skipped_exons {
    print "Linking skipped exons to respective transcripts...\n\n";
    my $k = 0;
    #my $output_file = $directory."\/SEUnannotated_".$basename.".FDR_".$FDR."\/test.tsv"; # For testing
    #open(OUT, ">$output_file") or die "couldn't open unannotated skipped exon statistics output file"; # For testing
    for (my $i=0; $i < scalar (@skipped_exon_array); $i++) {
        my @split_SE_element = split("\t", $skipped_exon_array[$i]);
        my $hit_count = 0;
        my $loop_after_hit_count = 0;
        for (my $j = $k; $j < scalar (@annotation_array); $j++) {
            my @split_ann_element = split("\t", $annotation_array[$j]);
            if ($hit_count > 0) {
                $loop_after_hit_count++;
            }
            if ($split_SE_element[0] eq $split_ann_element[0] and $split_SE_element[5] eq $split_ann_element[5] and $split_SE_element[1] > $split_ann_element[1] and $split_SE_element[2] < $split_ann_element[2]) {
                my @block_sizes = split("\,", $split_ann_element[10]);
                my @block_starts = split("\,", $split_ann_element[11]);
                my $SE_donor_count = 0;
                my $SE_acceptor_count = 0;
                for (my $z = 0; $z < scalar (@block_sizes); $z++) {
                    my $donor_length = $split_ann_element[1] + $block_starts[$z] + $block_sizes[$z];
                    if ($donor_length == $split_SE_element[1]) {
                        $SE_donor_count++;
                    }
                    my $acceptor_length = $split_ann_element[1] + $block_starts[$z];
                    if ($acceptor_length == $split_SE_element[2]) {
                        $SE_acceptor_count++;
                    }
                }
                if ($SE_donor_count == 1 and $SE_acceptor_count == 1) {
                    push(@merged_SE_ann_array, $annotation_array[$j]."\t\t".$skipped_exon_array[$i]); 
                    #print OUT $annotation_array[$j], "\t\t", $skipped_exon_array[$i], "\n"; # For testing
                }
                $hit_count++;
                $k = $j-500;
            }
            elsif ($hit_count > 0 and $loop_after_hit_count > 500) {
                last;
            }
        }
    }
    #close(OUT); # For testing
}

sub annotated_vs_unannotated_SE_analysis {
    print "Identifying annotated and unannotated skipped exons...\n";
    open(OUT, ">$outdir/1_annotated_unannotated_SE_counts.txt") or die "couldn't open skipped exon statistics output file";
    my $header = "chr\tUpExonDonor\tDoExon accept\tGeneID\tIncDiff\tStrand\n";
    open(OUTup, ">$outdir/2_unannotated_pos_IncDiff.tsv") or die "couldn't open unannotated skipped exon statistics output file";
    print OUTup $header;
    open(OUTun, ">$outdir/3_unannotated_neg_IncDiff.tsv") or die "couldn't open unannotated skipped exon statistics output file";
    print OUTun $header;
    open(OUTap, ">$outdir/4_annotated_pos_IncDiff.tsv") or die "couldn't open annotated skipped exon statistics output file";
    print OUTap $header;
    open(OUTan, ">$outdir/5_annotated_neg_IncDiff.tsv") or die "couldn't open annotated skipped exon statistics output file";
    print OUTan $header;
    my $SE_junc_match = 0;
    my $matched_sum_neg = 0;
    my $unmatched_sum_neg = 0;
    my $matched_sum_pos = 0;
    my $unmatched_sum_pos = 0;
    my $prev_SE_elem = "null";
    my $prev_amp = 0;
    my $prevSEOut;
    foreach my $elem(@merged_SE_ann_array) {
        my @split_elem = split("\t", $elem);
        my $trans_start = $split_elem[1];
        my $SE_elem = join("", @split_elem[13..16]);
        my $donor = $split_elem[14];
        my $acceptor = $split_elem[15];
        my @split_starts = split("\,", $split_elem[11]);
        my @split_sizes = split("\,", $split_elem[10]);
        if ($SE_elem ne $prev_SE_elem) {
            if ($SE_junc_match > 0) {
                if ($prev_amp < 0) {
                    $matched_sum_neg++;
                    print OUTan $prevSEOut;
                }
                elsif ($prev_amp > 0) {
                    $matched_sum_pos++;
                    print OUTap $prevSEOut;
                }
            }
            else {
                if ($prev_amp < 0) {
                    $unmatched_sum_neg++;
                    print OUTun $prevSEOut;
                }
                elsif ($prev_amp > 0) {
                    $unmatched_sum_pos++;
                    print OUTup $prevSEOut;
                }
            } 
            $prev_SE_elem = $SE_elem;
            $prev_amp = $split_elem[17];
            $prevSEOut = join("\t", @split_elem[13..18])."\n";
            $SE_junc_match = 0;
        }
        $SE_junc_match = junc_match(\@split_starts,\@split_sizes,$trans_start,$donor,$acceptor,$SE_junc_match); 
    }
    if ($SE_junc_match > 0) {
        if ($prev_amp < 0) {
            $matched_sum_neg++;
            print OUTan $prevSEOut;
        }
        elsif ($prev_amp > 0) {
            $matched_sum_pos++;
            print OUTap $prevSEOut;
        }
    }
    else {
        if ($prev_amp < 0) {
            $unmatched_sum_neg++;
            print OUTun $prevSEOut;
        }
        elsif ($prev_amp > 0) {
            $unmatched_sum_pos++;
            print OUTup $prevSEOut;
        }
    } 
    my $sum_neg = $matched_sum_neg + $unmatched_sum_neg;
    if ($sum_neg == 0) {
        $sum_neg = 1;
    }
    my $sum_pos = $matched_sum_pos + $unmatched_sum_pos;
    if ($sum_pos == 0) {
        $sum_pos = 1;
    }
    print OUT "Negative SE values\n\tAnnotated skipping events: ", $matched_sum_neg, "\n\tFraction annotated: ", $matched_sum_neg/($sum_neg), "\n\n\tUnannotated skipping events: ", $unmatched_sum_neg, "\n\tFraction unannotated: ", $unmatched_sum_neg/($sum_neg),"\n\nPositive SE values\n\tAnnotated skipping events: ", $matched_sum_pos, "\n\tFraction annotated: ", $matched_sum_pos/($sum_pos), "\n\n\tUnannotated skipping events: ", $unmatched_sum_pos, "\n\tFraction unannotated: ", $unmatched_sum_pos/($sum_pos);
    close(OUT);
    close(OUTup);
    close(OUTun);
    close(OUTap);
    close(OUTan);
}

sub junc_match {
    my $donor_hit = 0;
    for (my $i = 0; $i < scalar @{$_[0]}; $i++) {
        my $isoform_exon_donor = $_[2] + ${$_[0]}[$i] + ${$_[1]}[$i];
        my $isoform_acceptor = $_[2] + ${$_[0]}[$i];
        if ($donor_hit == 0) {
            if ($_[3] == $isoform_exon_donor) {
                $donor_hit++;
            }
        }
        elsif ($donor_hit == 1) {
            $donor_hit++;
            if ($_[4] == $isoform_acceptor) {
                $_[5]++;
            }
        }
    }
    return $_[5];
}

options;
qc;
print "\n\n***********************\nRunning SEUnannotated...\n";
makedirectories;
print_input_parameters;
prepare_skipped_exon_array;
prepare_annotation_file;
link_annotation_to_skipped_exons;
annotated_vs_unannotated_SE_analysis;

print "\nSEUnannotated Done!\n***********************\n\n";