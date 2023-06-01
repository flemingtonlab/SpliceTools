use warnings;
use strict;  
use File::Basename;

my $skipped_exon_file;
my $FDR;
my $annotation_file;

my @skipped_exon_array = ();
my @annotation_array = ();
my $sorted_merged_output_file;
my $outdir;
my $basename;
my @neg_IncDiff_counts_array = ();
my @pos_IncDiff_counts_array = ();

sub program_info {
    print "\n\tSENumberSkipped.pl will generate:\n\t\t- lists of all predicted SE isoforms with number of exons skipped\n\t\t- lists of isoforms with maximum number of skipped exons\n\t\t- statistics for number of exons skipped\n\n\tNote: Uses annotated exon information to determine number of exons skipped\n\t(0 exons skipped means that skipped exons were not in annotation file)\n\n\n\tUsage: perl SENumberSkipped.pl [OPTIONS] -s <skipped exon file (rMATS JCEC)> -a <bed12 annotation file> -f <FDR>\n\n\tRequired:\n\t\t-s <skipped exon file>\n\t\t-a <bed12 annotation file>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl SENumberSkipped.pl -s PATH/SEfile.txt -a PATH/bed12_annotation.bed -f 0.05\n\n";
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

sub makedirectories {
    $basename = basename($skipped_exon_file);
    $basename =~ s/\.txt$//g;
    $outdir = dirname($skipped_exon_file)."\/SENumberSkipped_".$basename."_FDR_".$FDR;
    `mkdir $outdir`;
    `mkdir $outdir"\/1_neg_IncDiff"`;
    `mkdir $outdir"\/2_pos_IncDiff"`;
}

sub print_output_parameters {
    open(OUTparam, ">$outdir\/1_1_input_parameters.txt") or die "couldn't open run parameters output file";

    print OUTparam "Skipped exon file: ", basename($skipped_exon_file), "\n", "Annotation file: ", basename($annotation_file), "\n", "FDR: ", $FDR, "\n\n";
    print "\nINPUT:\nSkipped exon file: ", basename($skipped_exon_file), "\n", "Annotation file: ", basename($annotation_file), "\n", "FDR: ", $FDR, "\n\n";
    close(OUTparam);
}

sub prepare_skipped_exon_array {
    print "Preparing skipped exon data and putting into array...\n\n";
    open(INF, "<$skipped_exon_file") or die "couldn't open skipped exon input file";
    open(OUT, ">$skipped_exon_file.temp") or die "couldn't open temporary skipped exon output file";
    while(my $line = <INF>) {
        chomp($line);
        next if ($. == 1);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        if ($split_line[19] <= $FDR) {
            print OUT $split_line[3], "\t", $split_line[8], "\t", $split_line[9], "\t", $split_line[2], "\t", $split_line[22], "\n";
        }
    }
    close(INF);
    close(OUT);

    `sort -k1,1 -k2,2n -k3,3n $skipped_exon_file.temp > $skipped_exon_file.sorted.tsv`;
    `rm $skipped_exon_file.temp`;

    open(INF, "<$skipped_exon_file.sorted.tsv") or die "couldn't open sorted output skipped exon file";
    open(OUT, ">$skipped_exon_file.sorted.deduplicated.tsv") or die "couldn't open deduplicated skipped exon output file";
    my @previous_amplitude_array = ();
    my $previous_line_minus_amplitude;
    my $line_minus_amplitude;
    my $amplitude;
    my $count = 1;
    my $average_amplitude;
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        $line_minus_amplitude = join("\t", @split_line[0..3]);
        $amplitude = $split_line[4];

        if ($. == 1) {
            push(@previous_amplitude_array, $amplitude);
            $previous_line_minus_amplitude = $line_minus_amplitude;
        }
        elsif ($line_minus_amplitude eq $previous_line_minus_amplitude) {
            push(@previous_amplitude_array, $amplitude);
            $count++;
        }
        else {
            my $total = eval join '+', @previous_amplitude_array;
            $average_amplitude = $total/$count;
            print OUT $previous_line_minus_amplitude, "\t", $average_amplitude, "\n";
            $previous_line_minus_amplitude = $line_minus_amplitude;
            @previous_amplitude_array = ();
            push(@previous_amplitude_array, $amplitude);
            $count = 1;
        }
    }
    my $total = eval join '+', @previous_amplitude_array;
    $average_amplitude = $total/$count;
    print OUT $previous_line_minus_amplitude, "\t", $average_amplitude, "\n";
    close(INF);
    close(OUT);
    `rm $skipped_exon_file.sorted.tsv`;

    open(INF, "<$skipped_exon_file.sorted.deduplicated.tsv") or die "couldn't open sorted output skipped exon file";
    while(my $line = <INF>) {
        chomp($line);
        push(@skipped_exon_array, $line);
    }
    close(INF);
    `rm $skipped_exon_file.sorted.deduplicated.tsv`;
}

sub prepare_annotation_file {
    print "Processing annotation file...\n\n"; #Sorting and putting into array
    `sort -k1,1 -k2,2n -k3,3n $annotation_file > $annotation_file.sorted.tsv`;

    open(INF, "<$annotation_file.sorted.tsv") or die "couldn't open sorted annotation file";
    my $count = 0;
    while(my $line = <INF>) {
        chomp($line);
        $line =~ s/\"//g;
        push(@annotation_array, $line);
    }
    `rm $annotation_file.sorted.tsv`;
    close(INF);
}

sub combine_annotation_to_skipped_exons {
    print "Linking skipped exons to annotated transcripts...\n\n";
    my $skipped_exon_array_size = @skipped_exon_array;
    my $annotation_array_size = @annotation_array;
    my $k = 0;
    open(OUT, ">$skipped_exon_file.merged.tsv") or die "couldn't open merged skipped exon-annotation output file";
    for (my $i=0; $i < $skipped_exon_array_size; $i++) {
        my $skipped_exon_line = $skipped_exon_array[$i];
        my @split_skipped_exon_line = split("\t", $skipped_exon_line);
        my $hit_count = 0;
        my $loop_after_hit_count = 0;
        for (my $j = $k; $j < $annotation_array_size; $j++) {
            my $annotation_array_line = $annotation_array[$j];
            my @split_annotation_array_line = split("\t", $annotation_array_line);
            if ($hit_count > 0) {
                $loop_after_hit_count++;
            }
            if ($split_skipped_exon_line[0] eq $split_annotation_array_line[0] and $split_skipped_exon_line[1] > $split_annotation_array_line[1] and $split_skipped_exon_line[2] < $split_annotation_array_line[2]) {
                my @block_sizes = split("\,", $split_annotation_array_line[10]);#
                my @block_starts = split("\,", $split_annotation_array_line[11]);#
                my $block_sizes_length = @block_sizes;#
                my $SE_donor_count = 0;
                my $SE_acceptor_count = 0;
                for (my $z = 0; $z < $block_sizes_length; $z++) {
                    my $donor_length = $split_annotation_array_line[1] + $block_starts[$z] + $block_sizes[$z];
                    if ($donor_length == $split_skipped_exon_line[1]) {
                        $SE_donor_count++;
                    }
                    my $acceptor_length = $split_annotation_array_line[1] + $block_starts[$z];
                    if ($acceptor_length == $split_skipped_exon_line[2]) {
                        $SE_acceptor_count++;
                    }
                }
                if ($SE_donor_count == 1 and $SE_acceptor_count == 1) {
                    print OUT $annotation_array_line, "\t\t", $skipped_exon_line, "\n";
                    
                }
                $hit_count++;
                $k = $j-1000;
            }
            elsif ($hit_count > 0 and $loop_after_hit_count > 500) {
                last;
            }
        }
    }
    close(OUT);
}

sub sort_merged_file {
    print "Sorting merged SE and annotation file...\n\n";
    my $merged_input_file = $skipped_exon_file.".merged.tsv";
    $sorted_merged_output_file = $skipped_exon_file.".merged.sorted.tsv";
    `sort -k14,14 -k15,15n -k16,16n $merged_input_file > $sorted_merged_output_file`;
    `rm $merged_input_file`;  
}

sub count_max_skipped_exons {
    print "Counting max number of skipped exons for each SE event...\n\n";
    my $all_isoform_skipping_counts_output_file_neg_IncDiff = $outdir."\/1_neg_IncDiff\/3_all_isoform_skipped_counts_neg_IncDiff.tsv";
    my $all_isoform_skipping_counts_output_file_pos_IncDiff = $outdir."\/2_pos_IncDiff\/3_all_isoform_skipped_counts_pos_IncDiff.tsv";
    my $max_skipping_counts_isoform_output_file_neg_IncDiff = $outdir."\/1_neg_IncDiff\/2_isoform_w_max_number_skipped_exons_neg_IncDiff.tsv";
    my $max_skipping_counts_isoform_output_file_pos_IncDiff = $outdir."\/2_pos_IncDiff\/2_isoform_w_max_number_skipped_exons_pos_IncDiff.tsv";
    
    open(INF, "<$sorted_merged_output_file") or die "couldn't open sorted merged output file file";
    open(OUT1, ">$all_isoform_skipping_counts_output_file_neg_IncDiff") or die "couldn't open output file";
    open(OUT2, ">$all_isoform_skipping_counts_output_file_pos_IncDiff") or die "couldn't open output file";
    open(OUT3, ">$max_skipping_counts_isoform_output_file_neg_IncDiff") or die "couldn't open output file";
    open(OUT4, ">$max_skipping_counts_isoform_output_file_pos_IncDiff") or die "couldn't open output file";

    print OUT1 "chr\tisoform start\tisoform end\tisoform_gene\tintensity\tstrand\tcoding start\tcoding end\tcolor\tblock number\tblock sizes\tblock starts\t\tSE chr\tSkipping start\tSkipping end\tgene\tIncDiff\tNumber of exons skipped\n";
    print OUT2 "chr\tisoform start\tisoform end\tisoform_gene\tintensity\tstrand\tcoding start\tcoding end\tcolor\tblock number\tblock sizes\tblock starts\t\tSE chr\tSkipping start\tSkipping end\tgene\tIncDiff\tNumber of exons skipped\n";
    print OUT3 "chr\tisoform start\tisoform end\tisoform_gene\tintensity\tstrand\tcoding start\tcoding end\tcolor\tblock number\tblock sizes\tblock starts\t\tSE chr\tSkipping start\tSkipping end\tgene\tIncDiff\tNumber of exons skipped\n";
    print OUT4 "chr\tisoform start\tisoform end\tisoform_gene\tintensity\tstrand\tcoding start\tcoding end\tcolor\tblock number\tblock sizes\tblock starts\t\tSE chr\tSkipping start\tSkipping end\tgene\tIncDiff\tNumber of exons skipped\n";
    
    my $previous_SE_ID = "null";
    my $max_counts_line;
    my $running_skipped_exons_count = 0;
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        my $chr = $split_line[0];
        my $transcript_start = $split_line[1];
        my @split_exon_starts = split("\,", $split_line[11]);
        my @split_exon_sizes = split("\,", $split_line[10]);
        my $exon_number = $split_line[9];

        my $SE_ID = join("\t", @split_line[13..16]);
        my $SE_donor = $split_line[14];
        my $SE_acceptor = $split_line[15];
        my $IncDiff = $split_line[17];

        my $donor_hit = 0;
        my $acceptor_hit = 0;
        my $skipped_exon_count = 0;
        for (my $i = 0; $i < $exon_number; $i++) {
            my $block_donor = $transcript_start + $split_exon_starts[$i] + $split_exon_sizes[$i];
            my $block_acceptor = $transcript_start + $split_exon_starts[$i];
            if ($block_donor == $SE_donor) {
                $donor_hit++;
            }
            elsif($block_acceptor == $SE_acceptor) {
                $acceptor_hit++;
            }
            elsif ($donor_hit == 1 and $acceptor_hit == 0) {
                $skipped_exon_count++;
            }
        }
        if ($. == 1) {
            $running_skipped_exons_count = $skipped_exon_count;
            $max_counts_line = $line;
            $previous_SE_ID = $SE_ID;
        }
        elsif ($SE_ID eq $previous_SE_ID) {
            if ($skipped_exon_count > $running_skipped_exons_count) {
                $running_skipped_exons_count = $skipped_exon_count;
                $max_counts_line = $line;
            }
            if ($IncDiff < 0) {
                print OUT1 $line, "\t", $skipped_exon_count, "\n";
            }
            elsif ($IncDiff > 0) {
                print OUT2 $line, "\t", $skipped_exon_count, "\n";
            }
        }
        elsif ($SE_ID ne $previous_SE_ID) {  
            my @split_max_counts_line = split("\t", $max_counts_line);
            if ($split_max_counts_line[17] < 0) {
                print OUT1 $max_counts_line, "\t", $running_skipped_exons_count, "\n";
                print OUT3 $max_counts_line, "\t", $running_skipped_exons_count, "\n";
            }
            elsif ($split_max_counts_line[17] > 0) {
                print OUT2 $max_counts_line, "\t", $running_skipped_exons_count, "\n";
                print OUT4 $max_counts_line, "\t", $running_skipped_exons_count, "\n";
            }
            if ($IncDiff < 0) {
                push (@neg_IncDiff_counts_array, $running_skipped_exons_count);
            }
            elsif ($IncDiff > 0) {
                push (@pos_IncDiff_counts_array, $running_skipped_exons_count);
            }
            $running_skipped_exons_count = $skipped_exon_count;
            $max_counts_line = $line;
            $previous_SE_ID = $SE_ID;
        }
    }
    close(INF);
    close(OUT1);
    close(OUT2);
    close(OUT3);
    close(OUT4);
}

sub neg_IncDiff_summary {
    print "Summary data for negative IncDiff events...\n\n";
    #my $a;
    #my $b;
    my @sorted_neg_IncDiff_counts_array = sort { $a <=> $b } @neg_IncDiff_counts_array;
    my $sorted_neg_IncDiff_counts_array_size = @sorted_neg_IncDiff_counts_array;

    my $statistics_output_file = $outdir."\/1_neg_IncDiff\/1_statistics_neg_IncDiff.tsv";
    open(OUT, ">$statistics_output_file") or die "couldn't open output file";
    print OUT "Skipped Exon #\tNumber of events\n";
    my $number_of_skipped_exons = 0;
    my $number_of_skipping_events_per_count = 0;
    for (my $a = 0; $a < $sorted_neg_IncDiff_counts_array_size; $a++) {
        if ($a == $sorted_neg_IncDiff_counts_array_size - 1) {
            if ($number_of_skipped_exons == $sorted_neg_IncDiff_counts_array[$a]) {
                $number_of_skipping_events_per_count++;
                print OUT $number_of_skipped_exons, "\t", $number_of_skipping_events_per_count, "\n";
            }
            elsif ($number_of_skipped_exons < $sorted_neg_IncDiff_counts_array[$a]) {
                print OUT $number_of_skipped_exons, "\t", $number_of_skipping_events_per_count, "\n";
                $number_of_skipped_exons = $sorted_neg_IncDiff_counts_array[$a];
                $number_of_skipping_events_per_count = 1;
                print OUT $number_of_skipped_exons, "\t", $number_of_skipping_events_per_count;
            }
        }
        elsif ($a == 0 || $number_of_skipped_exons == $sorted_neg_IncDiff_counts_array[$a]) {
            $number_of_skipped_exons = $sorted_neg_IncDiff_counts_array[$a];
            $number_of_skipping_events_per_count++;
        }
        elsif ($number_of_skipped_exons < $sorted_neg_IncDiff_counts_array[$a]) {
            print OUT $number_of_skipped_exons, "\t", $number_of_skipping_events_per_count, "\n";
            $number_of_skipped_exons = $sorted_neg_IncDiff_counts_array[$a];
            $number_of_skipping_events_per_count = 1;
        }
    }
    close(OUT);
}

sub pos_IncDiff_summary {
    print "Summary data for positive IncDiff events...\n\n";
    #my $a;
    #my $b;
    my @sorted_pos_IncDiff_counts_array = sort { $a <=> $b } @pos_IncDiff_counts_array;
    my $sorted_pos_IncDiff_counts_array_size = @sorted_pos_IncDiff_counts_array;

    my $statistics_output_file = $outdir."\/2_pos_IncDiff\/1_statistics_pos_IncDiff.tsv";
    open(OUT, ">$statistics_output_file") or die "couldn't open output file";
    print OUT "Skipped Exon #\tNumber of events\n";
    my $number_of_skipped_exons = 0;
    my $number_of_skipping_events_per_count = 0;
    for (my $a = 0; $a < $sorted_pos_IncDiff_counts_array_size; $a++) {
        if ($a == $sorted_pos_IncDiff_counts_array_size - 1) {
            if ($number_of_skipped_exons == $sorted_pos_IncDiff_counts_array[$a]) {
                $number_of_skipping_events_per_count++;
                print OUT $number_of_skipped_exons, "\t", $number_of_skipping_events_per_count, "\n";
            }
            elsif ($number_of_skipped_exons < $sorted_pos_IncDiff_counts_array[$a]) {
                print OUT $number_of_skipped_exons, "\t", $number_of_skipping_events_per_count, "\n";
                $number_of_skipped_exons = $sorted_pos_IncDiff_counts_array[$a];
                $number_of_skipping_events_per_count = 1;
                print OUT $number_of_skipped_exons, "\t", $number_of_skipping_events_per_count;
            }
        }
        elsif ($a == 0 || $number_of_skipped_exons == $sorted_pos_IncDiff_counts_array[$a]) {
            $number_of_skipped_exons = $sorted_pos_IncDiff_counts_array[$a];
            $number_of_skipping_events_per_count++;
        }
        elsif ($number_of_skipped_exons < $sorted_pos_IncDiff_counts_array[$a]) {
            print OUT $number_of_skipped_exons, "\t", $number_of_skipping_events_per_count, "\n";
            $number_of_skipped_exons = $sorted_pos_IncDiff_counts_array[$a];
            $number_of_skipping_events_per_count = 1;
        }
    }
    close(OUT);
}


options;
qc;
print "\n\n***********************\nRunning SENumberSkipped...\n";
makedirectories;
print_output_parameters;
prepare_skipped_exon_array;
prepare_annotation_file;
combine_annotation_to_skipped_exons;
sort_merged_file;
count_max_skipped_exons;
neg_IncDiff_summary;
pos_IncDiff_summary;

`rm $skipped_exon_file.merged.sorted.tsv`;

print "SENumberSkipped DONE!\n***********************\n\n";