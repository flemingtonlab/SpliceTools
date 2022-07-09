use warnings;
use strict;  
use File::Basename;

my $retained_intron_file;
my $FDR;
my $annotation_file;
my $directory;
my $basename;
my $annotation_file_basename;
my @neg_IncDiff_counts_array = ();
my @pos_IncDiff_counts_array = ();


sub program_info {
    print "\n\tRIIntronExonSizes.pl will generate:\n\t\t- lists of sizes for upstream exon, retained intron, and downstream exon for statistically\n\t\t  significant positive and negative IncDiff retained intron events\n\t\t- summary file with average sizes\n\n\tUsage: perl RIIntronExonSizes.pl [OPTIONS] -r <retained intron file (rMATS JCEC)> -a <bed12 annotation file> -f <FDR>\n\n\tRequired:\n\t\t-r <retained intron file>\n\t\t-a <bed12 annotation file>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl RIIntronExonSizes.pl -r PATH/RIfile.txt -a PATH/bed12_annotation_file -f 0.05\n\n";
    exit;
}

sub options {
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }    
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] eq "\-r") {
            $retained_intron_file = $ARGV[$i+1];
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
    if ($FDR !~ m/^\d/) {
        print "\nFDR is not numeric!\n";
        program_info;
        exit;
    } 
    elsif (not defined($retained_intron_file)) {
        print "\nRetained intron file not defined!\n";
        program_info;
        exit;
    }
    open(INF, "<$retained_intron_file") or die "couldn't open retained intron file";
    while(my $line = <INF>) {
        chomp($line);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        if ($. == 2) {
            
            if (scalar @split_line != 23 || $split_line[4] !~ m/[+-]/ || $split_line[5] !~ m/^\d/ || $split_line[6] !~ m/^\d/ || $split_line[7] !~ m/^\d/ || $split_line[8] !~ m/^\d/ || $split_line[9] !~ m/^\d/ || $split_line[10] !~ m/^\d/) {
                print "\nRetained intron file may not be an rMATS JCEC file format!\n\n";
                program_info;
                exit;
            }
            elsif ($split_line[3] !~ m/^chr/) {
                print "\nRetained intron file may not be an rMATS JCEC file format!\n     Is fourth column chromosome? Fourth column doesn't begin with \"chr\"\n\n";
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

my @annotation_array  = ();
my @annotation_counts_array = ();
sub makedirectories {
    $directory = dirname($retained_intron_file);
    $basename = basename($retained_intron_file);
    $annotation_file_basename = basename($annotation_file);
    `mkdir $directory"\/RIIntronExonSizes_"$basename"_FDR_"$FDR`;
}

sub print_input_parameters {
    my $output_param_file = $directory."\/RIIntronExonSizes_".$basename."_FDR_".$FDR."\/1_1_run_parameters.txt";
    open(OUTparam, ">$output_param_file") or die "couldn't open run parameters output file";
    print OUTparam "Retained intron file: ", $basename, "\nAnnotation file: ", $annotation_file_basename, "\nFDR: ", $FDR;
    print "\nINPUT:\nRetained intron file: ", $basename, "\nAnnotation file: ", $annotation_file_basename, "\nFDR: ", $FDR, "\n";
    close(OUTparam);
}

sub average_of_array {
    my $array_length = @_;
    if ($array_length == 0) {
        return "no significant RIs";
    }
    else {
        my $total = 0;
        for my $element(@_) {
            $total += $element;
        }
        return int($total/$array_length);
    }
}

sub RI_exon_intron_sizes {
    print "\nCalculating exon and intron sizes for neg IncDiff and for pos IncDiff RI events...\n\n";
    my @upstream_exon_size_array_neg_IncDiff = ();
    my @intron_size_array_neg_IncDiff = ();
    my @downstream_exon_size_array_neg_IncDiff = ();
    my @upstream_exon_size_array_pos_IncDiff = ();
    my @intron_size_array_pos_IncDiff = ();
    my @downstream_exon_size_array_pos_IncDiff = ();

    open(INF, "<$retained_intron_file") or die "couldn't open retained intron input file";
    while(my $line = <INF>) {
        chomp($line);
        next if ($. == 1);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        if ($split_line[19] < $FDR)  {
            my $strand = $split_line[4];
            my $amplitude = $split_line[22];
            if ($strand eq "+") {
                if ($amplitude < 0) {
                    push(@upstream_exon_size_array_neg_IncDiff, $split_line[8]-$split_line[7]+1);
                    push(@intron_size_array_neg_IncDiff, $split_line[9]-$split_line[8]-1);
                    push(@downstream_exon_size_array_neg_IncDiff, $split_line[10]-$split_line[9]+1);
                }
                elsif ($amplitude > 0) {
                    push(@upstream_exon_size_array_pos_IncDiff, $split_line[8]-$split_line[7]+1);
                    push(@intron_size_array_pos_IncDiff, $split_line[9]-$split_line[8]-1);
                    push(@downstream_exon_size_array_pos_IncDiff, $split_line[10]-$split_line[9]+1);
                }
            }
            elsif ($strand eq "-") {
                if ($amplitude < 0) {
                    push(@upstream_exon_size_array_neg_IncDiff, $split_line[10]-$split_line[9]+1);
                    push(@intron_size_array_neg_IncDiff, $split_line[9]-$split_line[8]-1);
                    push(@downstream_exon_size_array_neg_IncDiff, $split_line[8]-$split_line[7]+1);
                }
                elsif ($amplitude > 0) {
                    push(@upstream_exon_size_array_pos_IncDiff, $split_line[10]-$split_line[9]+1);
                    push(@intron_size_array_pos_IncDiff, $split_line[9]-$split_line[8]-1);
                    push(@downstream_exon_size_array_pos_IncDiff, $split_line[8]-$split_line[7]+1);
                }
            }
        }
    }
    close(INF);


    my $output_file_neg_IncDiff = $directory."\/RIIntronExonSizes_".$basename."_FDR_".$FDR."\/2_intron_exon_sizes_neg_IncDiff.tsv";
    open(OUTneg, ">$output_file_neg_IncDiff") or die "couldn't open output file";
    print OUTneg "Upstream exon\tRetained intron\tDownstream exon\n";
    my $neg_IncDiff_array_length = @upstream_exon_size_array_neg_IncDiff;
    for (my $i = 0; $i < $neg_IncDiff_array_length; $i++) {
        print OUTneg $upstream_exon_size_array_neg_IncDiff[$i], "\t", $intron_size_array_neg_IncDiff[$i], "\t", $downstream_exon_size_array_neg_IncDiff[$i], "\n";
    }
    close(OUTneg);


    my $output_file_pos_IncDiff = $directory."\/RIIntronExonSizes_".$basename."_FDR_".$FDR."\/3_intron_exon_sizes_pos_IncDiff.tsv";
    open(OUTpos, ">$output_file_pos_IncDiff") or die "couldn't open output file";
    print OUTpos "Upstream exon\tRetained intron\tDownstream exon\n";
    my $pos_IncDiff_array_length = @upstream_exon_size_array_pos_IncDiff;
    for (my $i = 0; $i < $pos_IncDiff_array_length; $i++) {
        print OUTpos $upstream_exon_size_array_pos_IncDiff[$i], "\t", $intron_size_array_pos_IncDiff[$i], "\t", $downstream_exon_size_array_pos_IncDiff[$i], "\n";
    }
    close(OUTpos);


    my $output_file = $directory."\/RIIntronExonSizes_".$basename."_FDR_".$FDR."\/1_intron_exon_size_SUMMARY.tsv";
    open(OUTstats, ">$output_file") or die "couldn't open output file";
    print OUTstats "Retained intron file: ", $basename, "\tFDR = ", $FDR, "\n\tAve upstream exon\tAve retained intron\tAve downstream exon\n";

    my @averages_pos_IncDiff = ();
    push(@averages_pos_IncDiff, average_of_array(@upstream_exon_size_array_pos_IncDiff));
    push(@averages_pos_IncDiff, average_of_array(@intron_size_array_pos_IncDiff));
    push(@averages_pos_IncDiff, average_of_array(@downstream_exon_size_array_pos_IncDiff));
    print OUTstats "Positive IncDiff\t", join("\t", @averages_pos_IncDiff), "\n";

    my @averages_neg_IncDiff = ();
    push(@averages_neg_IncDiff, average_of_array(@upstream_exon_size_array_neg_IncDiff));
    push(@averages_neg_IncDiff, average_of_array(@intron_size_array_neg_IncDiff));
    push(@averages_neg_IncDiff, average_of_array(@downstream_exon_size_array_neg_IncDiff));
    print OUTstats "Negative IncDiff\t", join("\t", @averages_neg_IncDiff), "\n\n";
}

sub prepare_annotation_file {
    print "Processing annotation file...\n\n"; #Sorting and putting into array
    `sort -k1,1 -k2,2n -k3,3n $annotation_file > $annotation_file.sorted.tsv`;

    open(INF, "<$annotation_file.sorted.tsv") or die "couldn't open sorted annotation file";
    my $count = 0;
    while(my $line = <INF>) {
        if ($count == 10)  {
            chomp($line);
            $line =~ s/\"//g;
            push(@annotation_array, $line);
            $count = 0;
        }
        else {
            $count++;
        }
    }
    `rm $annotation_file.sorted.tsv`;
    close(INF);
}

sub annotated_intron_exon_sizes {
    print "Calculating annotated intron/exon sizes at potential retained intron configurations...\n\n";
    my @upstream_exon_size_array = ();
    my @intron_size_array = ();
    my @downstream_exon_size_array = ();
    
    my $annotation_array_size = @annotation_array;
    for (my $j = 0; $j < $annotation_array_size; $j++) {
        my $annotation_array_line = $annotation_array[$j];
        my @split_annotation_array_line = split("\t", $annotation_array_line);
        my @block_sizes = split("\,", $split_annotation_array_line[10]);#
        my @block_starts = split("\,", $split_annotation_array_line[11]);#
        my $block_number = $split_annotation_array_line[9];#
        my $strand = $split_annotation_array_line[5];
        if ($strand eq "+") {
            if ($block_number >1) {
                for (my $z = 1; $z < $block_number; $z++) {
                    push(@upstream_exon_size_array, $block_sizes[$z-1]);
                    push(@intron_size_array, $block_starts[$z] - ($block_starts[$z-1] + $block_sizes[$z-1]));
                    push(@downstream_exon_size_array, $block_sizes[$z]);
                }
            }
        }
        elsif ($strand eq "-") {
            if ($block_number >1) {
                for (my $z = 2; $z < $block_number; $z++) {
                    push(@downstream_exon_size_array, $block_sizes[$z-1]);
                    push(@intron_size_array, $block_starts[$z] - ($block_starts[$z-1] + $block_sizes[$z-1]));
                    push(@upstream_exon_size_array, $block_sizes[$z]);
                }
            }
        }
    }
    
    my $output_file_annotation = $directory."\/RIIntronExonSizes_".$basename."_FDR_".$FDR."\/4_intron_exon_sizes_annotated_exons.tsv";
    open(OUTannot, ">$output_file_annotation") or die "couldn't open output file";
    print OUTannot "Upstream exon\tRetained intron\tDownstream exon\n";
    my $annotated_array_length = @upstream_exon_size_array;
    for (my $i = 0; $i < $annotated_array_length; $i++) {
        print OUTannot $upstream_exon_size_array[$i], "\t", $intron_size_array[$i], "\t", $downstream_exon_size_array[$i], "\n";
    }
    close(OUTannot);


    my @averages = ();
    push(@averages, average_of_array(@upstream_exon_size_array));
    push(@averages, average_of_array(@intron_size_array));
    push(@averages, average_of_array(@downstream_exon_size_array));
    print OUTstats "Annotation data\t", join("\t", @averages), "\n";
    
}
close(OUTstats);

options;
qc;
print "\n\n***********************\nRunning RIIntronExonSizes...\n";
makedirectories;
print_input_parameters;
RI_exon_intron_sizes;
prepare_annotation_file;
annotated_intron_exon_sizes;
print "RIIntronExonSizes DONE!\n***********************\n\n";