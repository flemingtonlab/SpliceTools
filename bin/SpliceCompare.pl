use warnings;
use strict;  
use File::Basename; 

my $path_to_files;
my $output_directory;
my $overlap_out = 1;
my $max_minus_log_10_hypergeo_val = "280";
my $FDR;

my $dirname = dirname(__FILE__);

my @files = ();
my $out_dir;

my @unique_sample_name_array = ();
my @splice_type_array = ();

sub program_info {
    print "\n\tSpliceCompare.pl takes groups of rMATS differential splicing files and outputs lists of\n\tstatistically significant events with positive IncDiff and negative IncDiff and summary files.\n\tIt then compares differential splicing changes across experiments to assess functional\n\trelationships (hypergeometric test), outputting pval matrices, -log10(pval) matrices,\n\tand .svg graphical representation files displaying -log10(pval) of all comparisons.\n\n\tImportant notes:\n\t\tInput file names must end in one or more of the following suffexes:\n\t\t\tA3SS.MATS.JCEC.txt\n\t\t\tA5SS.MATS.JCEC.txt\n\t\t\tMXE.MATS.JCEC.txt\n\t\t\tRI.MATS.JCEC.txt\n\t\t\tSE.MATS.JCEC.txt\n\n\t\tBy default, SpliceCompare.pl generates a file of common splicing changes for each comparison.\n\t\tFor large numbers of input files, this can lead to the generation of millions\n\t\tof comparison files. Use option -p 0 to suppress this output.\n\n\t\tBy default, -log10(pval) maximum is 280. For experiments with less significant\n\t\thypergeometric test values, max -log10(pval) can be decreased using -m option.\n\n\tUsage: perl SpliceCompare.pl [OPTIONS] -i <input files directory (rMATS JCEC)> -o <path to output directory> -f <FDR>\n\n\tRequired:\n\t\t-i <input files directory (rMATS JCEC)>\n\t\t-o <path to output directory>\n\t\t-f <FDR>\n\n\tOptional:\n\t\t-p <1=output overlap files (default), 0=no overlap files>\n\t\t-m <max cluster array value (default = 280)>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl SpliceCompare.pl -i PATH/input_files_directory -o PATH/output_files_directory -m 280 -p 1 -f 0.05\n\n";
    exit;
}

sub options {
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }  
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] eq "\-i") {
            $path_to_files = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-o") {
            $output_directory = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-p") {
            $overlap_out = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-m") {
            $max_minus_log_10_hypergeo_val = $ARGV[$i+1];
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
    if ($path_to_files eq "") {
        print "\nInput files directory not defined!\n";
        program_info;
        exit;
    }
    elsif ($output_directory eq "") {
        print "\nOutput directory not defined!\n";
        program_info;
        exit;
    }
    elsif ($FDR eq "") {
        print "\nFDR not defined!\n";
        program_info;
        exit;
    }
    elsif ($FDR !~ m/^\d/) {
        print "\nFDR is not numeric!\n";
        program_info;
        exit;
    }
    elsif ($max_minus_log_10_hypergeo_val !~ m/^\d/) {
        print "\nMax cluster map value is not numeric!\n";
        program_info;
        exit;
    }
}

sub read_input_dir {
    print "\nReading input directory...\n\n";
    opendir my $input_dir, "$path_to_files" or die "Can't open directory: $!";
    foreach my $g (sort readdir $input_dir) {
        next if ($g eq '.' || $g eq '..' || $g eq '.DS_Store');
        my $path_files = $path_to_files."/".$g;
        push(@files, $path_files);
    }
    closedir $input_dir;
}

sub make_output_dirs {
    print "Making output directories...\n\n";
    $out_dir = $output_directory."/SpliceCompare_out_FDR_".$FDR;
    `mkdir $out_dir`;
    `mkdir $out_dir/A3SS_files`;
    `mkdir $out_dir/A5SS_files`;
    `mkdir $out_dir/MXE_files`;
    `mkdir $out_dir/RI_files`;
    `mkdir $out_dir/SE_files`;

    my @subdir_array = ('A3SS_files', 'A5SS_files', 'MXE_files', 'RI_files', 'SE_files');

    foreach my $sub_dir(@subdir_array) {
        `mkdir $out_dir/$sub_dir/1_rMATS_sig`;
        `mkdir $out_dir/$sub_dir/1_rMATS_sig/1_original`;
        `mkdir $out_dir/$sub_dir/1_rMATS_sig/1_sig_all`;
        `mkdir $out_dir/$sub_dir/1_rMATS_sig/2_sig_pos_IncDiff`;
        `mkdir $out_dir/$sub_dir/1_rMATS_sig/3_sig_neg_IncDiff`;

        `mkdir $out_dir/$sub_dir/2_common_alt_splicing`;
        `mkdir $out_dir/$sub_dir/2_common_alt_splicing/2_data_matrix_pval`;
        `mkdir $out_dir/$sub_dir/2_common_alt_splicing/3_data_matrix_minus_log10_pval`;
        `mkdir $out_dir/$sub_dir/2_common_alt_splicing/4_clustermaps`;
        if ($overlap_out == 1) {
            `mkdir $out_dir/$sub_dir/2_common_alt_splicing/5_pairwise_overlap`;
            `mkdir $out_dir/$sub_dir/2_common_alt_splicing/5_pairwise_overlap/1_pos_IncDiff`;
            `mkdir $out_dir/$sub_dir/2_common_alt_splicing/5_pairwise_overlap/2_neg_IncDiff`;
        }
    }
}

sub process_JCEC_files {
    print "Processing JCEC files...\n\n";
    my @sample_name_array = ();
    foreach my $file(@files) {
        my $basename = basename($file);
        my @split_basename = split("_", $basename);
        my $splice_type = pop(@split_basename);
        my $sample_name = join("_", @split_basename);
        push (@splice_type_array, $splice_type);
        push (@sample_name_array, $sample_name);
        my $unique_path;
        if($splice_type eq "A3SS.MATS.JCEC.txt") {
            $unique_path = "A3SS_files";
        }
        elsif ($splice_type eq "A5SS.MATS.JCEC.txt") {
            $unique_path = "A5SS_files";
        }
        elsif ($splice_type eq "MXE.MATS.JCEC.txt") {
            $unique_path = "MXE_files";
        }
        elsif ($splice_type eq "RI.MATS.JCEC.txt") {
            $unique_path = "RI_files";
        }
        elsif ($splice_type eq "SE.MATS.JCEC.txt") {
            $unique_path = "SE_files";
        }
        open(INF, "<$file") or die "couldn't open input file";
        open(OUT, ">$out_dir/$unique_path/1_rMATS_sig/1_original/$basename") or die "couldn't open output file";
        open(OUT2, ">$out_dir/$unique_path/1_rMATS_sig/1_sig_all/$basename") or die "couldn't open output file";
        open(OUT3, ">$out_dir/$unique_path/1_rMATS_sig/2_sig_pos_IncDiff/$basename") or die "couldn't open output file";
        open(OUT4, ">$out_dir/$unique_path/1_rMATS_sig/3_sig_neg_IncDiff/$basename") or die "couldn't open output file";

        while(my $line = <INF>) {
            chomp($line);
            my @split_line = split("\t", $line);
            if ($. == 1) {
                print OUT $line, "\n";
                print OUT2 $line, "\n";
                print OUT3 $line, "\n";
                print OUT4 $line, "\n";
            }
            elsif (($splice_type ne "MXE.MATS.JCEC.txt" and $split_line[19] <= $FDR) || ($splice_type eq "MXE.MATS.JCEC.txt" and $split_line[21] <= $FDR)) {
                print OUT $line, "\n";
                print OUT2 $line, "\n";
                my $gene_name = $split_line[2];
                $gene_name =~ s/"//g;
                if ($splice_type eq "MXE.MATS.JCEC.txt") {
                    if ($split_line[24] > 0) {
                        print OUT3 $line, "\n";
                    }
                    elsif ($split_line[24] < 0) {
                        print OUT4 $line, "\n";
                    }
                }
                else {
                    if ($split_line[22] > 0) {
                        print OUT3 $line, "\n";
                    }
                    elsif ($split_line[22] < 0) {
                        print OUT4 $line, "\n";
                    }
                }
            }
            else {
                print OUT $line, "\n";
            }
        }
        close(INF);
        close(OUT);
        close(OUT2);
        close(OUT3);
        close(OUT4);
    }
    @unique_sample_name_array = unique(@sample_name_array);
    if (!@unique_sample_name_array) {
        print "ERROR!!\nInput directory doesn't contain alt splicing files with any of the following suffexes:\n\tA3SS.MATS.JCEC.txt\n\tA5SS.MATS.JCEC.txt\n\tMXE.MATS.JCEC.txt\n\tRI.MATS.JCEC.txt or\n\tSE.MATS.JCEC.txt!\n\n****************\n\n";
        program_info;
        exit;
    }
}

sub summarize {
    print "Summarizing alt splicing\n\n"; ####################
    open(OUT, ">$out_dir/1_1_splicing_summaries_all_samples.tsv") or die "couldn't open output file";
    
    foreach my $sample_name_file(@unique_sample_name_array) {
        print OUT $sample_name_file;
        my $A3SS_sample_output_name = $sample_name_file."_A3SS.MATS.JCEC.txt";
        my $A5SS_sample_output_name = $sample_name_file."_A5SS.MATS.JCEC.txt";
        my $MXE_sample_output_name = $sample_name_file."_MXE.MATS.JCEC.txt";
        my $RI_sample_output_name = $sample_name_file."_RI.MATS.JCEC.txt";
        my $SE_sample_output_name = $sample_name_file."_SE.MATS.JCEC.txt";
        my @total_event_array = ();
        my @induced_event_array = ();
        my @decreased_event_array = ();
        my @fraction_up_array = ();
        my @fraction_down_array = ();
        my @directory_tails = ("A3SS_files/1_rMATS_sig/1_sig_all/$A3SS_sample_output_name", "A5SS_files/1_rMATS_sig/1_sig_all/$A5SS_sample_output_name", "MXE_files/1_rMATS_sig/1_sig_all/$MXE_sample_output_name", "RI_files/1_rMATS_sig/1_sig_all/$RI_sample_output_name", "SE_files/1_rMATS_sig/1_sig_all/$SE_sample_output_name");       
        foreach my $directory_tail(@directory_tails) {
            next unless (-e $out_dir."/".$directory_tail);
            my $header = basename($directory_tail);
            my @split_header = split("_", $header);
            $header = $split_header[-1];
            $header =~ s/\.MATS\.JCEC\.txt//g;
            print OUT "\t", $header;
            open(INF, "<$out_dir/$directory_tail") or die "couldn't open input file";
            my $induced = 0;
            my $decreased = 0;
            while(my $line = <INF>) {
                chomp($line);
                my @split_line = split("\t", $line);
                next if ($. == 1);
                if ($directory_tail eq "MXE_files/1_rMATS_sig/1_sig_all/$MXE_sample_output_name") {
                    if ($split_line[24] > 0) {
                        $induced++;
                    }
                    elsif ($split_line[24] < 0) {
                    $decreased++;
                    }
                }
                else {
                    if ($split_line[22] > 0) {
                        $induced++;
                    }
                    elsif ($split_line[22] < 0) {
                        $decreased++;
                    }
                }
            }
            my $total = $induced + $decreased;
            my $fraction_up = $induced/($total+0.005);
            my $fraction_down = $decreased/($total+0.005);

            push(@total_event_array, $total);
            push(@induced_event_array, $induced);
            push(@decreased_event_array, $decreased);
            push(@fraction_up_array, $fraction_up);
            push(@fraction_down_array, $fraction_down);
            close(INF);
        }
        print OUT "\n";
        print OUT "Pos IncDiff\t", join("\t", @induced_event_array), "\n";
        print OUT "Neg IncDiff\t", join("\t", @decreased_event_array), "\n";
        print OUT "All\t", join("\t", @total_event_array), "\n\n";
        print OUT "Fraction pos IncDiff\t", join("\t", @fraction_up_array), "\n";
        print OUT "Fraction neg IncDiff\t", join("\t", @fraction_down_array), "\n\n\n";
    }
    close(OUT);
}

sub compare {
    print "Matching alt splicing events...\n\n";
    my @splice_junction_directories_array = ('A3SS_files', 'A5SS_files', 'MXE_files', 'RI_files', 'SE_files'); 
    foreach my $splice_junction_directory(@splice_junction_directories_array) {
        my @output_4_data_array = ();
        opendir my $dir, "$out_dir/$splice_junction_directory/1_rMATS_sig/1_original" or die "Can't open directory: $!";
        my @files = ();
        foreach my $f (sort readdir $dir) {
            next if ($f eq '.' || $f eq '..');
            push(@files, $f);
        }
        closedir $dir;
        if (!@files) {
            `rm -r $out_dir/$splice_junction_directory`;
            next;
        }
        my @analysis_data_array = ();
        my @array_starts = ();
        my $line_count = 0;

        my @pos = ();
        my @pos_starts = ();
        my $pos_line_count = 0;

        my @neg = ();
        my @neg_starts = ();
        my $neg_line_count = 0;

        foreach my $file(@files) {
            open(INF, "<$out_dir/$splice_junction_directory/1_rMATS_sig/1_original/$file") or die "couldn't open input file";
            push(@analysis_data_array, $file);
            push(@array_starts, $line_count);
            push(@pos, $file);
            push(@pos_starts, $pos_line_count);
            push(@neg, $file);
            push(@neg_starts, $neg_line_count);
            while(my $line = <INF>) {
                chomp($line);
                next if ($. == 1);
                $line =~ s/\"//g;
                push(@analysis_data_array, $line);
                $line_count++;
                my @split_line = split("\t", $line);
                if ($splice_junction_directory eq "MXE_files") {
                    if ($split_line[21] < $FDR) {
                        if ($split_line[24] > 0) {
                            push(@pos, join("\t", @split_line[1..12]));
                            $pos_line_count++;
                        }
                        elsif ($split_line[24] < 0) {
                            push(@neg, join("\t", @split_line[1..12]));
                            $neg_line_count++;
                        }
                    }
                }
                else {
                    if ($split_line[19] < $FDR) {
                        if ($split_line[22] > 0) {
                            push(@pos, join("\t", @split_line[1..10]));
                            $pos_line_count++;
                        }
                        elsif ($split_line[22] < 0) {
                            push(@neg, join("\t", @split_line[1..10]));
                            $neg_line_count++;
                        }
                    }
                }
            }
            $line_count++;
            $pos_line_count++;
            $neg_line_count++;
            close(INF);
        }
        my $data_file_number = @array_starts;
        push(@array_starts, $line_count+1);
        push(@analysis_data_array, "null");
        push(@pos_starts, $pos_line_count+1);
        push(@pos, "null");
        push(@neg_starts, $neg_line_count+1);
        push(@neg, "null");
        open(OUT1a, ">$out_dir/$splice_junction_directory/2_common_alt_splicing/1_2_overlap_summary.tsv") or die "couldn't open output file";

        print OUT1a "Summary statistics on splice junction overlaps\n\n";
        print OUT1a "Sample\ttotal splicing events across all samples\ttotal alt splicing events\t\tall significant alt splicing events\tall common significant alt splicing events\tall - fraction common alt splicing events\thypergeometric all\t\tpos IncDiff significant alt splicing events\tcommon pos IncDiff significant alt splicing events\tpos IncDiff - fraction common alt splicing events\thypergeometric pos IncDiff\t\tneg IncDiff significant alt splicing events\tcommon neg IncDiff significant alt splicing events\tneg IncDiff - fraction common alt splicing events\thypergeometric neg IncDiff\n";

        my $output_dir_path_2 = $out_dir."/".$splice_junction_directory."/2_common_alt_splicing";
        open(OUT1b, ">$output_dir_path_2/3_data_matrix_minus_log10_pval/1_data_matrix_all.tsv") or die "couldn't open output file";
        open(OUT1c, ">$output_dir_path_2/3_data_matrix_minus_log10_pval/2_data_matrix_pos_IncDiff.tsv") or die "couldn't open output file";
        open(OUT1d, ">$output_dir_path_2/3_data_matrix_minus_log10_pval/3_data_matrix_neg_IncDiff.tsv") or die "couldn't open output file";
        open(OUT1e, ">$output_dir_path_2/2_data_matrix_pval/1_data_matrix_pvalue_all.tsv") or die "couldn't open output file";
        open(OUT1f, ">$output_dir_path_2/2_data_matrix_pval/2_data_matrix_pvalue_pos_IncDiff.tsv") or die "couldn't open output file";
        open(OUT1g, ">$output_dir_path_2/2_data_matrix_pval/3_data_matrix_pvalue_neg_IncDiff.tsv") or die "couldn't open output file";
        my @hypergeo_output_matrix_name_array = ();
        my @hypergeo_values_all_array = ();
        my @hypergeo_vals_pos_array = ();
        my @hypergeo_vals_neg_array = ();

        open(OUT1, ">$output_dir_path_2/1_1_summary_statistics_each_sample.tsv") or die "couldn't open output file";
        print OUT1 "Sample name\tAll significant changes\tPos IncDiff signif changes\tFraction pos vs all signif IncDiff\tNeg IncDiff signif changes\tFraction neg vs all signif IncDiff\tTotal splice junctions detected\n";
        for (my $a=0; $a<$data_file_number-1; $a++) {
            for (my $b=$a+1; $b<$data_file_number; $b++) {
                next if ($a == $b);
                #print "\nSTART PAIRING...\n\n";
                my @array_1 = @analysis_data_array[$array_starts[$a]..($array_starts[$a+1]-1)];
                my @array_2 = @analysis_data_array[$array_starts[$b]..($array_starts[$b+1]-1)];
                my $array_1_name = shift(@array_1);
                my $array_2_name = shift(@array_2);
                my @combined_array_1_and_2 = (@array_1, @array_2);
                my @unique_array_1_and_2_elements = unique(@combined_array_1_and_2);
                my $total_jncts_detected_across_samples = @unique_array_1_and_2_elements;

                my @pos_array_1 = @pos[$pos_starts[$a]..($pos_starts[$a+1]-1)];
                my @pos_array_2 = @pos[$pos_starts[$b]..($pos_starts[$b+1]-1)];
                shift(@pos_array_1);
                shift(@pos_array_2);
                my @pos_combined_array_1_and_2 = (@pos_array_1, @pos_array_2);
                my @pos_unique_array_1_and_2_elements = unique(@pos_combined_array_1_and_2);

                my @neg_array_1 = @neg[$neg_starts[$a]..($neg_starts[$a+1]-1)];
                my @neg_array_2 = @neg[$neg_starts[$b]..($neg_starts[$b+1]-1)];
                shift(@neg_array_1);
                shift(@neg_array_2);
                my @neg_combined_array_1_and_2 = (@neg_array_1, @neg_array_2);
                my @neg_unique_array_1_and_2_elements = unique(@neg_combined_array_1_and_2);

                print "     - Comparing ", $array_1_name, " and ", $array_2_name, "\n";
                $array_1_name =~ s/_\w\w.MATS.JCEC.txt$|_\w\w\w.MATS.JCEC.txt$|_\w\w\w\w.MATS.JCEC.txt$//g;
                $array_2_name =~ s/_\w\w.MATS.JCEC.txt$|_\w\w\w.MATS.JCEC.txt$|_\w\w\w\w.MATS.JCEC.txt$//g;

                #print "\nMatching start";
                my %pos_array_1 = map { $_ => 1 } @pos_array_1;
                my @common_pos_junctions = grep { $pos_array_1{$_} } @pos_array_2;

                my %neg_array_1 = map { $_ => 1 } @neg_array_1;
                my @common_neg_junctions = grep { $neg_array_1{$_} } @neg_array_2;
                #print "\tMatching done\n";

                my @common_all_junctions = (@common_pos_junctions, @common_neg_junctions);
                
                my $detected_jncts_1 = scalar @array_1;
                my $detected_jncts_2 = scalar @array_2;
                my $sig_pos_IncDiff_array_1 = scalar @pos_array_1;
                @pos_array_2 = grep { $_ ne "null" } @pos_array_2;
                my $sig_pos_IncDiff_array_2 = scalar @pos_array_2;
                my $sig_neg_IncDiff_array_1 = scalar @neg_array_1;
                @neg_array_2 = grep { $_ ne "null" } @neg_array_2;
                my $sig_neg_IncDiff_array_2 = scalar @neg_array_2;
                my $significant_all_array_1 = $sig_pos_IncDiff_array_1 + $sig_neg_IncDiff_array_1;

                my $significant_all_array_2 = $sig_pos_IncDiff_array_2 + $sig_neg_IncDiff_array_2;

                if ($overlap_out == 1) {
                    open(OUT2, ">$output_dir_path_2/5_pairwise_overlap/1_pos_IncDiff/junctions_in.$array_1_name.and.$array_2_name.pos_IncDiff.tsv") or die "couldn't open input file";
                    open(OUT3, ">$output_dir_path_2/5_pairwise_overlap/2_neg_IncDiff/junctions_in.$array_1_name.and.$array_2_name.neg_IncDiff.tsv") or die "couldn't open input file";
                    if ($splice_junction_directory =~ m/A5SS_files/g || $splice_junction_directory =~ m/A3SS_files/g) {
                        my $header = "transcript\tgene\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\n";
                        print OUT2 $header, join("\n", @common_pos_junctions);
                        print OUT3 $header, join("\n", @common_neg_junctions);
                    }
                    elsif ($splice_junction_directory =~ m/MXE_files/g) {
                        my $header = "transcript\tgene\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";
                        print OUT2 $header, join("\n", @common_pos_junctions);
                        print OUT3 $header, join("\n", @common_neg_junctions);
                    }
                    elsif ($splice_junction_directory =~ m/RI_files/g) {
                        my $header = "transcript\tgene\tchr\tstrand\triExonStart_0base\triExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";
                        print OUT2 $header, join("\n", @common_pos_junctions);
                        print OUT3 $header, join("\n", @common_neg_junctions);
                    }
                    elsif ($splice_junction_directory =~ m/SE_files/g) {
                        my $header = "transcript\tgene\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";
                        print OUT2 $header, join("\n", @common_pos_junctions);
                        print OUT3 $header, join("\n", @common_neg_junctions);
                    }
                    close(OUT2);
                    close(OUT3);
                }

                my $common_all_junctions_count = scalar @common_all_junctions;
                #print "Hyper start\t";
                my $hypergeo_all = `python $dirname/accessory_scripts/hypergeometric_test $total_jncts_detected_across_samples $significant_all_array_1 $significant_all_array_2 $common_all_junctions_count`;

                my $common_pos_junctions_count = scalar @common_pos_junctions;
                my $hypergeo_pos = `python $dirname/accessory_scripts/hypergeometric_test $total_jncts_detected_across_samples $sig_pos_IncDiff_array_1 $sig_pos_IncDiff_array_2 $common_pos_junctions_count`;

                my $common_neg_junctions_count = scalar @common_neg_junctions;
                my $hypergeo_neg = `python $dirname/accessory_scripts/hypergeometric_test $total_jncts_detected_across_samples $sig_neg_IncDiff_array_1 $sig_neg_IncDiff_array_2 $common_neg_junctions_count`;
                chomp($hypergeo_all);
                chomp($hypergeo_pos);
                chomp($hypergeo_neg);
                #print "Hyper done\n";

                print OUT1a $array_1_name, "\t", $total_jncts_detected_across_samples, "\t", $detected_jncts_1, "\t\t", $significant_all_array_1, "\t", scalar (@common_all_junctions), "\t", (scalar(@common_all_junctions)+0.05)/($significant_all_array_1+0.05), "\t", $hypergeo_all, "\t\t", $sig_pos_IncDiff_array_1, "\t", scalar(@common_pos_junctions), "\t", (scalar(@common_pos_junctions)+0.05)/($sig_pos_IncDiff_array_1+0.05), "\t", $hypergeo_pos, "\t\t", $sig_neg_IncDiff_array_1, "\t", scalar(@common_neg_junctions), "\t", (scalar(@common_neg_junctions)+0.05)/($sig_neg_IncDiff_array_1+0.05), "\t", $hypergeo_neg, "\n";
                print OUT1a $array_2_name, "\t", $total_jncts_detected_across_samples, "\t", $detected_jncts_2 , "\t\t", $significant_all_array_2, "\t", scalar(@common_all_junctions), "\t", (scalar(@common_all_junctions)+0.05)/($significant_all_array_2+0.05), "\t", $hypergeo_all, "\t\t", $sig_pos_IncDiff_array_2, "\t", scalar(@common_pos_junctions), "\t", (scalar(@common_pos_junctions)+0.05)/($sig_pos_IncDiff_array_2+0.05), "\t", $hypergeo_pos, "\t\t", $sig_neg_IncDiff_array_2, "\t", scalar(@common_neg_junctions), "\t", (scalar(@common_neg_junctions)+0.05)/($sig_neg_IncDiff_array_2+0.05), "\t", $hypergeo_neg, "\n";

                print OUT1a "\n\n";

                push(@hypergeo_values_all_array, $hypergeo_all);
                push(@hypergeo_vals_pos_array, $hypergeo_pos);
                push(@hypergeo_vals_neg_array, $hypergeo_neg);

                if ($a == 0 and $b == 1) {
                    push(@hypergeo_output_matrix_name_array, $array_1_name, $array_2_name);
                }
                else {
                    push(@hypergeo_output_matrix_name_array, $array_2_name);    
                }

                my $output_4_data_array_1 = $array_1_name."\t".$significant_all_array_1."\t".$sig_pos_IncDiff_array_1."\t".($sig_pos_IncDiff_array_1+0.05)/($significant_all_array_1+0.05)."\t".$sig_neg_IncDiff_array_1."\t".($sig_neg_IncDiff_array_1+0.05)/($significant_all_array_1+0.05)."\t".scalar @array_1;

                if (grep (/^$output_4_data_array_1$/, @output_4_data_array)) {
                    next;
                }
                else {
                    push(@output_4_data_array, $output_4_data_array_1);
                }
                if ($a == $data_file_number-2 and $b == $data_file_number-1) {
                    my $output_4_data_array_2 = $array_2_name."\t".$significant_all_array_2."\t".$sig_pos_IncDiff_array_2."\t".($sig_pos_IncDiff_array_2+0.05)/($significant_all_array_2+0.05)."\t".$sig_neg_IncDiff_array_2."\t".($sig_neg_IncDiff_array_2+0.05)/($significant_all_array_2+0.05)."\t".scalar @array_2;
                    push(@output_4_data_array, $output_4_data_array_2);
                }
                #print "END CYCLE\n\n";
            }
        }
        close(OUT1a);
        my $files_number = @files;
        print OUT1b "\t", join("\t", @hypergeo_output_matrix_name_array[0..($files_number-1)]), "\n";
        print OUT1c "\t", join("\t", @hypergeo_output_matrix_name_array[0..($files_number-1)]), "\n";
        print OUT1d "\t", join("\t", @hypergeo_output_matrix_name_array[0..($files_number-1)]), "\n";
        print OUT1e "\t", join("\t", @hypergeo_output_matrix_name_array[0..($files_number-1)]), "\n";
        print OUT1f "\t", join("\t", @hypergeo_output_matrix_name_array[0..($files_number-1)]), "\n";
        print OUT1g "\t", join("\t", @hypergeo_output_matrix_name_array[0..($files_number-1)]), "\n";

        my $hypergeo_count = 0;
        for (my $y=0; $y<$files_number; $y++) {
            print OUT1b $hypergeo_output_matrix_name_array[$y];
            print OUT1c $hypergeo_output_matrix_name_array[$y];
            print OUT1d $hypergeo_output_matrix_name_array[$y];
            print OUT1e $hypergeo_output_matrix_name_array[$y];
            print OUT1f $hypergeo_output_matrix_name_array[$y];
            print OUT1g $hypergeo_output_matrix_name_array[$y];
        
            my $x = 0;
            for (my $z=0; $z<$files_number; $z++) {
                if ($z > $y) {
                    print OUT1e "\t", $hypergeo_values_all_array[$hypergeo_count];
                    print OUT1f "\t", $hypergeo_vals_pos_array[$hypergeo_count];
                    print OUT1g "\t", $hypergeo_vals_neg_array[$hypergeo_count];
                    if ($hypergeo_values_all_array[$hypergeo_count] > 0) {
                        print OUT1b "\t", -log10($hypergeo_values_all_array[$hypergeo_count]);
                    }
                    else {
                        print OUT1b "\t", $max_minus_log_10_hypergeo_val;
                    }
                    if ($hypergeo_vals_pos_array[$hypergeo_count] > 0) {
                        print OUT1c "\t", -log10($hypergeo_vals_pos_array[$hypergeo_count]);
                    }
                    else {
                        print OUT1c "\t", $max_minus_log_10_hypergeo_val;
                    }
                    if ($hypergeo_vals_neg_array[$hypergeo_count] > 0) {
                        print OUT1d "\t", -log10($hypergeo_vals_neg_array[$hypergeo_count]);
                    }
                    else {
                        print OUT1d "\t", $max_minus_log_10_hypergeo_val;
                    }
                    $hypergeo_count++;
                }
                elsif ($z == $y) {
                    print OUT1b "\t", $max_minus_log_10_hypergeo_val;
                    print OUT1c "\t", $max_minus_log_10_hypergeo_val;
                    print OUT1d "\t", $max_minus_log_10_hypergeo_val;
                    print OUT1e "\t0";
                    print OUT1f "\t0";
                    print OUT1g "\t0";
                }
                elsif ($z == 0) {
                    $x = $y-1;
                    print OUT1e "\t", $hypergeo_values_all_array[$x];
                    print OUT1f "\t", $hypergeo_vals_pos_array[$x];
                    print OUT1g "\t", $hypergeo_vals_neg_array[$x];

                    if ($hypergeo_values_all_array[$x] > 0) {
                        print OUT1b "\t", -log10($hypergeo_values_all_array[$x]);
                    }
                    else {
                        print OUT1b "\t", $max_minus_log_10_hypergeo_val;
                    }
                    if ($hypergeo_vals_pos_array[$x] > 0) {
                        print OUT1c "\t", -log10($hypergeo_vals_pos_array[$x]);
                    }
                    else {
                        print OUT1c "\t", $max_minus_log_10_hypergeo_val;
                    }
                    if ($hypergeo_vals_neg_array[$x] > 0) {
                        print OUT1d "\t", -log10($hypergeo_vals_neg_array[$x]);
                    }
                    else {
                        print OUT1d "\t", $max_minus_log_10_hypergeo_val;
                    }
                }
                else {
                    $x = $x+($files_number-($z+1));
                    print OUT1e "\t", $hypergeo_values_all_array[$x];
                    print OUT1f "\t", $hypergeo_vals_pos_array[$x];
                    print OUT1g "\t", $hypergeo_vals_neg_array[$x];
                    if ($hypergeo_values_all_array[$x] > 0) {
                        print OUT1b "\t", -log10($hypergeo_values_all_array[$x]);
                    }
                    else {
                        print OUT1b "\t", $max_minus_log_10_hypergeo_val;
                    }
                    if ($hypergeo_vals_pos_array[$x] > 0) {
                        print OUT1c "\t", -log10($hypergeo_vals_pos_array[$x]);
                    }
                    else {
                        print OUT1c "\t", $max_minus_log_10_hypergeo_val;
                    }
                    if ($hypergeo_vals_neg_array[$x] > 0) {
                        print OUT1d "\t", -log10($hypergeo_vals_neg_array[$x]);
                    }
                    else {
                        print OUT1d "\t", $max_minus_log_10_hypergeo_val;
                    }
                }
            }
            print OUT1b "\n";
            print OUT1c "\n";
            print OUT1d "\n";
            print OUT1e "\n";
            print OUT1f "\n";
            print OUT1g "\n";
        }
        close(OUT1b);
        close(OUT1c);
        close(OUT1d);
        close(OUT1e);
        close(OUT1f);
        close(OUT1g);
        print OUT1 join("\n", @output_4_data_array);
        close(OUT1);

        my $clustermap_all = $output_dir_path_2."/3_data_matrix_minus_log10_pval/1_data_matrix_all.tsv";
        my $clustermap_pos = $output_dir_path_2."/3_data_matrix_minus_log10_pval/2_data_matrix_pos_IncDiff.tsv";
        my $clustermap_neg = $output_dir_path_2."/3_data_matrix_minus_log10_pval/3_data_matrix_neg_IncDiff.tsv";
        my $clustmap_path = $dirname."/accessory_scripts/clustermap.py";

        `python $clustmap_path $clustermap_all`;
        `python $clustmap_path $clustermap_pos`;
        `python $clustmap_path $clustermap_neg`;

        `mv $output_dir_path_2/3_data_matrix_minus_log10_pval/\*label_order $output_dir_path_2/4_clustermaps/`;
        `mv $output_dir_path_2/3_data_matrix_minus_log10_pval/\*tsv.svg $output_dir_path_2/4_clustermaps/`;
    }
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

sub unique {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

options;
qc;
read_input_dir;
make_output_dirs;
process_JCEC_files;
summarize;
compare;

`rm -r $out_dir/*/1_rMATS_sig/1_original`;
print "\nDONE! with all analyses\n\n";


