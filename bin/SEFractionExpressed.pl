use warnings;
use strict; 
use File::Basename;       

my $skipped_exon_file;
my $expression_file;
my $TPM1;
my $TPM2;
my $cond1_ct;
my $cond2_ct;
my $FDR;

my $outdir;
my $basename;
my $expression_basename;

my %expressed_genes;
my @pos_IncDiff;
my @neg_IncDiff;
my @neg_IncDiff_genes = ();
my @pos_IncDiff_genes = ();
my @no_SE_genes = ();

sub program_info {
    print "\n\tSEFractionExpressed.pl determines the fraction of expressed genes with statistically\n\tsignificant SE events at genes with a minimum input control or test TPM value.\n\n\tAlso generates:\n\t\t- lists of statistically significant negative and positive IncDiff SE events in\n\t\t  expressed genes (includes gene expression values for each)\n\t\t- file with list of genes without statistcially significant exon skipping.\n\n\n\tUsage: perl SEFractionExpressed.pl [OPTIONS] -s <skipped exon file (rMATS JCEC)> -e <expression file> -TPM <min TPMs condition 1,min TPMs condition 2> -SN <sample number condition 1,sample number condition 2> -f <FDR>\n\n\tRequired:\n\t\t-s <skipped exon file>\n\t\t-e <expression file>\n\t\t-TPM <min TPMs condition 1,min TPMs condition 2>\n\t\t\tNote: if condition 1 or 2 to be not considered, enter \"\-\"\n\t\t\t(e.g. TPM 3,- or TPM -,3)\n\t\t-SN <sample number condition 1,sample number condition 2>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl SEFractionExpressed.pl -s PATH/SEfile.txt -e PATH/expression.tsv -TPM 2,2 -SN 3,3 -f 0.05\n\n";
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
        elsif ($ARGV[$i] eq "\-e") {
            $expression_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-TPM") {
            my $TPMs = $ARGV[$i+1];
            ($TPM1, $TPM2) = split("\,", $TPMs);
            if ($TPM1 eq "-") {
                $TPM1 = 1000001;
            }
            if ($TPM2 eq "-") {
                $TPM2 = 1000001;
            }
        }
        elsif ($ARGV[$i] eq "\-SN") {
            my $sample_counts = $ARGV[$i+1];
            ($cond1_ct, $cond2_ct) = split("\,", $sample_counts);
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
        print "\nFDR not defined!\n\n";
        program_info;
        exit;
    }
    elsif ($FDR !~ m/^\d/) {
        print "\nFDR is not numeric!\n\n";
        program_info;
        exit;
    } 
    elsif (not defined($skipped_exon_file)) {
        print "\nSkipped exon file not defined!\n\n";
        program_info;
        exit;
    }
    elsif (not defined($expression_file)) {
        print "\nExpression file not defined!\n\n";
        program_info;
        exit;
    }
    elsif (not defined($TPM1)) {
        print "\nMinimum TPMs for condition 1 not defined!\n\tIs -TPM value in correct format?\n\n";
        program_info;
        exit;
    }
    elsif (not defined($TPM2)) {
        print "\nMinimum TPMs for condition 2 not defined!\n\tIs -TPM value in correct format?\n\n";
        program_info;
        exit;
    }
    elsif (not defined($cond1_ct)) {
        print "\nCondition 1 sample count not defined!\n\tIs -SN value in correct format?\n\n";
        program_info;
        exit;
    }
    elsif (not defined($cond2_ct)) {
        print "\nCondition 2 sample count not defined!\n\tIs -SN value in correct format?\n\n";
        program_info;
        exit;
    }
    open(INF, "<$skipped_exon_file") or die "couldn't open skipped exon file\n\n";
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
                print "\nSkipped exon file may not be an rMATS JCEC file format!\n\tIs fourth column chromosome? Fourth column doesn't begin with \"chr\"\n\n";
                program_info;
                exit;
            }
        }
        elsif ($. == 3) {
            last;
        }
    }
    close(INF);

    open(INF, "<$expression_file") or die "couldn't open expression file";
    while(my $line = <INF>) {
        chomp($line);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        next if ($. == 1);
        if (scalar @split_line < (1+$cond1_ct+$cond2_ct) || $split_line[1] !~ m/^\d/ || $split_line[$cond1_ct+$cond2_ct] !~ m/^\d/) {
            print "\n     Expression file may not be correct format!\n\t\t- Also check input numbers of samples for each condition.\n\n";
            program_info;
            exit;
        }
    }
    close(INF);
}

sub makedirectory {
    $basename = basename($skipped_exon_file);
    $basename =~ s/\.txt$//g;
    $expression_basename = basename($expression_file);
    $outdir = dirname($skipped_exon_file)."\/SEFractionExpressed_".$basename."_cntl_TPM_".$TPM1."_test_TPM_".$TPM2."_min_rMATS_FDR_".$FDR;
    `mkdir $outdir`;
    print "\nINPUT:\nSkipped exon file: ", $basename, "\nExpression file: ", $expression_basename, "\nMin TPMs, condition 1 = ", $TPM1, "\nMin TPMs, condition 2 = ", $TPM2, "\nSample count 1 = ", $cond1_ct, "\nSample count 2 = ", $cond2_ct, "\nSE FDR = ", $FDR, "\n";
}

sub process_expression_file {
    print "\nFinding expressed genes...\n\n";
    open(INF, "<$expression_file") or die "couldn't open input file";

    while (my $line = <INF>) {
	    chomp($line);
        my @split_line = split("\t", $line);
        next if ($. == 1);
        my $sum_cond_1 = 0; 
        for (my $j = 0; $j < $cond1_ct; $j++) {
            $sum_cond_1 += $split_line[1+$j];
        }
        my $sum_cond_2 = 0;
        for (my $k = 0; $k < $cond2_ct; $k++) {
            #if (not defined $split_line[1+$cond1_ct+$k] || $split_line[1+$cond1_ct+$k] =~ m/"N"/g) {
                #print "\n\nERROR!!\nError reading expression data! Are the correct number of control and test samples inputted?\n Exiting SEFraction expressed.pl\n\n";
                #exit;
            #}
            $sum_cond_2 += $split_line[1+$cond1_ct+$k];
        }
        my $average_cond_1 = $sum_cond_1/$cond1_ct;
        my $average_cond_2 = $sum_cond_2/$cond2_ct;
        next if ($average_cond_1 < $TPM1 and $average_cond_2 < $TPM2);
        $expressed_genes{$split_line[0]} = $split_line[0]."\t".$average_cond_1."\t".$average_cond_2."\t".$average_cond_2/($average_cond_1+0.5);
    }
    close(INF);
}

sub process_SE_file {
    print "Identifying significant SE events...\n\n";
    open(INF, "<$skipped_exon_file") or die "couldn't open input file";
    while (my $line = <INF>) {
	    chomp($line);
        next if ($. == 1);  
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        if ($split_line[19] < $FDR) {
            if ($split_line[22] > 0) {
                push(@pos_IncDiff, $split_line[3]."\t".$split_line[8]."\t".$split_line[9]."\t".$split_line[1]."\_".$split_line[2]."\t1000\t".$split_line[4]."\t".$split_line[19]."\t".$split_line[22]);
            }
            elsif ($split_line[22] < 0) {
                push(@neg_IncDiff, $split_line[3]."\t".$split_line[8]."\t".$split_line[9]."\t".$split_line[1]."\_".$split_line[2]."\t1000\t".$split_line[4]."\t".$split_line[19]."\t".$split_line[22]);
            }
        }
    }
    close(INF);
}

sub expressed_with_SE {
    print "Matching SE data with expression data...\n\n";
    my $sig_neg_IncDiff_rMATS_out = $outdir."\/1_sig_neg_IncDiff.tsv";
    my $sig_pos_IncDiff_rMATS_out = $outdir."\/2_sig_pos_IncDiff.tsv";
    my $genes_with_no_sig_SEs_out = $outdir."\/3_genes_without_sig_SE.tsv";

    open(OUT1, ">$sig_neg_IncDiff_rMATS_out") or die "couldn't open output file";
    open(OUT2, ">$sig_pos_IncDiff_rMATS_out") or die "couldn't open output file";
    open(OUT3, ">$genes_with_no_sig_SEs_out") or die "couldn't open output file";
    print OUT1 "chr\tuSE_donor\tdSE_acceptor\tIsoform_Gene\tIntensity\tStrand\tFDR\tIncDiff\tGene\tAverage cntl expression\tAverage test expression\tCondition 1(TPM)/(Condition 2(TPM) + 0.5)\n";
    print OUT2 "chr\tuSE_donor\tdSE_acceptor\tIsoform_Gene\tIntensity\tStrand\tFDR\tIncDiff\tGene\tAverage cntl expression\tAverage test expression\tCondition 1(TPM)/(Condition 2(TPM) + 0.5)\n";
    print OUT3 "Gene\tCondition 1(TPM)\tCondition 2(TPM)\tCondition 1(TPM)/(Condition 2(TPM) + 0.5)\n";

    foreach my $neg_element(@neg_IncDiff) {

        my @split_element = split("\t", $neg_element);
        my @split_ID = split("_", $split_element[3]);
        if (exists $expressed_genes{$split_ID[1]}) {
            print OUT1 $neg_element, "\t", $expressed_genes{$split_ID[1]}, "\n";
            push(@neg_IncDiff_genes, $split_ID[1]);
        }
    }
    foreach my $pos_element(@pos_IncDiff) {
        my @split_element = split("\t", $pos_element);
        my @split_ID = split("_", $split_element[3]);
        if (exists $expressed_genes{$split_ID[1]}) {
            print OUT2 $pos_element, "\t", $expressed_genes{$split_ID[1]}, "\n";
            push(@pos_IncDiff_genes, $split_ID[1]);
        }
    }
    foreach my $neg_element(@neg_IncDiff) {

        my @split_element = split("\t", $neg_element);
        my @split_ID = split("_", $split_element[3]);
        if (exists $expressed_genes{$split_ID[1]}) {
            delete $expressed_genes{$split_ID[1]};
        }
    }  
    foreach my $pos_element(@pos_IncDiff) {
        my @split_element = split("\t", $pos_element);
        my @split_ID = split("_", $split_element[3]);
        if (exists $expressed_genes{$split_ID[1]}) {
            delete $expressed_genes{$split_ID[1]};
        }
    }
    foreach (keys %expressed_genes) {
        print OUT3 $expressed_genes{$_}, "\n";
        push(@no_SE_genes, $_);
    }
    close(OUT1);
    close(OUT2);
    close(OUT3);
}

sub statistics {
    print "Calculating fraction of genes with SE...\n";
    my @unique_neg_IncDiff_genes = unique(@neg_IncDiff_genes);
    my $num_unique_neg_IncDiff_genes = @unique_neg_IncDiff_genes;
    my @unique_pos_IncDiff_genes = unique(@pos_IncDiff_genes);
    my $num_unique_pos_IncDiff_genes = @unique_pos_IncDiff_genes;
    my @combined_neg_and_pos_IncDiff_genes = (@neg_IncDiff_genes, @pos_IncDiff_genes);
    my @unique_combined_neg_and_pos_IncDiff_genes = unique(@combined_neg_and_pos_IncDiff_genes);
    my $number_unique_combined_neg_and_pos_IncDiff_genes = @unique_combined_neg_and_pos_IncDiff_genes;
    my $number_no_SE_genes = @no_SE_genes;
    my $total_genes = $number_unique_combined_neg_and_pos_IncDiff_genes + $number_no_SE_genes; 

    my $stats_out = $outdir."\/1_1_fraction_of_genes_with_SE.txt";

    open(OUT, ">$stats_out") or die "couldn't open stats output file";
    print OUT "Fractions genes with SE relative to all expressed genes:\n\nSkipped exon file: ", $basename, "\nExpression file: ", $expression_basename, "\n\nMinimum expression condition 1 = ", $TPM1, " -or- Minimum expression condition 2 = ", $TPM2, "\nCondition 1 sample number = ", $cond1_ct, "\nCondition 2 sample number = ", $cond2_ct,  "\nSE FDR cutoff = ", $FDR, "\n\n\nExpressed genes with neg IncDiff = ", $num_unique_neg_IncDiff_genes,  "\nExpressed genes with pos IncDiff = ", $num_unique_pos_IncDiff_genes, "\nExpressed genes with no significant SE = ", $number_no_SE_genes, "\nTotal expressed genes = ", $total_genes, "\n\nFraction expressed genes with neg IncDiff = ", sprintf ("%.4f", $num_unique_neg_IncDiff_genes/$total_genes),
    "\nFraction expressed genes with pos IncDiff = ", sprintf("%.4f", $num_unique_pos_IncDiff_genes/$total_genes), "\n";
    close(OUT);
}

sub unique {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

options;
qc;
print "\n\n***********************\nRunning SEFractionExpressed...\n";
makedirectory;
process_expression_file;
process_SE_file;
expressed_with_SE;
statistics;
print "\nSEFractionExpressed DONE!\n***********************\n\n";