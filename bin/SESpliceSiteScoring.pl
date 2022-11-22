use warnings;
use strict;  
use File::Basename;

my $skipped_exon_file;
my $genome_fasta_file;
my $annotation_file = "none";
my $FDR;

my $script_dirname = dirname(__FILE__);
my @genome_array = ();
my $genome_array_length;
my @annotation_array = ();
my $outdir;
my $basename;

sub program_info {
    print "\n\tSESpliceSiteScoring.pl will generate:\n\t\t- lists of splice site scores (for plotting score distributions) for upstream donor,\n\t\t  skipped exon, and downstream acceptor sites for statistically significantly changed\n\t\t  exon skipping events.\n\t\t- summary file with average scores.\n\n\tNote: Inclusion of annotation file is optional but will generate data for annotated events\n\tfor comparison.\n\n\tUsage: perl SESpliceSiteScoring.pl [OPTIONS] -s <skipped exon file (rMATS JCEC)> -g <genome fasta file> -f <FDR> -a <bed12 annotation file>\n\n\tRequired:\n\t\t-s <skipped exon file>\n\t\t-g <genome fasta file>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-a <bed12 annotation file> (optional)\n\t\t-h help\n\n\tExample: perl SESpliceSiteScoring.pl -s PATH/SEfile.txt -g PATH/genome.fa -a PATH/bed12_annotation.bed -f 0.05 \n\n";
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
        elsif ($ARGV[$i] eq "\-g") {
            $genome_fasta_file = $ARGV[$i+1];
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
    elsif (not defined($genome_fasta_file)) {
        print "\nGenome fasta file not defined!\n\n";
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
    open(INF, "<$genome_fasta_file") or die "couldn't open genome fasta file";
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        if ($. == 1) {
            if (scalar @split_line != 1) {
                print "\nGenome fasta file appears to be incorrect format!\n\n";
                program_info;
                exit;
            }
            elsif ($line =~ m/^\>chr/ || $line =~ m/^\> chr/) {
                last;
            }
            else {
                print "\nGenome fasta file appears to be incorrect format!\n     First line doesn't match \"\>chr\" or \"\> chr\"\n\n";
                program_info;
                exit;
            }
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
    $outdir = dirname($skipped_exon_file)."\/SESpliceSiteScoring_".$basename."_FDR_".$FDR;
    `mkdir $outdir`;
    `mkdir $outdir"\/1_neg_IncDiff"`;
    `mkdir $outdir"\/2_pos_IncDiff"`;
    `mkdir $outdir"\/3_annotated"`;
}

sub print_input_parameters {
    open(OUT, ">$outdir\/1_1_input_parameters.txt") or die "couldn't open run parameters output file";

    print OUT "Skipped exon file: ", $basename, "\nAnnotation file: ", basename($annotation_file), "\nGenome fasta file: ", basename($genome_fasta_file), "\nFDR: ", $FDR, "\n\n";
    print "INPUT\nSkipped exon file: ", $basename, "\nAnnotation file: ", basename($annotation_file), "\nGenome fasta file: ", basename($genome_fasta_file), "\nFDR: ", $FDR, "\n";
    close(OUT);
}

sub genome_fasta_file_unwrapper {
    print "\nProcessing genome fasta file...\n\n"; #Unwrapping and putting into array
    open(INF, "<$genome_fasta_file") or die "couldn't open genome fasta file";
    my @line_array = ();
    while(my $line = <INF>) {
        chomp($line);
        if ($. == 1) {
            push(@genome_array, $line);
        }
        elsif ($line =~ m/^\>chr/) {
            push(@genome_array, join("", @line_array));
            my @split_line = split(" ", $line);
            push(@genome_array, $split_line[0]);
            @line_array = ();
        }
        elsif (eof(INF)) {
            push(@line_array, $line);
            push(@genome_array, join("", @line_array));
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
    my $upstream_exon_donor_neg_IncDiff = "\>upstream_exon_donor_neg_IncDiff\n";
    my $skipped_exon_acceptor_neg_IncDiff = "\>skipped_exon_acceptor_neg_IncDiff\n";
    my $skipped_exon_donor_neg_IncDiff = "\>skipped_exon_donor_neg_IncDiff\n";
    my $downstream_exon_acceptor_neg_IncDiff = "\>downstream_exon_acceptor_neg_IncDiff\n";

    my $upstream_exon_donor_pos_IncDiff = "\>upstream_exon_donor_pos_IncDiff\n";
    my $skipped_exon_acceptor_pos_IncDiff = "\>skipped_exon_acceptor_pos_IncDiff\n";
    my $skipped_exon_donor_pos_IncDiff = "\>skipped_exon_donor_pos_IncDiff\n";
    my $downstream_exon_acceptor_pos_IncDiff = "\>downstream_exon_acceptor_pos_IncDiff\n";

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
                            $upstream_exon_donor_neg_IncDiff = $upstream_exon_donor_neg_IncDiff.uc(substr $genome_array[$i], $split_line[8]-3, 9)."\n";
                            $skipped_exon_acceptor_neg_IncDiff = $skipped_exon_acceptor_neg_IncDiff.uc(substr $genome_array[$i], $split_line[5]-20, 23)."\n";
                            $skipped_exon_donor_neg_IncDiff = $skipped_exon_donor_neg_IncDiff.uc(substr $genome_array[$i], $split_line[6]-3, 9)."\n";
                            $downstream_exon_acceptor_neg_IncDiff = $downstream_exon_acceptor_neg_IncDiff.uc(substr $genome_array[$i], $split_line[9]-20, 23)."\n";
                        }
                        elsif ($amplitude > 0) {
                            $upstream_exon_donor_pos_IncDiff = $upstream_exon_donor_pos_IncDiff.uc(substr $genome_array[$i], $split_line[8]-3, 9)."\n";
                            $skipped_exon_acceptor_pos_IncDiff = $skipped_exon_acceptor_pos_IncDiff.uc(substr $genome_array[$i], $split_line[5]-20, 23)."\n";
                            $skipped_exon_donor_pos_IncDiff = $skipped_exon_donor_pos_IncDiff.uc(substr $genome_array[$i], $split_line[6]-3, 9)."\n";
                            $downstream_exon_acceptor_pos_IncDiff = $downstream_exon_acceptor_pos_IncDiff.uc(substr $genome_array[$i], $split_line[9]-20, 23)."\n";
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
                            $upstream_exon_donor_neg_IncDiff = $upstream_exon_donor_neg_IncDiff.$upstream_donor."\n";
                            $skipped_exon_acceptor_neg_IncDiff = $skipped_exon_acceptor_neg_IncDiff.$SE_acceptor."\n";
                            $skipped_exon_donor_neg_IncDiff = $skipped_exon_donor_neg_IncDiff.$SE_donor."\n";
                            $downstream_exon_acceptor_neg_IncDiff = $downstream_exon_acceptor_neg_IncDiff.$downstream_acceptor."\n";
                        }
                        elsif ($amplitude > 0) {
                            $upstream_exon_donor_pos_IncDiff = $upstream_exon_donor_pos_IncDiff.$upstream_donor."\n";
                            $skipped_exon_acceptor_pos_IncDiff = $skipped_exon_acceptor_pos_IncDiff.$SE_acceptor."\n";
                            $skipped_exon_donor_pos_IncDiff = $skipped_exon_donor_pos_IncDiff.$SE_donor."\n";
                            $downstream_exon_acceptor_pos_IncDiff = $downstream_exon_acceptor_pos_IncDiff.$downstream_acceptor."\n";
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

    my $neg_IncDiff_directory = $outdir."\/1_neg_IncDiff\/neg_IncDiff";
    my $pos_IncDiff_directory = $outdir."\/2_pos_IncDiff\/pos_IncDiff";

    open(OUT1, ">$neg_IncDiff_directory.upstream_donor_neg") or die "couldn't open upstream_donor_neg file";
    open(OUT2, ">$neg_IncDiff_directory.SE_acceptor_neg") or die "couldn't open SE_acceptor_neg file"; 
    open(OUT3, ">$neg_IncDiff_directory.SE_donor_neg") or die "couldn't open SE_donor_neg file"; 
    open(OUT4, ">$neg_IncDiff_directory.downstream_acceptor_neg") or die "couldn't open downstream_acceptor_neg file";  
    open(OUT5, ">$pos_IncDiff_directory.upstream_donor_pos") or die "couldn't open upstream_donor_pos file";
    open(OUT6, ">$pos_IncDiff_directory.SE_acceptor_pos") or die "couldn't open SE_acceptor_pos file"; 
    open(OUT7, ">$pos_IncDiff_directory.SE_donor_pos") or die "couldn't open SE_donor_pos file"; 
    open(OUT8, ">$pos_IncDiff_directory.downstream_acceptor_pos") or die "couldn't open downstream_acceptor_pos file";  

    print OUT1 $upstream_exon_donor_neg_IncDiff;
    print OUT2 $skipped_exon_acceptor_neg_IncDiff;
    print OUT3 $skipped_exon_donor_neg_IncDiff;
    print OUT4 $downstream_exon_acceptor_neg_IncDiff;
    print OUT5 $upstream_exon_donor_pos_IncDiff;
    print OUT6 $skipped_exon_acceptor_pos_IncDiff;
    print OUT7 $skipped_exon_donor_pos_IncDiff;
    print OUT8 $downstream_exon_acceptor_pos_IncDiff;

    close(OUT1);
    close(OUT2);
    close(OUT3);
    close(OUT4);
    close(OUT5);
    close(OUT6);
    close(OUT7);
    close(OUT8);

    print OUTsummary "\tupstream exon donor ave\tSE acceptor ave\tSE donor ave\tdownstream exon acceptor ave\n";
    print OUTsummary "neg IncDiff\t", score5($neg_IncDiff_directory.".upstream_donor_neg"), "\t", score3($neg_IncDiff_directory.".SE_acceptor_neg"), "\t", score5($neg_IncDiff_directory.".SE_donor_neg"), "\t", score3($neg_IncDiff_directory.".downstream_acceptor_neg"), "\npos IncDiff\t", score5($pos_IncDiff_directory.".upstream_donor_pos"), "\t", score3($pos_IncDiff_directory.".SE_acceptor_pos"), "\t", score5($pos_IncDiff_directory.".SE_donor_pos"), "\t", score3($pos_IncDiff_directory.".downstream_acceptor_pos"), "\n";
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

    my $annotated_directory = $outdir."\/3_annotated\/annotated";

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

    print OUTsummary "\nannotated\t", score5($annotated_directory.".upstream_donor"), "\t", score3($annotated_directory.".middle_acceptor"), "\t", score5($annotated_directory.".middle_donor"), "\t", score3($annotated_directory.".downstream_acceptor"); 
}

sub score3 {
    my $score3 = $script_dirname."\/score3_mod.pl";
    foreach my $input(@_) {
        my $output_file = $input."_scores.tsv";
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
    my $score5 = $script_dirname."\/score5_mod.pl";
    foreach my $input(@_) {
        my $output_file = $input."_scores.tsv";
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
    if ($array_length == 0) {
        return "no significant SEs";
    }
    else {
        return $total/$array_length;
    }
}

options;
qc;
print "\n\n***********************\nRunning SESpliceSiteScoring...\n\n";
makedirectories;
print_input_parameters;
genome_fasta_file_unwrapper;

my $summary_output_file = $outdir."\/1_1_summary.tsv";
open(OUTsummary, ">$summary_output_file") or die "couldn't open upstream_donor_neg file";

SE_scoring;
if (-e $annotation_file) {
    prepare_annotation_file;
    annotated_scoring;
}
close(OUTsummary);

print "\nSESpliceSiteScoring DONE!\n***********************\n\n";