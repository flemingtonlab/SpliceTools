use warnings;
use strict;  
use File::Basename;

my $skipped_exon_file;
my $annotation_file;
my $genome_fasta_file;
my $FDR;
my @skipped_exon_array = ();
my @annotation_array = ();
my @merged_SE_ann_array = ();
my @genome_array = ();
my $script_dirname = dirname(__FILE__);

my %genetic_code;
my $basename;
my $outdir;

sub program_info {
    print "\n\tSETranslateNMD.pl will:\n\t\t- generate predicted RNA(DNA) isoform sequences based on annotated exon\n\t\t  structures upstream and downstream from skipping event\n\t\t- generate translated isoform sequences for all such events\n\t\t- generate separate listings of frameshifted protein isoforms\n\t\t- generate candidate neopeptides produced by frameshifts\n\t\t- generate lists of skipped exon protein sequences for in-frame\n\t\t  events for BATCH submission to NCBI conserved domain search\n\t\t  (https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi)\n\t\t- Identify events predicted to undergo Nonsense Mediated RNA\n\t\t  Decay (NMD) and generate statistics on skipping events predicted\n\t\t  to undergo and not undergo NMD\n\n\n\tUsage: perl SETranslateNMD.pl [OPTIONS] -s <skipped exon file (rMATS)> -a <bed12 annotation file> -g <genome fasta file> -f <FDR>\n\n\tRequired:\n\t\t-s <skipped exon file>\n\t\t-a <bed12 annotation file>\n\t\t-g <genome fasta file>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl SETranslateNMD.pl -s PATH/SEfile.txt -a PATH/bed12_annotation.bed -g PATH/genome.fa -f 0.05\n\n";
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
    if ($FDR !~ m/^\d/) {
        print "\nFDR is not numeric!\n\n";
        program_info;
        exit;
    } 
    elsif (not defined($skipped_exon_file)) {
        print "\nSkipped exon file directory not defined!\n\n";
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
    $outdir = dirname($skipped_exon_file)."\/SETranslateNMD_".$basename."_FDR_".$FDR;
    `mkdir $outdir`;
    `mkdir $outdir"\/1_nucleotide_sequences"`;
    `mkdir $outdir"\/2_SE_isoforms_all"`;
    `mkdir $outdir"\/2_SE_isoforms_all/temp"`;
    `mkdir $outdir"\/3_frameshifted_isoforms"`;
    `mkdir $outdir"\/3_frameshifted_isoforms/temp"`;
    `mkdir $outdir"\/4_neopeptides"`;
    `mkdir $outdir"\/4_neopeptides/temp"`;
    `mkdir $outdir"\/5_skipped_exon_protein_sequences"`;
    `mkdir $outdir"\/5_skipped_exon_protein_sequences/temp"`;
    `mkdir $outdir"\/6_NMD"`;
    `mkdir $outdir"\/6_NMD/SE_NMD_lists.bed"`;
    `mkdir $outdir"\/6_NMD/gene_lists"`;
}

sub print_input_parameters {
    open(OUT, ">$outdir\/1_1_input_parameters.txt") or die "couldn't open run parameters output file";

    print OUT "Skipped exon file: ", $basename, "\nAnnotation file: ", basename($annotation_file), "\nGenome fasta file: ", basename($genome_fasta_file), "\nFDR: ", $FDR, "\n\n";
    print "\nINPUT:\nSkipped exon file: ", $basename, "\nAnnotation file: ", basename($annotation_file), "\nGenome fasta file: ", basename($genome_fasta_file), "\nFDR: ", $FDR, "\n\n";
    close(OUT);
}

sub process_skipped_exon_array {
    print "Pocessing skipped exon file...\n\n";
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

    open(INF, "<$skipped_exon_file.sorted.tsv") or die "couldn't open sorted output skipped exon file";
    my @amp_arr = ();
    my $prev_line_minus_amp;
    my $prev_strand;
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        my $line_minus_amp = join("\t", @split_line[0..3]);
        my $amp = $split_line[4];
        my $strand = $split_line[5];
        if ($. == 1) {
            push(@amp_arr, $amp);
            $prev_line_minus_amp = $line_minus_amp;
            $prev_strand = $strand;
        }
        elsif ($line_minus_amp eq $prev_line_minus_amp) {
            push(@amp_arr, $amp);
        }
        else {
            push(@skipped_exon_array, $prev_line_minus_amp."\t".(eval join '+', @amp_arr)/scalar (@amp_arr)."\t".$prev_strand);
            $prev_line_minus_amp = $line_minus_amp;
            @amp_arr = ();
            push(@amp_arr, $amp);
            $prev_strand = $strand;
        }
    }
    close(INF);
    `rm $skipped_exon_file.sorted.tsv`;
}

sub process_annotation_file {
    print "Processing annotation file...\n\n"; #Sorting and putting into array
    `sort -k1,1 -k2,2n -k3,3n $annotation_file > $annotation_file.sorted.tsv`;

    open(INF, "<$annotation_file.sorted.tsv") or die "couldn't open sorted annotation file";
    while(my $line = <INF>) {
        chomp($line);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        next if ($split_line[1] == $split_line[6] || $split_line[2] == $split_line[7]);
        push(@annotation_array, $line);
    }
    `rm $annotation_file.sorted.tsv`;
    close(INF);
}

sub link_annotation_to_skipped_exons {
    print "Linking skipped exons to annotated transcripts...\n\n";
    my $k = 0;
    for (my $i=0; $i < scalar (@skipped_exon_array); $i++) {
        my @split_SE_element = split("\t", $skipped_exon_array[$i]);
        my $hit_count = 0;
        my $loop_after_hit_count = 0;
        for (my $j = $k; $j < scalar (@annotation_array); $j++) {
            my @split_ann_element = split("\t", $annotation_array[$j]);
            next if ($split_ann_element[6] == $split_ann_element[7]); # Exclude non-coding transcripts
            if ($hit_count > 0) {
                $loop_after_hit_count++;
            }
            if ($split_SE_element[0] eq $split_ann_element[0] and $split_SE_element[5] eq $split_ann_element[5] and $split_SE_element[1] > $split_ann_element[6] and $split_SE_element[2] < $split_ann_element[7]) {
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
                }
                $hit_count++;
                $k = $j-500;
            }
            elsif ($hit_count > 0 and $loop_after_hit_count > 500) {
                last;
            }
        }
    }
}

sub genome_fasta_file_unwrapper {
    print "Processing genome fasta file...\n\n"; #Unwrapping and putting into array
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
}

sub bed2fasta {
    print "Extracting nucleotide sequences for each transcript...\n";
    open(OUT1, ">$outdir/1_nucleotide_sequences/transcripts.fa") or die "couldn't open transcript sequence output file";
    open(OUT2, ">$outdir/1_nucleotide_sequences/coding_region.fa") or die "couldn't open coding region sequence output file";
    open(OUT3, ">$outdir/1_nucleotide_sequences/SE_transcripts_neg_IncDiff.fa") or die "couldn't open SE negative IncDiff output nucleotide fasta file";
    open(OUT4, ">$outdir/1_nucleotide_sequences/SE_transcripts_pos_IncDiff.fa") or die "couldn't open SE positive IncDiff output nucleotide fasta file";
    open(OUT5, ">$outdir/1_nucleotide_sequences/SKIPPED_EXONs_neg_IncDiff.fa") or die "couldn't open skipped exons fasta output file";
    open(OUT6, ">$outdir/1_nucleotide_sequences/SKIPPED_EXONs_pos_IncDiff.fa") or die "couldn't open skipped exons fasta output file";
    
    foreach my $elem(@merged_SE_ann_array) {
        next if $elem =~ m/^track/;
        my @split_elem = split("\t", $elem);
        my $chr = $split_elem[0];
        my $gene_name = $chr."\,".$split_elem[14]."\,".$split_elem[15]."\,".$split_elem[3]."\,1000\,".$split_elem[5];
        my $strand = $split_elem[5];
        my $trans_start = $split_elem[1];
        my $trans_end = $split_elem[2];
        my $coding_start = $split_elem[6];
        my $coding_end = $split_elem[7];
        my $block_number = $split_elem[9];
        my @block_sizes = split("\,", $split_elem[10]);
        my @block_starts = split("\,", $split_elem[11]);
        my $chr_hit_test = 0;

        my $SE_donor = $split_elem[14];
        my $SE_acceptor = $split_elem[15];
        my $IncDiff = $split_elem[17];
        my @skipped_exons_size_array = ();
        my @skipped_exons_start_array;
        my $inframe_skipped_exons_sequence = "";
        my $skipped_exons_length = 0;
        
        for (my $ga = 0; $ga < scalar @genome_array; $ga++) {
            my $gelem = $genome_array[$ga];
            if ($chr_hit_test > 1) {
                last;
            }
            if ($gelem =~ m/^\>/) {
                if ($gelem =~ m/$chr$/g and $chr_hit_test == 0) {
                    $chr_hit_test++;
                }
            }
            elsif ($chr_hit_test == 1) {
                my $transcript_sequence = "";
                my $coding_region_sequence = "";
                my $SE_transcript_sequence = "";
                my $SE_coding_start = 0;
                my @SE_block_sizes = ();
                my $SE_donor_junction = 0;
                $chr_hit_test++;
                for (my $i = 0; $i < $block_number; $i++) {
                    my $tblock_start = $trans_start + $block_starts[$i];
                    my $tblock_end = $tblock_start + $block_sizes[$i];

                    my $block_sequence = substr $gelem, $tblock_start, $block_sizes[$i];
                    $transcript_sequence = $transcript_sequence.$block_sequence;
                    #
                    if (($strand eq "+" and $coding_start > $SE_donor and $coding_start < $SE_acceptor) || ($strand eq "-" and $coding_end > $SE_donor and $coding_end < $SE_acceptor)) {
                        $SE_transcript_sequence = "INIT CODON SKIPPED:DEPPIKS NODOC TINI - ";
                    }
                    elsif ($tblock_end <= $SE_donor || $tblock_start >= $SE_acceptor) {
                        $SE_transcript_sequence = $SE_transcript_sequence.$block_sequence;
                        push(@SE_block_sizes, $block_sizes[$i]);
                        if ($strand eq "+" and $tblock_end < $coding_start) {
                            $SE_coding_start = $SE_coding_start + $block_sizes[$i];
                        }
                        elsif ($strand eq "-" and $tblock_end < $coding_end) {
                            $SE_coding_start = $SE_coding_start + $block_sizes[$i];
                        }
                        elsif ($strand eq "+" and $tblock_end >= $coding_start and $tblock_start <= $coding_start) {
                            my $distance_in_block_to_coding_start = $coding_start - $tblock_start;
                            $SE_coding_start = $SE_coding_start + $distance_in_block_to_coding_start;
                        }
                        elsif ($strand eq "-" and $tblock_end >= $coding_end and $tblock_start <= $coding_end) {
                            my $distance_in_block_to_coding_start = $coding_end - $tblock_start;
                            $SE_coding_start = $SE_coding_start + $distance_in_block_to_coding_start;
                        }
                        if ($strand eq "+" and $tblock_end == $SE_donor) {
                            my $running_sum = 0;
                            foreach my $j(@SE_block_sizes) {
                                $running_sum += $j;
                            }
                            $SE_donor_junction = $running_sum;
                        }
                        elsif ($strand eq "-" and $tblock_start == $SE_acceptor) {
                            my $running_sum = 0;
                            for (my $j = 0; $j < (scalar(@SE_block_sizes) - 1); $j++) {
                                $running_sum = $running_sum + $SE_block_sizes[$j];
                            }
                            $SE_donor_junction = $running_sum+1;
                        }
                        else {
                        }
                    }
                    else {
                        push(@skipped_exons_size_array, $block_sizes[$i]);
                        push(@skipped_exons_start_array, $block_starts[$i]);
                    }

                    if ($coding_start > $tblock_start and $coding_start < $tblock_end) {
                        if ($coding_end > $tblock_end) {
                            my $first_coding_size = ($tblock_end - $coding_start) + 1;
                            my $coding_block_sequence = substr $gelem, $coding_start, $first_coding_size;
                            $coding_region_sequence = $coding_region_sequence.$coding_block_sequence;
                        }
                        else {
                            my $within_exon_coding_size = ($coding_end - $coding_start);
                            my $coding_block_sequence = substr $gelem, $coding_start, $within_exon_coding_size;
                            $coding_region_sequence = $coding_region_sequence.$coding_block_sequence;
                        }
                    }
                    elsif ($tblock_start < ($coding_end-1) and $tblock_end > ($coding_end-1)) {
                        my $last_coding_size = ($coding_end - $tblock_start);
                        my $coding_block_sequence = substr $gelem, $tblock_start, $last_coding_size;
                        $coding_region_sequence = $coding_region_sequence.$coding_block_sequence;
                    }
                    elsif ($coding_start < $tblock_start and ($coding_end-1) > $tblock_end) {
                        $coding_region_sequence = $coding_region_sequence.$block_sequence;
                    }
                }
                my $SE_junction_donors = 0;
                my $cumulative_length = 0;
                foreach my $block(@SE_block_sizes) {
                    $cumulative_length += $block;
                    $SE_junction_donors = $SE_junction_donors."\,".$cumulative_length;
                }
                my $skipped_exons_length = 0;
                my $frameshift_parameter;

                foreach my $skipped_block(@skipped_exons_size_array) {
                    $skipped_exons_length += $skipped_block;
                }

                if (int($skipped_exons_length/3) != $skipped_exons_length/3) {
                    $frameshift_parameter = "FRAMESHIFT";
                }
                else {
                    $frameshift_parameter = "INFRAME";
                    for (my $c = 0; $c < scalar(@skipped_exons_size_array); $c++) {
                        my $block_sequence = substr $gelem, ($trans_start + $skipped_exons_start_array[$c]), ($skipped_exons_size_array[$c]);
                        $inframe_skipped_exons_sequence = $inframe_skipped_exons_sequence.$block_sequence;
                    }
                }
                if ($strand eq "+") {
                    print OUT1 "\>", $gene_name, "\n", $transcript_sequence, "\n";
                    
                    if ($coding_start == $coding_end) {
                        print OUT2 "\>", $gene_name, "\nNoncoding RNA\n";
                    }
                    elsif ($coding_start != $coding_end) {
                        print OUT2 "\>", $gene_name, "\n", $coding_region_sequence, "\n";
                        next if ($SE_donor_junction == 0);
                        if ($IncDiff < 0) { 
                            print OUT3 "\>", $gene_name, "\;", $SE_junction_donors, "\;", $SE_coding_start + 1, "\;". $SE_donor_junction, "\;", $frameshift_parameter, "\n", $SE_transcript_sequence, "\n";#
                            if ($inframe_skipped_exons_sequence) {
                                my $frame_offset = offset_calculator($SE_donor_junction - $SE_coding_start);
                                print OUT5 "\>", $gene_name, "\;", $SE_junction_donors, "\;", $SE_coding_start + 1, "\;". $SE_donor_junction, "\;", $frameshift_parameter, "\n", $inframe_skipped_exons_sequence, "\t", $frame_offset, "\n";#
                            }
                        }
                        elsif ($IncDiff > 0) {
                            print OUT4 "\>", $gene_name, "\;", $SE_junction_donors, "\;", $SE_coding_start + 1, "\;". $SE_donor_junction, "\;", $frameshift_parameter, "\n", $SE_transcript_sequence, "\n";#
                            if ($inframe_skipped_exons_sequence) {
                                my $frame_offset = offset_calculator($SE_donor_junction - $SE_coding_start);
                                print OUT6 "\>", $gene_name, "\;", $SE_junction_donors, "\;", $SE_coding_start + 1, "\;". $SE_donor_junction, "\;", $frameshift_parameter, "\n", $inframe_skipped_exons_sequence, "\t", $frame_offset, "\n";#
                            }
                        }
                    }
                }
                elsif ($strand eq "-") {
                    my $reverse_transcript_sequence = scalar reverse($transcript_sequence);
                    $reverse_transcript_sequence =~ tr/ACGTacgt/TGCATGCA/;
                    my $reverse_SE_transcript_sequence = scalar reverse($SE_transcript_sequence);#
                    $reverse_SE_transcript_sequence =~ tr/ACGTacgt/TGCATGCA/;#
                    my $reverse_inframe_skipped_exons_sequence;
                    if ($inframe_skipped_exons_sequence) {
                        $reverse_inframe_skipped_exons_sequence = scalar reverse($inframe_skipped_exons_sequence);
                        $reverse_inframe_skipped_exons_sequence =~ tr/ACGTacgt/TGCATGCA/;#
                    }
                    print OUT1 "\>", $gene_name, "\n", $reverse_transcript_sequence, "\n";
                    
                    if ($coding_start == $coding_end) {
                        print OUT2 "\>", $gene_name, "\nNoncoding RNA\n";
                    }
                    elsif ($coding_start != $coding_end) {
                        my $reverse_coding_region_sequence = scalar reverse ($coding_region_sequence);
                        $reverse_coding_region_sequence =~ tr/ACGTacgt/TGCATGCA/;
                        print OUT2 "\>", $gene_name, "\n", $reverse_coding_region_sequence, "\n";
                        next if ($SE_donor_junction == 0);

                        my $reverse_SE_coding_start = $cumulative_length - $SE_coding_start;
                        my $reverse_SE_donor_junction = $cumulative_length - $SE_donor_junction + 1;
                        my @split_SE_junction_donors = split("\,", $SE_junction_donors);
                        my @subtracted_SE_junction_donors = ();
                        for (my $k = 0; $k < scalar(@split_SE_junction_donors); $k++) {
                            my $reverse_length = $cumulative_length - $split_SE_junction_donors[$k];
                            push(@subtracted_SE_junction_donors, $reverse_length);
                        }
                        my @reverse_subtracted_SE_junction_donors = reverse @subtracted_SE_junction_donors;
                        my $joined_reverse_subtracted_SE_junction_donors = join("\,", @reverse_subtracted_SE_junction_donors);

                        if ($IncDiff < 0) { 
                            print OUT3 "\>", $gene_name, "\;", $joined_reverse_subtracted_SE_junction_donors, "\;", $reverse_SE_coding_start + 1, "\;". $reverse_SE_donor_junction,   "\;", $frameshift_parameter, "\n", $reverse_SE_transcript_sequence,"\n";
                            if ($reverse_inframe_skipped_exons_sequence) {
                                my $frame_offset = offset_calculator($reverse_SE_donor_junction - $reverse_SE_coding_start);
                                print OUT5 "\>", $gene_name, "\;", $joined_reverse_subtracted_SE_junction_donors, "\;", $reverse_SE_coding_start + 1, "\;". $reverse_SE_donor_junction,   "\;", $frameshift_parameter, "\n", $reverse_inframe_skipped_exons_sequence, "\t", $frame_offset, "\n";
                            }
                        }
                        elsif ($IncDiff > 0) {
                            print OUT4 "\>", $gene_name, "\;", $joined_reverse_subtracted_SE_junction_donors, "\;", $reverse_SE_coding_start + 1, "\;". $reverse_SE_donor_junction,  "\;", $frameshift_parameter, "\n", $reverse_SE_transcript_sequence, "\n";
                            if ($reverse_inframe_skipped_exons_sequence) {
                                my $frame_offset = offset_calculator($reverse_SE_donor_junction - $reverse_SE_coding_start);
                                print OUT6 "\>", $gene_name, "\;", $joined_reverse_subtracted_SE_junction_donors, "\;", $reverse_SE_coding_start + 1, "\;". $reverse_SE_donor_junction,   "\;", $frameshift_parameter, "\n", $reverse_inframe_skipped_exons_sequence, "\t", $frame_offset, "\n";
                            }
                        }  
                    }
                }
            }
        }
        
    }
    close(INF); 
    close(OUT1);
    close(OUT2);
    close(OUT3);
    close(OUT4);
    close(OUT5);
    close(OUT6);
}

sub translate {
    print "\nTranslating skipped exon isoforms and NMD determinations...\n";
    my @files = ($outdir."\/1_nucleotide_sequences/SE_transcripts_pos_IncDiff.fa", $outdir."\/1_nucleotide_sequences/SE_transcripts_neg_IncDiff.fa");
    foreach my $file(@files) {
        my $file_ID;
        if (basename($file) =~ m/pos_IncDiff/g) {
            $file_ID = "pos_IncDiff";
        }
        elsif (basename($file) =~ m/neg_IncDiff/g) {
            $file_ID = "neg_IncDiff";
        }
        my $file_basename = basename($file);
        print "\t", $file_basename, "\n\t\tSE Isoform protein sequences...\n\t\tSE Neopeptides...\n";
        open(INF, "<$file") or die "couldn't open bed file";
        open(OUT1, ">$outdir/2_SE_isoforms_all/SE_isoforms_all_$file_ID.fa") or die "couldn't open translated fasta output fasta file";
        open(OUT1b, ">$outdir/2_SE_isoforms_all/temp/SE_isoforms_all_$file_ID.temp") or die "couldn't open translated temporary output file";
        open(OUT2, ">$outdir/3_frameshifted_isoforms/temp/Frameshifted_SE_isoforms_$file_ID.fa") or die "couldn't open translated frameshifted output fasta file";
        open(OUT2b, ">$outdir/3_frameshifted_isoforms/temp/Frameshifted_SE_isoforms_$file_ID.temp") or die "couldn't open translated frameshifted temporary output fasta file";
        open(OUT3,  ">$outdir/4_neopeptides/temp/SE_NEOPEPTIDES_$file_ID.fa") or die "couldn't open neopeptide output file";
        open(OUT3b,  ">$outdir/4_neopeptides/temp/SE_NEOPEPTIDES_$file_ID.temp") or die "couldn't open neopeptide output file";

        my $gene_info;
        my $count = 0;
        while(my $line = <INF>) {
            chomp($line);
            if ($line =~ m/^\>/) {
                $gene_info = $line;
                $count++;
            }
            elsif ($count == 1) {
                $count = 0;
                my @split_gene_info = split("\;", $gene_info);
                my @split_splice_donors = split("\,", $split_gene_info[1]);
                my @split_gene_ID_SE_junction_info = split("\,", $split_gene_info[0]);

                my $transcript_length = length($line);
                my $coding_start = $split_gene_info[2];
                my $skipped_exon_donor_start = $split_gene_info[3];
                my $frameshift_call = $split_gene_info[4];
                my $frameshift_call_relative_to_coding_start = $skipped_exon_donor_start - $coding_start - 1;
                my $frameshift_call_relative_to_coding_start_in_codons = int($frameshift_call_relative_to_coding_start/3);
                my $protein_sequence = "";
                for (my $a = ($coding_start - 1); $a < ($transcript_length - 4); $a += 3) {
                    my $codon = substr $line, $a, 3;
                    my $amino_acid = codon2aa($codon);
                    if ($amino_acid eq "_") {
                        last;
                    }
                    else {
                        $protein_sequence = $protein_sequence.$amino_acid;
                    }
                }
                my $protein_sequence_length = length($protein_sequence);

                next if ($protein_sequence_length - $frameshift_call_relative_to_coding_start_in_codons < 0);
                my $stop_codon_position_on_transcript = $coding_start + ($protein_sequence_length*3);
                my $neopeptide = substr $protein_sequence, $frameshift_call_relative_to_coding_start_in_codons, ($protein_sequence_length - $frameshift_call_relative_to_coding_start_in_codons + 1);

                if ($frameshift_call eq "INFRAME") {
                    print OUT1 $gene_info, "\;", $stop_codon_position_on_transcript, "\;no_NMD\;", $frameshift_call_relative_to_coding_start_in_codons, "\n", $protein_sequence, "\n";
                    print OUT1b $gene_info, "\;", $stop_codon_position_on_transcript, "\;no_NMD\;", $frameshift_call_relative_to_coding_start_in_codons, "\tno_NMD\t", $protein_sequence, "\n";
                }
                elsif ($frameshift_call eq "FRAMESHIFT") {
                    my $NMD;
                    if ($split_splice_donors[-2] - $stop_codon_position_on_transcript >= 50) {
                        $NMD = "NMD";
                    }
                    else {
                        $NMD = "no_NMD";
                    }

                    print OUT1 $gene_info, "\;", $stop_codon_position_on_transcript, "\;", $NMD, "\;", $frameshift_call_relative_to_coding_start_in_codons, "\n", $protein_sequence, "\n";
                    print OUT1b $gene_info, "\;", $stop_codon_position_on_transcript, "\;", $NMD, "\;", $frameshift_call_relative_to_coding_start_in_codons, "\t", $NMD, "\t",$protein_sequence, "\n";
                    print OUT2 $gene_info, "\;", $stop_codon_position_on_transcript, "\;", $NMD, "\;", $frameshift_call_relative_to_coding_start_in_codons, "\n", $protein_sequence, "\n";
                    print OUT2b $gene_info, "\;", $stop_codon_position_on_transcript, "\;", $NMD, "\;", $frameshift_call_relative_to_coding_start_in_codons, "\t", $protein_sequence, "\n";
                    print OUT3 $gene_info, "\;", $stop_codon_position_on_transcript, "\;", $NMD, "\;", $frameshift_call_relative_to_coding_start_in_codons, "\n", $neopeptide, "\n";
                    print OUT3b $gene_info, "\;", $stop_codon_position_on_transcript, "\;", $NMD, "\;", $frameshift_call_relative_to_coding_start_in_codons, "\t", $neopeptide, "\n";
                }
            }
        }
        print OUT1b "e,n,d,x_y\tend\tend";
        print OUT2b "e,n,d,x_y\tend\tend";
        print OUT3b "e,n,d,x_y\tend\tend";
        close(INF);
        close(OUT1);
        close(OUT1b);
        close(OUT2);
        close(OUT2b);
        close(OUT3);
        close(OUT3b);
        close(OUT4);

        my $presorted_TRANSLATED_Frameshifted_file = $outdir."\/3_frameshifted_isoforms/temp/Frameshifted_SE_isoforms_".$file_ID.".temp";
        my $sorted_TRANSLATED_Frameshifted_file = $outdir."\/3_frameshifted_isoforms/temp/Frameshifted_SE_isoforms_".$file_ID.".sorted";
        `sort -k2,2 $presorted_TRANSLATED_Frameshifted_file > $sorted_TRANSLATED_Frameshifted_file`;

        open(INF, "<$sorted_TRANSLATED_Frameshifted_file") or die "couldn't open translated frameshifted input file";
        open(OUT, ">$outdir/3_frameshifted_isoforms/Frameshifted_SE_isoforms_deduped_$file_ID.fa") or die "couldn't open deduplicated translated frameshifted output file";
        my $previous_frameshifted_sequence;
        my @ID_frameshift_array = ();
        my @amino_acid_divergence_junction_array;
        while(my $line = <INF>) {
            chomp($line);
            my @split_line = split("\t", $line);
            my @split_ID = split("\;", $split_line[0]);
            my $transcript_ID = $split_ID[0];
            my $amino_acid_divergence = $split_ID[-1];
            
            if ($. == 1) {
                $previous_frameshifted_sequence = $split_line[1];
                push(@ID_frameshift_array, $transcript_ID);
                push(@amino_acid_divergence_junction_array, $amino_acid_divergence);
            }
            elsif ($split_line[1] eq $previous_frameshifted_sequence) {
                push(@ID_frameshift_array, $transcript_ID);
                push(@amino_acid_divergence_junction_array, $amino_acid_divergence);
            }
            else {
                my $new_ID = join("\,", @ID_frameshift_array);
                $new_ID =~ s/>//g;
                my $new_amino_acid_divergence =  join("\,", @amino_acid_divergence_junction_array);

                print OUT "\>", $new_ID, "\t", $new_amino_acid_divergence, "\n", $previous_frameshifted_sequence, "\n";

                $previous_frameshifted_sequence = $split_line[1];
                @ID_frameshift_array = ();
                push(@ID_frameshift_array, $transcript_ID);
                @amino_acid_divergence_junction_array = ();
                push(@amino_acid_divergence_junction_array, $amino_acid_divergence);
            }
        }
        close(INF);
        close(OUT); 
        my $presorted_neopeptides_file = $outdir."\/4_neopeptides/temp/SE_NEOPEPTIDES_".$file_ID.".temp";
        my $sorted_neopeptides_file = $outdir."\/4_neopeptides/temp/SE_NEOPEPTIDES.$file_ID.sorted";
        `sort -k2,2 $presorted_neopeptides_file > $sorted_neopeptides_file`;
        open(INF, "<$sorted_neopeptides_file") or die "couldn't open neopeptides input file";
        open(OUT, ">$outdir/4_neopeptides/SE_NEOPEPTIDES.deduped.$file_ID.fa") or die "couldn't open deduplicated neopeptides output file";
        my $previous_neopeptide_sequence;
        my @ID_neopeptide_array = ();
        my $counter = 0;
        while(my $line = <INF>) {
            chomp($line);
            my @split_line = split("\t", $line);
            my @split_ID = split("\;", $split_line[0]);
            my $transcript_ID = $split_ID[0];
            next if (not defined $split_line[1]);
            if ($counter == 0) {
                $previous_neopeptide_sequence = $split_line[1];
                push(@ID_neopeptide_array, $transcript_ID);
                $counter++;
            }
            elsif ($split_line[1] eq $previous_neopeptide_sequence) {
                push(@ID_neopeptide_array, $transcript_ID);
            }
            else {
                my $new_ID = join("\,", @ID_neopeptide_array);
                $new_ID =~ s/>//g;

                print OUT "\>", $new_ID, "\n", $previous_neopeptide_sequence, "\n";

                $previous_neopeptide_sequence = $split_line[1];
                @ID_neopeptide_array = ();
                push(@ID_neopeptide_array, $transcript_ID);
            }
        }
        close(INF);
        close(OUT); 

        my $presorted_file = $outdir."\/2_SE_isoforms_all/temp/SE_isoforms_all_$file_ID.temp";
        my $sorted_file = $outdir."\/2_SE_isoforms_all/temp/SE_isoforms_all_$file_ID.sorted";
        `sort -k3,3 -k2,2 $presorted_file > $sorted_file`;

        my @file_lines = ();
        open(INF, "<$sorted_file") or die "couldn't open translated temp input file";
        while(my $line = <INF>) {
            chomp($line);
            push(@file_lines, $line);
        }
        close(INF);

        open(OUT, ">$outdir/6_NMD/SE_protein_isoforms_$file_ID.deduped.fa") or die "couldn't open translated decuplicated output file";
        my $foreach_count = 0;
        my $file_lines_length = @file_lines;
        my $previous_sequence = "";
        my @ID_array = ();
        my $previous_NMD = "";
        foreach my $line(@file_lines) {
            my @split_line = split("\t", $line);
            my @split_ID = split("\;", $split_line[0]);
            my $transcript_ID = $split_ID[0];
            if ($foreach_count == 0) {
                $previous_sequence = $split_line[2];
                push(@ID_array, $transcript_ID);
                $previous_NMD = $split_line[1];
            }
            elsif ($split_line[2] eq $previous_sequence and $split_line[1] eq $previous_NMD) {
                push(@ID_array, $transcript_ID);
                if ($foreach_count == $file_lines_length-1) {
                    my $new_ID = join("\,", @ID_array);
                    $new_ID =~ s/>//g;
                    print OUT "\>", $new_ID, "\t", $previous_NMD, "\n", $previous_sequence, "\n";
                }
            }
            else {
                my $new_ID = join("\,", @ID_array);
                $new_ID =~ s/>//g;
                print OUT "\>", $new_ID, "\t", $previous_NMD, "\n", $previous_sequence, "\n";

                $previous_sequence = $split_line[2];
                @ID_array = ();
                push(@ID_array, $transcript_ID);
                $previous_NMD = $split_line[1];
            }
            $foreach_count++;
        }
        close(OUT);

        my $sorted_for_NMD = $outdir."\/6_NMD/$file_basename.TRANSLATED.sorted_for_NMD";
        `sort -k1,1 -k2,2 $presorted_file > $sorted_for_NMD`;
        my @NMD_file_lines = ();
        open(INF, "<$sorted_for_NMD") or die "couldn't open translated temp input file";
        while(my $line = <INF>) {
            chomp($line);
            push(@NMD_file_lines, $line);
        }
        close(INF);
        `rm $sorted_for_NMD`;

        print "\t\tNMD calculations...\n\n";
        open(OUT, ">$outdir/6_NMD/SE_NMD_statistics.$file_ID.txt") or die "couldn't open NMD counts output file";
        open(OUTa, ">$outdir/6_NMD/SE_NMD_lists.bed/SE_NMD.$file_ID.bed") or die "couldn't open NMD counts output file";
        open(OUTb, ">$outdir/6_NMD/SE_NMD_lists.bed/SE_no_NMD.$file_ID.bed") or die "couldn't open NMD counts output file";
        open(OUTnote, ">$outdir/6_NMD/SE_NMD_lists.bed/1_note_about_loading_bed_files_on_IGV.txt") or die "couldn't open NMD counts output file";
        print OUTnote "NOTE: Header line of bed files must be deleted to load on IGV.\n";
        close(OUTnote);
        print OUTa "chr\tSE_donor\tSE_acceptor\tNMD_isoforms\tno_NMD_isoforms\tstrand\n";
        print OUTb "chr\tSE_donor\tSE_acceptor\tno_NMD_isoforms\tno_NMD_isoforms\tstrand\n";
        $foreach_count = 0;
        my $NMD = 0;
        my $no_NMD = 0;
        my $gene_NMD = 0;
        my $gene_no_NMD = 0;
        my $previous_gene;
        my @NMD_genes = ();
        my @no_NMD_genes = ();
        my @NMD_isoforms = ();
        my @no_NMD_isoforms = ();
        my $previous_chr_junction_coords = "";
        my $previous_strand;
        foreach my $elem(@NMD_file_lines) {
            my @split_elem = split("\t", $elem);
            my @split_info = split("\;", $split_elem[0]);       
            my @split_info_0 = split("\,", $split_info[0]);
            my $gene_isoform_ID = $split_info_0[3];
            my @split_gene_isoform_ID = split("_", $gene_isoform_ID);
            my $isoform_ID = $split_gene_isoform_ID[0];
            my $chr_junction_coords = join("\t", @split_info_0[0..2]);
            $chr_junction_coords =~ s/^\>//;
            my $strand = $split_info_0[5];
            if ($foreach_count == 0) {
                $previous_chr_junction_coords = $chr_junction_coords;
                $previous_strand = $strand;
                $previous_gene = $split_gene_isoform_ID[1];
                if($split_elem[1] eq "NMD") {
                    push(@NMD_isoforms, $isoform_ID);
                    $gene_NMD++;
                }
                elsif($split_elem[1] eq "no_NMD") {
                    push(@no_NMD_isoforms, $isoform_ID);
                    $gene_no_NMD++;
                }
            }
            elsif ($foreach_count == scalar(@NMD_file_lines)-1) {
                if ($gene_NMD > 0 and $gene_no_NMD > 0) {
                    print OUTa $previous_chr_junction_coords, "\t", $previous_gene, "_", join("_", @NMD_isoforms), "\t", join("_", @no_NMD_isoforms), "\t", $previous_strand, "\n";
                    my @split_gene_ID = split("_", $NMD_isoforms[0]);
                    push(@NMD_genes, $previous_gene);
                    $NMD++;
                }
                elsif ($gene_NMD == 0 and $gene_no_NMD > 0){
                    print OUTb $previous_chr_junction_coords, "\t", $previous_gene, "_\.\t", join("_", @no_NMD_isoforms), "\t", $previous_strand, "\n";
                    my @split_gene_ID = split("_", $no_NMD_isoforms[0]);
                    push(@no_NMD_genes, $previous_gene);
                    $no_NMD++;
                }
                elsif ($gene_NMD > 0 and $gene_no_NMD == 0){
                    print OUTa $previous_chr_junction_coords, "\t", $previous_gene, "_", join("_", @NMD_isoforms), "\t\.\t", $previous_strand, "\n";
                    my @split_gene_ID = split("_", $NMD_isoforms[0]);
                    push(@NMD_genes, $previous_gene);
                    $NMD++;
                }
            }
            elsif ($chr_junction_coords eq $previous_chr_junction_coords) {
                if($split_elem[1] eq "NMD") {
                    push(@NMD_isoforms, $isoform_ID);
                    $gene_NMD++;
                }
                elsif($split_elem[1] eq "no_NMD") {
                    push(@no_NMD_isoforms, $isoform_ID);
                    $gene_no_NMD++;
                }
                
            }
            else {
                if ($gene_NMD > 0 and $gene_no_NMD > 0) {
                    print OUTa $previous_chr_junction_coords, "\t", $previous_gene, "_", join("_", @NMD_isoforms), "\t", join("_", @no_NMD_isoforms), "\t", $previous_strand, "\n";
                    my @split_gene_ID = split("_", $NMD_isoforms[0]);
                    push(@NMD_genes, $previous_gene);
                    $NMD++;
                }
                elsif ($gene_NMD == 0 and $gene_no_NMD > 0){
                    print OUTb $previous_chr_junction_coords, "\t", $previous_gene, "\.\t", join("_", @no_NMD_isoforms), "\t", $previous_strand, "\n";
                    my @split_gene_ID = split("_", $no_NMD_isoforms[0]);
                    push(@no_NMD_genes, $previous_gene);
                    $no_NMD++;
                }
                elsif ($gene_NMD > 0 and $gene_no_NMD == 0){
                    print OUTa $previous_chr_junction_coords, "\t", $previous_gene, "_", join("_", @NMD_isoforms), "\t\.\t", $previous_strand, "\n";
                    my @split_gene_ID = split("_", $NMD_isoforms[0]);
                    push(@NMD_genes, $previous_gene);
                    $NMD++;
                }
                $previous_chr_junction_coords = $chr_junction_coords;
                $previous_strand = $strand;
                @NMD_isoforms = ();
                @no_NMD_isoforms = ();
                $gene_NMD = 0;
                $gene_no_NMD = 0;
                $previous_gene = $split_gene_isoform_ID[1];
                if($split_elem[1] eq "NMD") {
                    push(@NMD_isoforms, $isoform_ID);
                    $gene_NMD++;
                }
                elsif($split_elem[1] eq "no_NMD") {
                    push(@no_NMD_isoforms, $isoform_ID);
                    $gene_no_NMD++;
                }
            }
            $foreach_count++;
        }
        my $tot = $NMD+$no_NMD;
        if ($tot == 0) {
            $tot = 1;
        }
        print OUT "Number predicted to undergo NMD: ", $NMD, "\n\tFraction predicted to undergo NMD: ", $NMD/($tot), "\n\nNumber predicted to not undergo NMD: ", $no_NMD, "\n\tFraction predicted to not undergo NMD: ", $no_NMD/($tot), "\n";
        close(OUT); 
        close(OUTa);
        close(OUTb);

        open(OUTc, ">$outdir/6_NMD/gene_lists/SE_NMD.$file_ID.txt") or die "couldn't open NMD counts output file";
        open(OUTd, ">$outdir/6_NMD/gene_lists/SE_no_NMD.$file_ID.txt") or die "couldn't open NMD counts output file";
        my @unique_NMD_genes = unique(@NMD_genes);
        my @unique_no_NMD_genes = unique(@no_NMD_genes);
        print OUTc join("\n", @unique_NMD_genes);
        print OUTd join("\n", @unique_no_NMD_genes);
        close(OUTc);
        close(OUTd);
        
    }
}

sub translate_skipped_exons {
    print "Translating in-frame skipped exons...\n";
    my @files = ($outdir."\/1_nucleotide_sequences/SKIPPED_EXONs_pos_IncDiff.fa", $outdir."\/1_nucleotide_sequences/SKIPPED_EXONs_neg_IncDiff.fa");
    foreach my $file(@files) {
        my $file_ID;
        if (basename($file) =~ m/pos_IncDiff/g) {
            $file_ID = "pos_IncDiff";
        }
        elsif (basename($file) =~ m/neg_IncDiff/g) {
            $file_ID = "neg_IncDiff";
        }

        print "\t", basename($file), "\n";
        open(INF, "<$file") or die "couldn't open bed file";
        open(OUT, ">$outdir/5_skipped_exon_protein_sequences/temp/SKIPPED_EXONs_$file_ID.fa") or die "couldn't open translated fasta output file";
        open(OUT2, ">$outdir/5_skipped_exon_protein_sequences/temp/SKIPPED_EXONs_$file_ID.temp") or die "couldn't open translated temporary fasta output file";
        my $gene_info;
        my $count = 0;
        while(my $line = <INF>) {
            chomp($line);
            if ($line =~ m/^\>/) {
                $gene_info = $line;
                $count++;
            }
            elsif ($count == 1) {
                my @split_line = split("\t", $line);
                my $sequence = $split_line[0];
                my $offset = $split_line[1];
                my $transcript_length = length($sequence);
                my $protein_sequence = "";
                my $max;
                if ($offset == 3) {
                    $max = $transcript_length;
                }
                else {
                    $max = $transcript_length - 2;
                }
                for (my $a = (3 - $offset); $a < $max; $a += 3) {
                    my $codon = substr $line, $a, 3;
                    my $amino_acid = codon2aa($codon);
                    if ($amino_acid eq "_") {
                        last;
                    }
                    else {
                        $protein_sequence = $protein_sequence.$amino_acid;
                    }
                }
                my @split_gene_info = split("\;", $gene_info);
                my $transcript_ID = $split_gene_info[0];
                print OUT $transcript_ID, "\n", $protein_sequence, "\n";
                print OUT2 $transcript_ID, "\t", $protein_sequence, "\n";
                $count = 0;
            }  
        }
        print OUT2 "null\tnull";
        close(INF);
        close(OUT);
        close(OUT2);

        my $presorted_file = $outdir."\/5_skipped_exon_protein_sequences\/temp\/SKIPPED_EXONs_$file_ID.temp";
        my $sorted_file = $outdir."\/5_skipped_exon_protein_sequences\/temp\/SKIPPED_EXONs_$file_ID.sorted";
        `sort -k2,2 $presorted_file > $sorted_file`;
        open(INF, "<$sorted_file") or die "couldn't open sorted input file";
        open(OUT1, ">$outdir/5_skipped_exon_protein_sequences/SKIPPED_EXONs_deduped_$file_ID.sorted.fa") or die "couldn't open deduped output file";
        open(OUT2, ">$outdir/5_skipped_exon_protein_sequences/SKIPPED_EXONs_deduped_BATCH_$file_ID.sorted.fa") or die "couldn't open deduped batch submission output file";

        my $previous_sequence = "null";
        my @ID_array = ();
        while(my $line = <INF>) {
            chomp($line);
            my @split_line = split("\t", $line);
            my $transcript_ID = $split_line[0];

            next if (!defined($split_line[1]));
            if ($. == 1 || scalar(@ID_array) == 0) {
                $previous_sequence = $split_line[1];
                push(@ID_array, $transcript_ID);
            }
            elsif ($split_line[1] eq $previous_sequence) {
                push(@ID_array, $split_line[0]);
            }
            else {
                my $new_ID = join("\,", @ID_array);
                $new_ID =~ s/>//g;
                my @split_new_ID = split(",", $new_ID);
                my @split_gene = split("_", $split_new_ID[3]);
                my $batch_new_ID = $split_gene[1]." ".$split_new_ID[0]."\:".$split_new_ID[1]."\-".$split_new_ID[2]." ".join(" ", @split_new_ID[4..5]);
                print OUT1 "\>", $new_ID, "\n", $previous_sequence, "\n";
                print OUT2 "\>", $batch_new_ID, "\n", $previous_sequence, "\n";

                $previous_sequence = $split_line[1];
                @ID_array = ();
                push(@ID_array, $transcript_ID);
            }
        }
        close(INF);
        close(OUT1);
        close(OUT2);    
        #`rm $presorted_file`;
        #`rm $sorted_file`;
    }
}
sub rm_temp {
    `rm -r $outdir"\/2_SE_isoforms_all/temp"`;
    `rm -r $outdir"\/3_frameshifted_isoforms/temp"`;
    `rm -r $outdir"\/4_neopeptides/temp"`;
    `rm -r $outdir"\/5_skipped_exon_protein_sequences/temp"`;
}

# Accessory subroutines
sub loadcode {
    open(INF, "<$script_dirname/genetic_code.tsv") or die "couldn't open genetic code file";
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        $genetic_code{$split_line[0]} = $split_line[1];
    }
    close(INF);
}

sub offset_calculator {
    my $x = ($_[0]/3) - int($_[0]/3);
    if ($x > 0.2 and $x < 0.4) {
        $x = 1;
    }
    elsif ($x > 0.5 and $x < 0.7) {
        $x = 2;
    }
    elsif ($x == 0) {
        $x = 3;
    }
    return $x;
}

sub codon2aa {
    my($codon) = @_;
    $codon = uc $codon;
    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }
    else{
        die "Bad codon '$codon'!!\n";
    }
}

sub unique {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

options;
qc;
print "\n\n***********************\nRunning SETranslateNMD...\n";
makedirectories;
print_input_parameters;
process_skipped_exon_array;
process_annotation_file;
link_annotation_to_skipped_exons;
genome_fasta_file_unwrapper;
bed2fasta;
loadcode;
translate;
translate_skipped_exons;
rm_temp;
print "\nSETranslateNMD DONE!\n***********************\n\n";