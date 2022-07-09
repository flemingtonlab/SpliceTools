use warnings;
use strict;  
use File::Basename;

my $retained_intron_file;
my $annotation_file;
my $genome_fasta_file;
my $expression_file;
my $TPM1;
my $TPM2;
my $cond1_ct;
my $cond2_ct;
my $FDR;

my $dirname = dirname(__FILE__);

sub program_info {
    print "\n\tRIMedley.pl will automatically run the following on an rMATS RI JCEC file:\n\t\t- RIFractionExpressed.pl\n\t\t- RIIntronExonSizes.pl\n\t\t- RISpliceSiteScoring.pl\n\n\tUsage: perl RIMedley.pl [OPTIONS] -r <retained intron file (rMATS JCEC)> -a <bed12 annotation file> -g <genome fasta file> -e <expression file> -TPM <min TPMs condition 1,min TPMs condition 2> -SN <sample number condition 1,sample number condition2> -f <FDR>\n\n\tRequired:\n\t\t-r <retained intron file>\n\t\t-a <bed12 annotation file>\n\t\t-g <genome fasta file>\n\t\t-e <expression file>\n\t\t-TPM <min TPMs condition 1,min TPMs condition 2>\n\t\t-SN <sample number condition 1,sample number condition 2>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl RIMedley.pl -r PATH/RIfile.txt -a PATH/bed12_annotation.bed -g PATH/genome.fa -e PATH/expression.tsv -TPM 2,2 -SN 3,3 -f 0.05\n\n";
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
        elsif ($ARGV[$i] eq "\-g") {
            $genome_fasta_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-a") {
            $annotation_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-f") {
            $FDR = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-e") {
            $expression_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-TPM") {
            my $TPMs = $ARGV[$i+1];
            ($TPM1, $TPM2) = split("\,", $TPMs);
        }
        elsif ($ARGV[$i] eq "\-SN") {
            my $sample_counts = $ARGV[$i+1];
            ($cond1_ct, $cond2_ct) = split("\,", $sample_counts);
        }
        elsif ($ARGV[$i] eq "\-h") {
            program_info;
            exit;
        }
    }
}

sub qc {
    if ($FDR !~ m/^\d/) {
        print "\nFDR is not numeric!\n";
        program_info;
        exit;
    } 
    elsif ($retained_intron_file eq "") {
        print "\nRetained intron file not defined!\n";
        program_info;
        exit;
    }
    elsif ($genome_fasta_file eq "") {
        print "\nGenome fasta file not defined!\n";
        program_info;
        exit;
    }
    elsif ($expression_file eq "") {
        print "\nExpression file not defined!\n";
        program_info;
        exit;
    }
    elsif ($FDR eq "") {
        print "\nFDR not defined!\n";
        program_info;
        exit;
    }
    elsif ($TPM1 eq "") {
        print "\nMinimum TPMs for condition 1 not defined!\n     Is -TPM value in correct format?";
        program_info;
        exit;
    }
    elsif ($TPM2 eq "") {
        print "\nMinimum TPMs for condition 2 not defined!\n     Is -TPM value in correct format?";
        program_info;
        exit;
    }
    elsif ($cond1_ct eq "") {
        print "\nCondition 1 sample count not defined!\n     Is -SN value in correct format?";
        program_info;
        exit;
    }
    elsif ($cond2_ct eq "") {
        print "\nCondition 2 sample count not defined!\n     Is -SN value in correct format?";
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

sub print_input_parameters {
    print "\n\nINPUT:\nRetained intron file: ", basename($retained_intron_file), "\nAnnotation file: ", basename($annotation_file), "\nGenome fasta file: ", basename($genome_fasta_file),  "\nExpression file: ", basename($expression_file), "\nMin TPMs, condition 1 = ", $TPM1, "\nMin TPMs, condition 2 = ", $TPM2, "\nSample number condition 1 = ", $cond1_ct, "\nSample number condition 2 = ", $cond2_ct, "\nSE FDR = ", $FDR, "\n";
}

sub RISpliceSiteScoring {
    my $RISpliceSiteScoring = $dirname."\/RISpliceSiteScoring.pl";
    system ("perl $RISpliceSiteScoring -r $retained_intron_file -g $genome_fasta_file -a $annotation_file -f $FDR");
}

sub RIIntronExonSizes {
    my $RIIntronExonSizes = $dirname."\/RIIntronExonSizes.pl";
    system ("perl $RIIntronExonSizes -r $retained_intron_file -a $annotation_file -f $FDR");
}

sub RIFractionExpressed {
    my $RIFractionExpressed = $dirname."\/RIFractionExpressed.pl";
    system ("perl $RIFractionExpressed -r $retained_intron_file -e $expression_file -SN $cond1_ct,$cond2_ct -TPM $TPM1,$TPM2 -f $FDR");
}


options;
qc;
print "\n\nRunning RIKitchenSink...";
print_input_parameters;
RISpliceSiteScoring;
RIIntronExonSizes;
RIFractionExpressed;
