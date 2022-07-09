use warnings;
use strict; 
use File::Basename;  

my $retained_intron_directory;
my $annotation_file;
my $expression_file;
my $TPM1;
my $TPM2;
my $cond1_ct;
my $cond2_ct;
my $FDR;
my @RI_files;
my $program_path_name = dirname(__FILE__)."\/RIFractionExpressed.pl";

sub program_info {
    print "\n\tRIFractionExpressedBatch.pl determines the fraction of expressed genes with statistically\n\tsignificant RI events at genes with a minimum input control or test TPM value\n\tfor all RI files in an input directory.\n\n\tAlso generates:\n\t\t- lists of statistically significant negative and positive IncDiff RI events in\n\t\t  expressed genes (includes gene expression values for each)\n\t\t- file with list of genes without statistcially significant change in RI.\n\n\tUsage: perl RIFractionExpressedBatch.pl [OPTIONS] -r <retained intron files directory (rMATS JCEC)> -e <expression file> -TPM <min TPMs condition 1,min TPMs condition 2> -SN <sample number condition 1,sample number condition 2> -f <FDR>\n\n\tRequired:\n\t\t-r <retained intron files directory>\n\t\t-e <expression file>\n\t\t-TPM <min TPMs condition 1,min TPMs condition 2>\n\t\t\tNote: if condition 1 or 2 to be not considered, enter \"\-\"\n\t\t\t(e.g. TPM 3,- or TPM -,3)\n\t\t-SN <sample number condition 1,sample number condition 2>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl RIFractionExpressedBatch.pl -r PATH/RI_files_directory -e PATH/expression.tsv -TPM 2,2 -SN 3,3 -f 0.05\n\n";
    exit;
}



sub options {  
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] eq "\-r") {
            $retained_intron_directory = $ARGV[$i+1];
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
    if (not defined ($FDR)) {
        print "\nFDR not defined!\n\n";
        program_info;
        exit;
    }
    elsif ($FDR !~ m/^\d/) {
        print "\nFDR is not numeric!\n\n";
        program_info;
        exit;
    } 
    elsif (not defined($retained_intron_directory)) {
        print "\nRetained intron file directory not defined!\n\n";
        program_info;
        exit;
    }
    elsif (not defined($expression_file)) {
        print "\nExpression file not defined!\n\n";
        program_info;
        exit;
    }
    elsif (not defined($TPM1)) {
        print "\nMinimum TPMs for condition 1 not defined!\n     Is -TPM value in correct format?\n\n";
        program_info;
        exit;
    }
    elsif (not defined($TPM2)) {
        print "\nMinimum TPMs for condition 2 not defined!\n     Is -TPM value in correct format?\n\n";
        program_info;
        exit;
    }
    elsif (not defined($cond1_ct)) {
        print "\nCondition 1 sample count not defined!\n     Is -SN value in correct format?\n\n";
        program_info;
        exit;
    }
    elsif (not defined($cond2_ct)) {
        print "\nCondition 2 sample count not defined!\n     Is -SN value in correct format?\n\n";
        program_info;
        exit;
    }
}

sub read_RI_file_dir {
    opendir my $input_dir, "$retained_intron_directory" or die "Can't open skipped exon directory: $!";
    foreach my $g (sort readdir $input_dir) {
        next if ($g eq '.' || $g eq '..' || $g eq '.DS_Store');
        my $path_files = $retained_intron_directory."/".$g;
        push(@RI_files, $path_files);
    }
    closedir $input_dir;
}

sub run_program {
    if (scalar @RI_files == 0) {
        print "NO INPUT FILES IN DIRECTORY!!\n\n";
        program_info;
        exit;
    }
    else {
        foreach my $RI_file(@RI_files) {
            print "\t- Processing RI file: ", basename($RI_file), "\n";
            `perl $program_path_name -r $RI_file -e $expression_file -TPM $TPM1,$TPM2 -SN $cond1_ct,$cond2_ct -f $FDR`;
        }
    }
}
options;
qc;
read_RI_file_dir;
print "\n\n";
run_program;
print "\n";