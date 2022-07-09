use warnings;
use strict; 
use File::Basename;  

my $skipped_exon_directory;
my $annotation_file;
my $expression_file;
my $TPM1;
my $TPM2;
my $cond1_ct;
my $cond2_ct;
my $FDR;
my @SE_files;
my $program_path_name = dirname(__FILE__)."\/SEFractionExpressed.pl";

sub program_info {
    print "\n\tSEFractionExpressedBatch.pl determines the fraction of expressed genes with statistically\n\tsignificant SE events at genes with a minimum input control or test TPM value for all\n\tSE files in an input directory.\n\n\tAlso generates:\n\t\t- lists of statistically significant negative and positive IncDiff SE events in\n\t\t  expressed genes (includes gene expression values for each)\n\t\t- file with list of genes without statistcially significant exon skipping.\n\n\n\tUsage: perl SEFractionExpressedBatch.pl [OPTIONS] -s <skipped exon files directory (rMATS JCEC)> -e <expression file> -TPM <min TPMs condition 1,min TPMs condition 2> -SN <sample number condition 1,sample number condition 2> -f <FDR>\n\n\tRequired:\n\t\t-s <skipped exon files directory>\n\t\t-e <expression file>\n\t\t-TPM <min TPMs condition 1,min TPMs condition 2>\n\t\t\tNote: if condition 1 or 2 to be not considered, enter \"\-\"\n\t\t\t(e.g. TPM 3,- or TPM -,3)\n\t\t-SN <sample number condition 1,sample number condition 2>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl SEFractionExpressedBatch.pl -s PATH/SE_input_files_directory -e PATH/expression.tsv -TPM 2,2 -SN 3,3 -f 0.05\n\n";
    exit;
}

sub options { 
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    } 
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] eq "\-s") {
            $skipped_exon_directory = $ARGV[$i+1];
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
    elsif (not defined($skipped_exon_directory)) {
        print "\nSkipped exon directory not defined!\n\n";
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
}

sub read_SE_file_dir {
    opendir my $input_dir, "$skipped_exon_directory" or die "Can't open skipped exon directory: $!";
    foreach my $g (sort readdir $input_dir) {
        next if ($g eq '.' || $g eq '..' || $g eq '.DS_Store');
        my $path_files = $skipped_exon_directory."/".$g;
        push(@SE_files, $path_files);
    }
    closedir $input_dir;
}

sub run_program {
    if (scalar @SE_files == 0) {
        print "NO INPUT FILES IN DIRECTORY!!\n\n";
        program_info;
        exit;
    }
    else {
        foreach my $SE_file(@SE_files) {
            print "\t- Processing SE file: ", basename($SE_file), "\n";
            `perl $program_path_name -s $SE_file -e $expression_file -TPM $TPM1,$TPM2 -SN $cond1_ct,$cond2_ct -f $FDR`;
        }
    }
}
options;
qc;
read_SE_file_dir;
print "\n\n";
run_program;
print "\n";