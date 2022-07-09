use warnings;
use strict; 
use File::Basename;  

my $skipped_exon_directory;
my $annotation_file;
my $FDR;
my @SE_files;
my $program_path_name = dirname(__FILE__)."\/SEIntronExonSizes.pl";

sub program_info {
    print "\n\tSEIntronExonSizesBatch.pl will generate the following for all SE files in an input files directory:\n\t\t- lists of sizes for upstream exon, upstream intron, skipped exon, downstream intron,\n\t\t  and downstream exon for statistically significant positive and\n\t\t  negative IncDiff skipped exon events.\n\t\t- summary file with average sizes (for comparison, includes average sizes for potential\n\t\t  SE events based on input annotation file).\n\n\tUsage: perl SEIntronExonSizesBatch.pl [OPTIONS] -s <skipped exon files directory (rMATS JCEC)> -a <bed12 annotation file> -f <FDR>\n\n\tRequired:\n\t\t-s <skipped exon files directory>\n\t\t-a <bed12 annotation file>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl SEIntronExonSizesBatch.pl -s PATH/SE_input_files_directory -a PATH/bed12_annotation.bed -f 0.05\n\n";
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
    elsif (not defined($skipped_exon_directory)) {
        print "\nSkipped exon directory not defined!\n\n";
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
            `perl $program_path_name -s $SE_file -a $annotation_file -f $FDR`;
        }
    }
}
options;
qc;
read_SE_file_dir;
print "\n\n";
run_program;
print "\n";