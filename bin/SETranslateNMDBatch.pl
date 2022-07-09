use warnings;
use strict; 
use File::Basename;  

my $skipped_exon_directory;
my $annotation_file;
my $genome_fasta_file;
my $FDR;
my @SE_files;
my $program_path_name = dirname(__FILE__)."\/SETranslateNMD.pl";

sub program_info {
    print "\n\tFor all SE files in an input SE files directory, SETranslateNMDBatch.pl will:\n\t\t- generate predicted RNA(DNA) isoform sequences based on annotated exon\n\t\t  structures upstream and downstream from skipping event\n\t\t- generate translated isoform sequences for all such events\n\t\t- generate separate listings of frameshifted protein isoforms\n\t\t- generate candidate neopeptides produced by frameshifts\n\t\t- generate lists of skipped exon protein sequences for in-frame\n\t\t  events for BATCH submission to NCBI conserved domain search\n\t\t  (https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi)\n\t\t- Identify events predicted to undergo Nonsense Mediated RNA\n\t\t  Decay (NMD) and generate statistics on skipping events predicted\n\t\t  to undergo and not undergo NMD\n\n\n\tUsage: perl SETranslateNMDBatch.pl [OPTIONS] -s <skipped exon files directory (rMATS)> -a <bed12 annotation file> -g <genome fasta file> -f <FDR>\n\n\tRequired:\n\t\t-s <skipped exon files directory>\n\t\t-a <bed12 annotation file>\n\t\t-g <genome fasta file>\n\t\t-f <FDR>\n\n\tAdditional:\n\t\t-h help\n\n\tExample: perl SETranslateNMDBatch.pl -s PATH/SE_input_files_directory -a PATH/bed12_annotation.bed -g PATH/genome.fa -f 0.05\n\n";
    exit;
}

sub options { 
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    } 
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] eq "\-p") {
            $program_path_name = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-s") {
            $skipped_exon_directory = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-a") {
            $annotation_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-g") {
            $genome_fasta_file = $ARGV[$i+1];
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
    elsif (not defined($skipped_exon_directory)) {
        print "\nSkipped exon file directory not defined!\n\n";
        program_info;
        exit;
    }
    elsif (not defined($genome_fasta_file)) {
        print "\nGenome fasta file not defined!\n\n";
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
            `perl $program_path_name -s $SE_file -a $annotation_file -g $genome_fasta_file -f $FDR`;
        }
    }
}
options;
qc;
read_SE_file_dir;
print "\n\n";
run_program;
print "\n";