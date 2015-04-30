#!/usr/bin/perl

use strict;
use warnings;

#-------------------------------------------------------------------------------------------------------------------------------------#
# STEP 1: OPEN USER-DEFINED FASTA FILE, SKIP FIRST LINE, REMOVE WHITESPACE
#-------------------------------------------------------------------------------------------------------------------------------------#

print "\n..........................................................\n        FASTA Nucleotide to Amino Acid Converter        \n..........................................................\n
	\nPlease enter the .txt FASTA file that you wish to convert.\n";
	


my $stdin = <STDIN>; 			# defines my $stdin variable as the user input
	chomp ($stdin);		

open (FILE, "<", $stdin) or die $!; 	#opens the user-defined input and will reference it as 'FILE'

while (<FILE>)  {
	next if $. == 1;		# $. means the current line number, so in this case
					# it says to skip if you are reading line number 1
		chomp;					
	my $sequence = $_;		# defines NTinput as what it's reading from the file, line by line
		$sequence =~ s/\s//g;	# removes whitespace from the file's input

#-------------------------------------------------------------------------------------------------------------------------------------#
# STEP 2: VERIFY THAT THE TOTAL # NUCs IS: A) DIVISIBLE BY 3, B) NUC CHARACTERS, C) NOT NUMBERS, D) NOT PUNCTUATION
#-------------------------------------------------------------------------------------------------------------------------------------#

my $length = length ($sequence);	# defines $length as the actual length of $sequence (which was read from FASTA file)
if (0 ne $length % 3) {			# this condition will give an error messages and exit the program if the length of $sequence 
					# cannot equally break into codons (i.e. it must be divisible by three)
	print "The length of your sequence does not break down equally into codons. Please fix and try again.\n";
	exit;
	}	

if ($sequence =~ /[B|D|E|F|H|I|J|K|L|M|M|O|P|Q|R|S|U|V|W|X|Y|Z]/i) { # If sequence contains a symbol or other invalid character an
								     # error message shows up and the program exits.
	print "You've use a non-nucleotide character. Please fix and try again.\n";
	exit;
	}

if ($sequence =~ /[^a-zA-z0-9]/) {	# if sequence contains punctuation, it will show an error message and exit the program
	print "You've used punctuation. Please fix and try again.\n";
	exit;
	}	

if ($sequence =~ /[0-9]/) { 		# If sequence contains a number, it will show an error message and exit the program
	print "You've used a number. Please fix and try again.\n";
	exit;
	}	

#-------------------------------------------------------------------------------------------------------------------------------------#
# STEP 3: DEFINE HASH, WHICH WILL BE USED TO TRANSLATE EACH NUC. CODON INTO ITS AMINO ACID SYMBOL
#-------------------------------------------------------------------------------------------------------------------------------------#
	
my %codontoaminoacid = (
    'TCA' => 'S', 
    'TCC' => 'S', 
    'TCG' => 'S', 
    'TCT' => 'S',
    'TTC' => 'F', 
    'TTT' => 'F',
    'TTA' => 'L', 
    'TTG' => 'L', 
    'TAC' => 'Y', 
    'TAT' => 'Y', 
    'TGC' => 'C', 
    'TGT' => 'C', 
    'TGG' => 'W', 
    'CTA' => 'L', 
    'CTC' => 'L', 
    'CTG' => 'L', 
    'CTT' => 'L',
    'CCA' => 'P', 
    'CCC' => 'P', 
    'CCG' => 'P', 
    'CCT' => 'P', 
    'CAC' => 'H',
    'CAT' => 'H', 
    'CAA' => 'Q', 
    'CAG' => 'Q', 
    'CGA' => 'R', 
    'CGC' => 'R', 
    'CGG' => 'R', 
    'CGT' => 'R',
    'ATA' => 'I',
    'ATC' => 'I', 
    'ATT' => 'I',
    'ATG' => 'M',
    'ACA' => 'T', 
    'ACC' => 'T',
    'ACG' => 'T',
    'ACT' => 'T',
    'AAC' => 'N', 
    'AAT' => 'N', 
    'AAA' => 'K',
    'AAG' => 'K',
    'AGC' => 'S',
    'AGT' => 'S',
    'AGA' => 'R', 
    'AGG' => 'R',
    'GTA' => 'V', 
    'GTC' => 'V', 
    'GTG' => 'V', 
    'GTT' => 'V', 
    'GCA' => 'A', 
    'GCC' => 'A',
    'GCG' => 'A',
    'GCT' => 'A',
    'GAC' => 'D',
    'GAT' => 'D', 
    'GAA' => 'E', 
    'GAG' => 'E',
    'GGA' => 'G',
    'GGC' => 'G', 
    'GGG' => 'G',
    'GGT' => 'G',
    'TAA' => '<STOP>', 
    'TAG' => '<STOP>', 
    'TGA' => '<STOP>', 
    );

# These are the specifications for my substition and processing.
# Reads string three values at a time, and translates them to the
# values as specified in the hash.

# Defines the translated product as $protein.

#-------------------------------------------------------------------------------------------------------------------------------------#
# STEP 4: INITIALIZE $PROTEIN VARIABLE AS BLANK, CREATE FOR LOOP TO TRANSLATE CODON
#-------------------------------------------------------------------------------------------------------------------------------------#

my ($protein) = "";			# Defines $protein as empty variable;

for (my $i=0; $i<length($sequence)-2; $i+=3)	# for loop, stating that: starting at place 0, so long as (where you are) is less than
						# the length of the imported $sequence - 2, you want to continue to the next set of 
						# three nucleotides
{
my $codon = substr($sequence,$i,3); 	# this is a substring, which allows you to define codon relative to the starting place and
					# which value relative to the first you wish to include. (in this one
					# the substring reads as following:
					# (reading $sequence, start at the $i-th letter from the left, continue to the third letter)

$protein .= $codontoaminoacid {$codon}; # this uses the hash to convert each codon, as defined in the substring above, and add it 
					# to the previously empty variable, $protein.
}


# Prints $protein to STDOUT.

sleep 1;

print "\nHow should I ouput the converted FASTA file?\n...save at the end of your original txt file. (enter '1')\n...as its own txt file.                       (enter '2')\n...just print it here on the command line.    (enter '3')\n";

my $oquestion=<STDIN>;
chomp $oquestion;




if ($oquestion eq 1){	
	print "\nAll done!\n\nYou can check out your protein sequence at the bottom of\nyour original file, $stdin.\n\nPlease press 'enter' to finish.\n";
	open (OUTPUT, ">>", $stdin) or die 
		"\nSorry, unable to open the specified file. 
		\nPlease try again. \n Reason: $! \n\n"; 

	print OUTPUT "\n\nYour translated protein is:\n$protein\n";
	<>;

	close OUTPUT or die 
		"\nSomething went wrong. \nReason: $! \n\n"; 
	

	
	}

if ($oquestion eq 2){	
	print "\nGreat! What should I name this new file?\n";
	my $newfile = <STDIN>;
	chomp $newfile;
	print "\nPlease press 'enter' to finish.\n";
	open (OUTPUT, ">", $newfile) or die 
		"\nSorry, unable to open the specified file. 
		\nPlease try again. \n Reason: $! \n\n"; 

	print OUTPUT "$protein\n";
	<>;


	
	close OUTPUT or die 
		"\nSomething went wrong. \nReason: $! \n\n"; 
	}

if ($oquestion eq 3){	
	print "\nYour nucleotide sequence was:\n $sequence\n";
	print "\nYour codon to amino acid conversion:\n $protein\n";
	}

if ($oquestion > 3){	
	print "\nThat's not a valid choice.  Please enter either 1, 2, or 3 next time.\n";
	}	
	
	

}


	
	

#if ($originalornot eq "N"||"No"||"no"){
#	print "You either entered no or put a different character.\nWhat should we save the file as?\n";
#}

#else {}

__END__

my $savefile = <STDIN>;
	chomp ($savefile);

if (length $savefile <1){	#again, if blank, open the original file and overwrite it with the converted file
	open (FASTAOUTPUT, ">", $originalsave) or die 
		"\nSorry, unable to open the specified file. 
		\nPlease try again. \n Reason: $! \n\n"; 

	print FASTAOUTPUT "$FASTAprotein\n";
	<>;

	close FASTAOUTPUT or die 
		"\nSomething went wrong. \nReason: $! \n\n"; 
}


else {	#if not blank, open the file specified, if it exists. if not, die.
	open (FASTAOUTPUT, "+<", $savefile) or die 
		"\nSorry, unable to open the specified file. 
		\nPlease try again. \n Reason: $! \n\n";

	print FASTAOUTPUT "$FASTAprotein\n" or die $!;
	<>;

	close FASTAOUTPUT or die 
		"\nSomething went wrong. \nReason: $! \n\n"; 
}


__END__








