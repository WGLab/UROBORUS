#!/usr/bin/perl
# This programm find circular RNA in RNA-seq;
use strict;
use warnings;
use Getopt::Long;

####################################################################################
# USAGE
#					UROBORUS.pl 1.0.0
#	circRNA identification tool in total RNA-seq data (or Ploy(A)-minus RNA-seq data)
my $USAGE =<<USAGE;
	UROBORUS.pl 0.1.3         (July.11,2016)
	usage:
	UROBORUS.pl -index /path/genome -gtf /path/genes.gtf -fasta /path/ unmapped.sam accepted_hits.bam
	-index:	genome index(use bowtie1 index);
	-gtf:	gene anotation file(*.gtf file);
	-fasta:	path for genome sequence in fasta file (*.fa) in separate chromosome;
	-p:	threads (Integer,default = 6);
	-temp:	keeping the temperary file;
	-help:	usage help;

USAGE
#contact: xfsong@nuaa.edu.cn
####################################################################################
#Version 0.1.3 modification:
#correct a bug: filtered the circRNA list, the results of circRNA cadidates were more than two reads supported.
#####################################################################################
	my ( $opt_index, $opt_gtf, $opt_fasta,$opt_p,$opt_help,$opt_temp);
	$opt_p = 6;
	$opt_help = 0;
	$opt_temp = 0;

	GetOptions(
		'index=s' => \$opt_index,
		'gtf=s' => \$opt_gtf,
		'fasta=s' => \$opt_fasta,
		'p=i' => \$opt_p,
		'help' => \$opt_help,
		'temp' => \$opt_temp,
	);
	
	if ($opt_help) { print "$USAGE\n"; exit 0; }
	
	if (!$opt_index ) { print "UROBORUS warning: missing parameters -index!","\n"; exit 0; }
	
	if (!$opt_gtf ) { print "UROBORUS warning: missing parameters -gtf!","\n"; exit 0; }
	
	if ( !$opt_fasta ) {print "UROBORUS warning: missing parameters -fasta!","\n"; exit 0;}
	
	if ( !$ARGV[0] ) { print "UROBORUS warning: missing unmapped.sam and accepted_hits.bam file!","\n";exit 0;}

        if ( !$ARGV[1] ) { print "UROBORUS warning: missing accepted_hits.bam or unmapped.sam file!","\n";exit 0;}
	
	my $then = time();
	my $when;
	$when = localtime();

	#set read_length_label = 0 if <100bp; read_length_label = 1 if >100bp
	open (RLL_IN, "$ARGV[0]") or die ("The input file should be *.sam file,please check the input file!\n");#input unmapped.sam file
	my $read_length_label = 0;
	my $main_line = <RLL_IN>;
	chomp $main_line;
	my @main_array = split("\t",$main_line);
	my $main_read_length = length($main_array[9]);
	
	if ( $main_read_length > 200 ) {die "Now UROBORUS can only support below 200bp read length total RNA-seq data!";}
	
	if ( $main_read_length > 100 )
	{
		$read_length_label = 1;
	}
	# end
	close(RLL_IN);
	
	
	
	#print STDOUT "[",$when,"]";i
	print STDERR "[",$when,"]";
	print STDERR " Preparing the reads..........";
	if ($read_length_label == 0) {&trimmer0};
	if ($read_length_label == 1) {&trimmer1};
	print STDERR "finished!","\n";

	$when = localtime();
	print STDERR "[",$when,"]";
	print STDERR " Now tophat..........";
	system("tophat --bowtie1 -p 6 -o seed_mapped_out --no-coverage-search $opt_index R20_1.fastq R20_2.fastq");
	print STDERR "finished!","\n";

	my $pwd_str1 = `pwd`;
	chomp $pwd_str1;
	my $pwd_str2 = join("",$pwd_str1,"/seed_mapped_out");
	chdir($pwd_str2); 
	system("samtools view accepted_hits.bam > accepted_hits.sam");
	chdir($pwd_str1); 

	$when = localtime();
	print STDERR "[",$when,"]";
	print STDERR " Filtering the mapped seeds...........";
	&mapped_seeds_filter;
	print STDERR "finished!","\n";

	$when = localtime();
	print STDERR "[",$when,"]";
	print STDERR " Locating the seeds in the genome...........";
	&locating("filtered_mapped_paired_seeds.sam","paired_seeds_location");
	if ($read_length_label == 0) {&locating("filtered_mapped_single_seed.sam","single_seed_location");}
	print STDERR "finished!","\n";

	$when = localtime();
	print STDERR "[",$when,"]";
	print STDERR " Extending the seeds...........";
	&paired_seeds_extension;#("paired_seeds_location","paired_seeds_extension");
	exit;

	if ($read_length_label == 0) 
		{&single_seed_extension;} #(input:single_seed_location.sam;      output:IMR_1.fastq,IMR_2.fastq);
	print STDERR "finished!","\n";

	if ($read_length_label == 0)
	{
		$when = localtime();
		print STDERR "[",$when,"]";
		print STDERR " Now bowtie..........";
	 	system("bowtie -p $opt_p -n 0 -M 20 -I 100 -X 100000 $opt_index --ff -1 IMR_1.fastq -2 IMR_2.fastq -S --sam-nohead > mapped_IMR.sam");
		print STDERR "finished!","\n";

		$when = localtime();
		print STDERR "[",$when,"]";
		print STDERR " Now filtering mapped IMR..........";
		&mapped_IMR_filter;
		print STDERR "finished!","\n";

	
		$when = localtime();
		print STDERR "[",$when,"]";
		print STDERR " Now sorting mapped IMR sam file..........";
		&sorting_samfile("filtered_mapped_IMR.sam","sorted_filtered_mapped_IMR.sam");
		print STDERR "finished!","\n";

		$when = localtime();
		print STDERR "[",$when,"]";
		print STDERR " Now locating mapped IMR..........";
		&locating("sorted_filtered_mapped_IMR.sam","mapped_IMR_location");
		print STDERR "finished!","\n";
	}
	
	$when = localtime();
	print STDERR "[",$when,"]";
	print STDERR " Now mating the reads...........";
	&paired_seeds_mating;   #input file: paired_seeds_extension.txt; output file: paired_seeds_to_reads;
	if ($read_length_label == 0) { &mapped_IMR_mating; }  #input file: mapped_IMR_location.txt; output file: IMR_to_reads.txt;
	print STDERR "finished!","\n";

	$when = localtime();
        print STDERR "[",$when,"]";
        print STDERR " Now filtering the accepted reads...........";
        &accepted_reads_processing; #("accepted_hits.bam","accepted_reads.txt");
	$when = localtime();
        print STDERR "[",$when,"]\n";


	if ($read_length_label == 0)
	{
		$when = localtime();
		print STDERR "[",$when,"]";
		print STDERR " Merging the results...........";
		unlink "mapped_reads_list.txt" if -e "mapped_reads_list.txt";
		system("cat paired_seeds_to_reads.txt IMR_to_reads.txt accepted_reads.txt >> mapped_reads_list.txt");
		print STDERR "finished!","\n";
	} elsif ($read_length_label == 1)
        {
                $when = localtime();
                print STDERR "[",$when,"]";
                print STDERR " Merging the results...........";
                unlink "mapped_reads_list.txt" if -e "mapped_reads_list.txt";
                system("cat paired_seeds_to_reads.txt accepted_reads.txt >> mapped_reads_list.txt");
                print STDERR "finished!","\n";
        }

	$when = localtime();
	print STDERR "[",$when,"]";
	print STDERR " Making the statistics...........";
	&cR_statistics(); 
	print STDERR "finished!","\n";

	$when = localtime();
	print STDERR "[",$when,"]";
	print STDERR " Output the circular RNA list!","\n";

	#deleting temp file!
	if (!$opt_temp)
	{
		unlink "filtered_mapped_paired_seeds.sam" if -e "filtered_mapped_paired_seeds.sam";
		unlink "filtered_mapped_single_seed.sam" if -e "filtered_mapped_single_seed.sam";
		unlink "paired_seeds_location.sam" if -e "paired_seeds_location.sam";
		unlink "single_seed_location.sam" if -e "single_seed_location.sam";
		unlink "R20_1.fastq" if -e "R20_1.fastq";
		unlink "R20_2.fastq" if -e "R20_2.fastq";
		unlink "BMR_1.fastq" if -e "BMR_1.fastq";
		unlink "BMR_2.fastq" if -e "BMR_2.fastq";
		unlink "paired_seeds_extension.sam" if -e "paired_seeds_extension.sam";
		unlink "IMR_1.fastq" if -e "IMR_1.fastq";
		unlink "IMR_2.fastq" if -e "IMR_2.fastq";
		unlink "filtered_mapped_IMR.sam" if -e "filtered_mapped_IMR.sam";
		unlink "mapped_IMR.sam" if -e "mapped_IMR.sam";
		unlink "mapped_IMR_location.sam" if -e "mapped_IMR_location.sam";
		unlink "paired_seeds_to_reads.txt" if -e "paired_seeds_to_reads.txt";
		unlink "accepted_reads.txt" if -e "accepted_reads.txt";
		unlink "sorted_filtered_mapped_IMR.sam" if -e "sorted_filtered_mapped_IMR.sam";
		unlink "circRNA_list_temp.txt" if -e "circRNA_list_temp.txt";
		unlink "IMR_to_reads.txt" if -e "IMR_to_reads.txt";
		unlink "mapped_reads_list.txt" if -e "mapped_reads_list.txt";
		unlink "circRNA_list_temp2.txt" if -e "circRNA_list_temp2.txt";
	}
	
	$when = localtime();
	print STDERR "[",$when,"]";
	print STDERR " Completed successfully!  The time elapsed: about ",0.01*int ( (time()-$then)/36 )," hours.\n";

	
	

sub trimmer0
{
	open (TR_IN, "$ARGV[0]") or die ("The input file should be *.sam file,please check the input file!\n");#input unmapped.sam file
	open (TR_OUT1, ">R20_1.fastq");
	open (TR_OUT2, ">R20_2.fastq");

	my $line;
	my @array;

	my $read_length;

	while ($line = <TR_IN>) {
	# generate head 20 of reads; 
		chomp $line;
		@array = split("\t",$line);
	
		$read_length = length($array[9]);
		if ( $read_length < 40 ) {next;}
		if ( $array[0] =~ /\._\./){die "The read label contains unexpected string '._.' please remove or replace it!";}
	# generate first 20bp of reads;	
		$array[10] =~ s/\//~/g;
		my $seq_label_A = join ("._.",join ("","@",$array[0]),"A",$array[9],$array[10] );
		my $seq20_A = substr ($array[9],0,20);
		my $qual_label_A = "+";
		my $qual20_A = substr ($array[10],0,20);
		print TR_OUT1 $seq_label_A, "\n", $seq20_A, "\n", $qual_label_A, "\n", $qual20_A, "\n";
	
	# generate tail 20bp of reads;
		my $seq_label_B = join ("._.",join ("","@",$array[0]),"B",$array[9],$array[10] );
		my $seq20_B = reverse substr ( reverse ($array[9]),0,20);
		my $qual_label_B = "+";
		my $qual20_B = reverse substr ( reverse ($array[10]),0,20);
		print TR_OUT2 $seq_label_B, "\n", $seq20_B, "\n", $qual_label_B, "\n", $qual20_B, "\n";
	}
	close(TR_IN);close(TR_OUT1);close(TR_OUT2);
}
	
	
sub trimmer1
{
	open (TR_IN, "$ARGV[0]") or die ("The input file should be *.sam file,please check the input file!\n");#input unmapped.sam file
	open (TR_OUT1, ">R20_1.fastq");
	open (TR_OUT2, ">R20_2.fastq");

	my $line;
	my @array;

	my $read_length;

	while ($line = <TR_IN>) {
	# generate head 20 of reads; 
		chomp $line;
		@array = split("\t",$line);
	
		$read_length = length($array[9]);
		if ( $read_length < 40 ) {next;}
		if ( $array[0] =~ /\._\./){die "The read label contains unexpected string '._.' please remove or replace it!";}
	# generate first 20bp of reads;	
		$array[10] =~ s/\//~/g;
		my $seq_label_A = join ("._.",join ("","@",$array[0]),"A",$array[9] );
		my $seq20_A = substr ($array[9],0,20);
		my $qual_label_A = "+";
		my $qual20_A = substr ($array[10],0,20);
		print TR_OUT1 $seq_label_A, "\n", $seq20_A, "\n", $qual_label_A, "\n", $qual20_A, "\n";
	
	# generate tail 20bp of reads;
		my $seq_label_B = join ("._.",join ("","@",$array[0]),"B",$array[9] );
		my $seq20_B = reverse substr ( reverse ($array[9]),0,20);
		my $qual_label_B = "+";
		my $qual20_B = reverse substr ( reverse ($array[10]),0,20);
		print TR_OUT2 $seq_label_B, "\n", $seq20_B, "\n", $qual_label_B, "\n", $qual20_B, "\n";
	}
	close(TR_IN);close(TR_OUT1);close(TR_OUT2);
}

sub mapped_seeds_filter
{
	my $input_file = `pwd`;

	chomp $input_file;
	$input_file = join("",$input_file,"/seed_mapped_out/accepted_hits.sam");
#	$input_file = "seeds_accepted_hits.sam";

	open(IN,$input_file) or die ("Running error! Please run the tophat using paired seeds first!\n");#input accepted_hits.sam of Tophat results
	open(PAIRED_OUT,">filtered_mapped_paired_seeds.sam");
	open(SINGLE_OUT,">filtered_mapped_single_seed.sam");#for OneAnchor finding

	my $paired_count = 0;
	my $single_count = 0;
	my $line;
	my @array;
	my $bit;
	my @flag;
	my @head;
	my $flag_length;

	while ( $line = <IN> )
	{
		@array = split(/\t/, $line);
		@head = split ("", $array[0]);
	
		if ( $head[0] ne "@" )
		{
			if ( $array[2] ne "chrM")
			{
				$bit = reverse ( sprintf("%b",$array[1])+0 ); #decimal to binary; 
				@flag = split ("", $bit);
				$flag_length = @flag;
				# keeping the flag in 9 bit;
				if ( $flag_length == 8 ) { $flag[8] = 0;}
				if ( $flag_length == 7 ) { $flag[8] = 0; $flag[7] = 0;}
				if ( $flag_length == 6 ) { $flag[8] = 0; $flag[7] = 0; $flag[6] = 0;}
				if ( $flag_length == 5 ) { $flag[8] = 0; $flag[7] = 0; $flag[6] = 0; $flag[5] = 0;}  

				if  ( 
					($array[6] eq "=")   # the paired reads are all mapped in the same chromsome;
					and ($flag[8] == 0)   # the alignment is primary;
					and (abs($array[3] - $array[7]) > 20) 
					and (abs($array[3] - $array[7]) < 3000000) # within the one gene; maximam length of gene;
					)
				{
					if ( $flag[4] == "0" and $flag[5] == "0" ) # the strand is forward;
					{
						if ( $flag[6] == "1" and $flag[7] == "0" ) # read is the first read;
						{
							if ( $array[3] > $array[7] ) #the read is in the downstream of the mate
							{
								print PAIRED_OUT $line;
								$paired_count +=1;
							}
						}
						if ( $flag[6] == "0" and $flag[7] == "1" ) # read is the second read;
						{
							if ( $array[3] < $array[7] ) #the read is in the upstream of the mate
							{
								print PAIRED_OUT $line;
								$paired_count +=1;
							}
						}
					}
			
					# the strand is reverse;
					if ( $flag[4] == "1" and $flag[5] == "1" ) # the strand is reverse;
					{
						if ( $flag[6] == "1" and $flag[7] == "0" ) # read is the first read;
						{
							if ( $array[3] < $array[7] ) #the read is in the uptream of the mate
							{
								print PAIRED_OUT $line;
								$paired_count +=1;
							}
						}
						if ( $flag[6] == "0" and $flag[7] == "1" ) # read is the second read;
						{
							if ( $array[3] > $array[7] ) #the read is in the downstream of the mate
							{
								print PAIRED_OUT $line;
								$paired_count +=1;
							}
						}
					}
				}
				
				#The following is to find single mapped seed in the unmapped.sam.
				if ( ($flag[8] == 0) ) # deleted ($array[6] eq "*") and in cR7_noFQ version
				{
					print SINGLE_OUT $line;
					$single_count +=1;
				}
				#The above is to find single mapped seed in the unmapped.sam.		
			}
		}
	}
	close(IN);
	close(PAIRED_OUT);
	close(SINGLE_OUT);
}

sub locating
{
	my $gtf_file = $opt_gtf;
	open (GTF,"$gtf_file");
	my $out_file = join("",$_[1],".sam");

	open (IN, "$_[0]");
	open (OUT,">$out_file");

	my $count;
	my $sam_line;
	my $gtf_line;
	my @sam_array;
	my @gtf_array;
	my $bit;
	my @flag;
	my @head;
	my $flag_length;
	my $strand;
	my $last_sam_chr = "chr1";
	my $index;
	$count = 0;

	$gtf_line = <GTF>;
	@gtf_array = split(/\s+/, $gtf_line);

	while ( $sam_line = <IN> )
	{
		chomp $sam_line;
		@sam_array = split(/\t/, $sam_line);
		if ( $sam_array[2] eq "chrM" ) 	{ next; }
		
		LABEL:
		
		if ( $sam_array[2] eq $gtf_array[0] )
		{
			if ( $gtf_array[2] eq "exon" ) 
			{
				if ( ( $sam_array[3] >= $gtf_array[3] - 4 ) and ($sam_array[3] < $gtf_array[4] + 4)  )
				{
					$count += 1;
					#the strand;below is not important;
					$bit = reverse ( sprintf("%b",$sam_array[1])+0 ); #decimal to binary; 
					@flag = split ("", $bit);
					$flag_length = @flag;
					# keeping the flag in 9 bit;
					if ( $flag_length == 8 ) { $flag[8] = 0;}
					if ( $flag_length == 7 ) { $flag[8] = 0; $flag[7] = 0;}
					if ( $flag_length == 6 ) { $flag[8] = 0; $flag[7] = 0; $flag[6] = 0;}
					if ( $flag_length == 5 ) { $flag[8] = 0; $flag[7] = 0; $flag[6] = 0; $flag[5] = 0;}  
					#If the strand is forward or reverse;
					if ( $sam_array[6] eq "=")# for 2 mapped seeds
					{
						if ( $flag[4] eq "0" and $flag[5] eq "0" ) { $strand = "+"; }# the strand is forward;
						if ( $flag[4] eq "1" and $flag[5] eq "1" ) { $strand = "-"; }# the strand is reverse;
					}
					if ( $sam_array[6] ne "=")# for 1 mapped seed
					{
						if ( $flag[4] eq "0" ){ $strand = "+"; }# the strand is forward;
						if ( $flag[4] eq "1" ){ $strand = "-"; }# the strand is reverse;
					}
					#the strand;above is not important;
					$gtf_array[11] =~ s/^\"|\";$//g;
					if ( $gtf_array[12] eq "transcript_id")
					{
						$index = 13;
					}
					else
					{
						$index = 15;
					}
					$gtf_array[$index] =~ s/^\"|\";$//g;
					print OUT $sam_line,"\t","XC:Z:",join(" ",$gtf_array[0],$gtf_array[3], "*", $gtf_array[4], "*", join("","strand",$gtf_array[6]), $gtf_array[11], join("","read",$strand), $gtf_array[$index] ),"\t","\n";
				}
				else
				{
					if ( $sam_array[3] > $gtf_array[4] )
					{
						$gtf_line = <GTF>;
			
						if ( defined $gtf_line ) 
						{
							@gtf_array = split(/\s+/, $gtf_line);
							goto LABEL;
						}
						else { last; }
					}
				}
			}
			else
			{
				$gtf_line = <GTF>;
				if ( defined $gtf_line )
				{
					@gtf_array = split(/\s+/, $gtf_line);
					goto LABEL;
				}
				else { last; }
			}
		}
		else
		{
			if ( $sam_array[2] ne $last_sam_chr )
			{
				$gtf_line = <GTF>;
				if ( defined $gtf_line ) 
				{
					@gtf_array = split(/\s+/, $gtf_line);
					goto LABEL;
				}
				else {last;}
			}
		}
		$last_sam_chr = $sam_array[2];
	}
	close(IN);
	close(GTF);
	close(OUT);
}

sub paired_seeds_extension
{
	open (SAM_IN,"paired_seeds_location.sam");
	open (SAM_OUT,">paired_seeds_extension.sam");
#	open (FQ_OUT1,">BMR_1.fastq");
#	open (FQ_OUT2,">BMR_2.fastq");

	my $mismatch_num;
	my $Lower_exon;
	my $Uper_exon;
	my $extension_length;
	my $element_num;
	my $new_sam_line;
	my $exon_extension;
	my $read_extension;
#	my $qual_extension;

	my $read_length;
	my $line;
	my $sam_line;
	my @array;
	my @sam_array;
	my @read_label;
	my $read_strand;
	my $i;

	while ( $sam_line = <SAM_IN> )
	{
		chomp $sam_line;
		if ( $sam_line =~ /XC:Z:([^\t]+)/ )
		{
			$line = $1;
			@array = split (/\s/,$line);
		}
		
		@sam_array = split (/\t/,$sam_line);
		@read_label = split (/\._\./,$sam_array[0]);
		$read_strand = $array[7];
		$read_strand =~ s/read//;
		
		$read_length = length($read_label[2]);
		
		if ( $read_strand eq "+" )#forward(+)
		{	
			if ( ($read_label[1] eq "B") and (($sam_array[3] - $array[1]) < ($read_length - 20)) )
	#			or ( ($sam_array[6] eq "*") and ($read_label[1] eq "B") and (($sam_array[3] - $array[4]) < ($read_length - 20)) and (($sam_array[3] - $array[4]) > 20)) ) 
			{
				$Lower_exon = &fetch_Lexon_seq($sam_array[2],$array[1]);
				$extension_length = $sam_array[3]-$array[1];
				$exon_extension = substr ( $Lower_exon,0,$extension_length);
				$read_extension = reverse substr ( reverse($read_label[2]),20,$extension_length );
#				$qual_extension = reverse substr ( reverse($read_label[3]),20,$extension_length );
				$mismatch_num=&seq_comp($exon_extension,$read_extension);
				
				#for output BMR_1.fastq file
#				$sam_array[9] = join ("",$read_extension,$sam_array[9]);
#				$sam_array[10] = join ("",$qual_extension,$sam_array[10]);
				
#				print FQ_OUT1 join("","@",$read_label[0]),"\n";
#				print FQ_OUT1 $sam_array[9],"\n";
#				print FQ_OUT1 "+","\n";
#				print FQ_OUT1 $sam_array[10],"\n";
				
				#for output BMR_1.fastq file
				
				if ( $mismatch_num < 3 )
				{
					$array[2] = join("",$extension_length+20,"m");
					$array[0] = join("","XC:Z:",$array[0]);
					$line = join (" ",$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],$array[7],$array[8]);
					$sam_line =~ s/XC:Z:[^\t]+/$line/;
					print SAM_OUT $sam_line,"\t","\n";
				}
			}
			elsif ( ($read_label[1] eq "A") and (($array[3]-$sam_array[3]) < ($read_length - 1)) )
			{
				$Uper_exon = &fetch_Rexon_seq($sam_array[2],$array[3]);
				$extension_length = $array[3] -$sam_array[3]-20+1;
				if ($extension_length < 0) { $extension_length = 0; }
				$exon_extension = reverse substr ( reverse($Uper_exon),0,$extension_length );
				$read_extension = substr($read_label[2],20,$extension_length);
#				$qual_extension = substr($read_label[3],20,$extension_length);
				$mismatch_num=&seq_comp($exon_extension,$read_extension);
				
				#for output BMR_2.fastq file
#				$sam_array[9] = join ("",$sam_array[9],$read_extension);
#				$sam_array[10] = join ("",$sam_array[10],$qual_extension);
				
#				print FQ_OUT2 join("","@",$read_label[0]),"\n";
#				print FQ_OUT2 $sam_array[9],"\n";
#				print FQ_OUT2 "+","\n";
#				print FQ_OUT2 $sam_array[10],"\n";
				#for output BMR_2.fastq file
				
				if ( $mismatch_num < 3 )
				{
					$array[4] = join("",$extension_length+20,"m");
					$array[0] = join("","XC:Z:",$array[0]);
					$line = join (" ",$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],$array[7],$array[8]);
					$sam_line =~ s/XC:Z:[^\t]+/$line/;
					print SAM_OUT $sam_line,"\t","\n";
				}
			}
		}
		else # reverse (-)
		{
			if ( ( $read_label[1] eq "A" ) and (($sam_array[3] - $array[1]) < ($read_length - 20)) )
			{
				$Lower_exon = &fetch_Lexon_seq($sam_array[2],$array[1]);
				$extension_length = $sam_array[3]-$array[1];
				$exon_extension = substr($Lower_exon,0,$extension_length);
				$read_label[2] =~ tr/ATCGatcg/TAGCtagc/;
				
				$read_extension = reverse substr ( $read_label[2],20,$extension_length );
#				$qual_extension = reverse substr ( $read_label[3],20,$extension_length );
				$mismatch_num=&seq_comp($exon_extension,$read_extension);
				
				#for output BMR_1.fastq file
#				$sam_array[9] = join ("",$read_extension,$sam_array[9]);
#				$sam_array[10] = join ("",$qual_extension,$sam_array[10]);
#				print FQ_OUT1 join("","@",$read_label[0]),"\n";
#				print FQ_OUT1 $sam_array[9],"\n";
#				print FQ_OUT1 "+","\n";
#				print FQ_OUT1 $sam_array[10],"\n";
				#for output BMR_1.fastq file
	
				if ( $mismatch_num < 3 )
				{
					$array[2] = join("",$extension_length+20,"m");
					$array[0] = join("","XC:Z:",$array[0]);
					$line = join (" ",$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],$array[7],$array[8]);
					$sam_line =~ s/XC:Z:[^\t]+/$line/;
					print SAM_OUT $sam_line,"\t","\n";
				}
			}
			elsif ( ($read_label[1] eq "B") and (($array[3]-$sam_array[3]) < ($read_length - 1)) )
			{
				$Uper_exon = &fetch_Rexon_seq($sam_array[2],$array[3]);
				$extension_length = $array[3] -$sam_array[3]-20+1;
				if ($extension_length < 0) { $extension_length = 0; }
				$exon_extension = reverse substr ( reverse($Uper_exon),0,$extension_length );
				$read_label[2] =~ tr/ATCGatcg/TAGCtagc/;
				
				$read_extension = substr ( reverse ($read_label[2]),20,$extension_length );
#				$qual_extension = substr ( reverse ($read_label[3]),20,$extension_length );
				$mismatch_num=&seq_comp($exon_extension,$read_extension);
				
				#for output BMR_2.fastq file
#				$sam_array[9] = join ("",$sam_array[9],$read_extension);
#				$sam_array[10] = join ("",$sam_array[10],$qual_extension);
#				print FQ_OUT2 join("","@",$read_label[0]),"\n";
#				print FQ_OUT2 $sam_array[9],"\n";
#				print FQ_OUT2 "+","\n";
#				print FQ_OUT2 $sam_array[10],"\n";
				#for output BMR_2.fastq file
				
				if ( $mismatch_num < 3 )
				{				
					$array[4] = join("",$extension_length+20,"m");
					$array[0] = join("","XC:Z:",$array[0]);
					$line = join (" ",$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],$array[7],$array[8]);
					$sam_line =~ s/XC:Z:[^\t]+/$line/;
					print SAM_OUT $sam_line,"\t","\n";
				}
			}			
		}
	}
	close(SAM_IN);
	close(SAM_OUT);
#	close(FQ_OUT1);
#	close(FQ_OUT2);
}

sub paired_seeds_mating
{
	open (IN, "paired_seeds_extension.sam");
	open (OUT,">paired_seeds_to_reads.txt");
	my @file= <IN>;
	my $line_num = @file;
	my $sam_line_i;
	my $sam_line_j;
	my $line_i;
	my $line_j;
	my @sam_array_i;
	my @sam_array_j;
	my @array_i;
	my @array_j;
	my $i;
	my $j;

	my $read_site;
	my @read_label;
	my $read_gene;
	my $read_length;
	my $mate_site;
	my @mate_label;
	my $mate_gene;

	for ( $i = 0; $i < $line_num; $i +=1 )
	{
		chomp $file[$i];
		$sam_line_i = $file[$i];
		
		if ( $sam_line_i =~ /XC:Z:([^\t]+)/ )
		{
			$line_i = $1;
			@array_i = split (/\s/,$line_i);
		}
		
		@sam_array_i = split(/\t/,$file[$i]);
		
		$mate_site = $sam_array_i[7]; 

		@read_label = split (/\._\./,$sam_array_i[0] );
		$read_gene = $array_i[6];
		$read_length = length($read_label[2]);
		
		if ( $sam_array_i[3] > $mate_site ) {next;} # if mate site is in the downstream of read site, then skip;
		
		for ( $j = $i+1; $j < $line_num; $j +=1 )
		{
			chomp $file[$j];
			$sam_line_j = $file[$j];
		
			if ( $sam_line_j =~ /XC:Z:([^\t]+)/ )
			{
				$line_j = $1;
				@array_j = split (/\s/,$line_j);
			}
			
			@sam_array_j = split(/\t/,$file[$j]);
			
			$read_site = $sam_array_j[7];
			@mate_label = split (/\._\./,$sam_array_j[0] );
			$mate_gene = $array_j[6];
			
			if ( $sam_array_j[2] ne $sam_array_i[2] ) {last;} # if the read and mate are not in the same chromosome, then ending mate seaching loop;
			
			if ( $sam_array_j[3] < $read_site ) {next;} # if read site is in the uptream of read site, then skip;
			
			#if read and mate are real pair-end seq, then output circRNA junction info;
			if ( ( $sam_array_j[3] eq $mate_site ) and ( $sam_array_i[3] eq $read_site) ) 
			{
				if ( ($mate_gene eq $read_gene) and ($mate_label[0] eq $read_label[0]) )
				{
					#output circRNA_Junctions.txt
					#/Chr0/Start1/End2/Strand3/GeneName4/GenomicLength5/read number6/matched Transcript_id7/
					$array_i[5] =~ s/strand//; 
					$array_i[1] = $array_i[1]-1; # 0-based 
					print OUT join ("\t",$array_i[0],$array_i[1],$array_j[3],$array_i[5],$array_i[6],$array_j[3]-$array_i[1],"1",$array_i[8],);
					print OUT "\n";
					last;
				}
			}
		}
	}
	close(IN);
	close(OUT);
}

sub single_seed_extension
{
	open (SAM_IN,"single_seed_location.sam");
	open (FQ_OUT1,">IMR_1.fastq");
	open (FQ_OUT2,">IMR_2.fastq");
	my $extension_length;
	my $read_extension;
	my $qual_extension;
	my $read_length;

	my $line;
	my $sam_line;
	my @array;
	my @sam_array;
	my @read_label;
	my $read_strand;
	my $mate_length;
	my $mate_label;
	my $mate_seq;
	my $mate_qual;

	while ( $sam_line = <SAM_IN> )
	{
		chomp $sam_line;
		
		if ( $sam_line =~ /XC:Z:([^\t]+)/ )
		{
			$line = $1;
			@array = split (/\s/,$line);
		}
		
		@sam_array = split (/\t/,$sam_line);
		@read_label = split (/\._\./,$sam_array[0]);
		$read_strand = $array[7];
		$read_strand =~ s/read//;
		
		$read_length = length($read_label[2]);
		$read_label[3] =~ s/\~/\//g; #restrore the quality score bit trimmed by tophat ;
		
		if ( length($read_label[3]) < $read_length  ) { next;}
		
		if ( $read_strand eq "+" )#forward(+)
		{	
			if	( ($read_label[1] eq "B") and (($sam_array[3] - $array[1]) < ($read_length - 19)) and (($sam_array[3] - $array[1]) > ($read_length - 40) ) )
			{
				$extension_length = $sam_array[3]-$array[1];
				$read_extension = reverse substr ( reverse($read_label[2]),20,$extension_length );
				$qual_extension = reverse substr ( reverse($read_label[3]),20,$extension_length );
				
				$sam_array[9] = join ("",$read_extension,$sam_array[9]);
				$sam_array[10] = join ("",$qual_extension,$sam_array[10]);
				
				$mate_length = $read_length - length($sam_array[9]);
				if ( $mate_length > 3 )
				{
					$mate_label = join("","@",$read_label[0]);
					$mate_seq = substr($read_label[2],0,$mate_length);
					
					$mate_qual =substr($read_label[3],0,$mate_length);
					print FQ_OUT2 $mate_label,"\n";
					print FQ_OUT2 $mate_seq,"\n";
					print FQ_OUT2 "+","\n";
					print FQ_OUT2 $mate_qual,"\n";
						
					print FQ_OUT1 join("","@",$read_label[0]),"\n";
					print FQ_OUT1 $sam_array[9],"\n";
					print FQ_OUT1 "+","\n";
					print FQ_OUT1 $sam_array[10],"\n";
				}
			}
			elsif ( ($read_label[1] eq "A") and (($array[3]-$sam_array[3]) < ($read_length - 1)) and (($array[1]-$sam_array[3]) > ($read_length - 21)) )
			{

				$extension_length = $array[3] -$sam_array[3]-20+1;
				if ($extension_length < 0) { $extension_length = 0; }

				$read_extension = substr($read_label[2],20,$extension_length);
				$qual_extension = substr($read_label[3],20,$extension_length);
				
				$sam_array[9] = join ("",$sam_array[9],$read_extension);
				$sam_array[10] = join ("",$sam_array[10],$qual_extension);
			
				$mate_length = $read_length - length($sam_array[9]);
				if ( $mate_length > 3 )	
				{	
					$mate_label = join("","@",$read_label[0]);
					$mate_seq = substr($read_label[2],0,$mate_length);
					$mate_qual =substr($read_label[3],0,$mate_length);
					print FQ_OUT1 $mate_label,"\n";
					print FQ_OUT1 $mate_seq,"\n";
					print FQ_OUT1 "+","\n";
					print FQ_OUT1 $mate_qual,"\n";
						
					print FQ_OUT2 join("","@",$read_label[0]),"\n";
					print FQ_OUT2 $sam_array[9],"\n";
					print FQ_OUT2 "+","\n";
					print FQ_OUT2 $sam_array[10],"\n";
				}
			}
		}
		else # reverse (-)
		{
			if ( ($read_label[1] eq "A") and (($sam_array[3] - $array[1]) < ($read_length - 19)) and (($sam_array[3] - $array[1]) > ($read_length - 40)) )
			{

				$extension_length = $sam_array[3]-$array[1];
				$read_label[2] =~ tr/ATCGatcg/TAGCtagc/;
				
				$read_extension = reverse substr ( $read_label[2],20,$extension_length );
				$qual_extension = reverse substr ( $read_label[2],20,$extension_length );
				
				$sam_array[9] = join ("",$read_extension,$sam_array[9]);
				$sam_array[10] = join ("",$qual_extension,$sam_array[10]);
		
				$mate_length = $read_length - length($sam_array[9]);
				if ( $mate_length > 3 )	
				{	
					$mate_label = join("","@",$read_label[0]);
					$mate_seq = substr(reverse($read_label[2]),0,$mate_length);
					$mate_qual =substr(reverse($read_label[3]),0,$mate_length);
					print FQ_OUT2 $mate_label,"\n";
					print FQ_OUT2 $mate_seq,"\n";
					print FQ_OUT2 "+","\n";
					print FQ_OUT2 $mate_qual,"\n";
					
					print FQ_OUT1 join("","@",$read_label[0]),"\n";
					print FQ_OUT1 $sam_array[9],"\n";
					print FQ_OUT1 "+","\n";
					print FQ_OUT1 $sam_array[10],"\n";
				}
			}
			elsif ( ($read_label[1] eq "B") and (($array[3]-$sam_array[3]) < ($read_length - 1)) and (($array[3]-$sam_array[3]) > ($read_length - 21)) )
			{

				$extension_length = $array[3] -$sam_array[3]-20+1;
				if ($extension_length < 0) { $extension_length = 0; }
				$read_label[2] =~ tr/ATCGatcg/TAGCtagc/;
				
				$read_extension = substr ( reverse ($read_label[2]),20,$extension_length );
				$qual_extension = substr ( reverse ($read_label[3]),20,$extension_length );
				
				$sam_array[9] = join ("",$sam_array[9],$read_extension);
				$sam_array[10] = join ("",$sam_array[10],$qual_extension);
				
				$mate_length = $read_length - length($sam_array[9]);
				if ( $mate_length > 3 )	
				{	
					$mate_label = join("","@",$read_label[0]);
					$mate_seq = reverse substr($read_label[2],0,$mate_length);
					$mate_qual =reverse substr($read_label[3],0,$mate_length);
					print FQ_OUT1 $mate_label,"\n";
					print FQ_OUT1 $mate_seq,"\n";
					print FQ_OUT1 "+","\n";
					print FQ_OUT1 $mate_qual,"\n";
						
					print FQ_OUT2 join("","@",$read_label[0]),"\n";
					print FQ_OUT2 $sam_array[9],"\n";
					print FQ_OUT2 "+","\n";
					print FQ_OUT2 $sam_array[10],"\n";
				}
			}			
		}
	}

	close(SAM_IN);
	close(FQ_OUT1);
	close(FQ_OUT2);
}

sub mapped_IMR_filter
{
	open(IN,"mapped_IMR.sam") or die ("Running error! Please run the Bowtie using IMR_1.fastq and IMR_2.fastq first!\n");#input mapped_IMR.sam of Bowtie results
	open(OUT,">filtered_mapped_IMR.sam");

	my $line;
	my @array;
	my $bit;
	my @flag;
	my @head;
	my $flag_length;

	while ( $line = <IN> )
	{
		@array = split(/\t/, $line);
		@head = split ("", $array[0]);
		
		if ( $array[2] ne "chrM" || $array[2] ne "*")
		{
			$bit = reverse ( sprintf("%b",$array[1])+0 ); #decimal to binary; 
			
			@flag = split ("", $bit);
		
			$flag_length = @flag;
			# keeping the flag in 9 bit;
			if ( $flag_length == 8 ) { $flag[8] = 0;}
			if ( $flag_length == 7 ) { $flag[8] = 0; $flag[7] = 0;}
			if ( $flag_length == 6 ) { $flag[8] = 0; $flag[7] = 0; $flag[6] = 0;}
			if ( $flag_length == 5 ) { $flag[8] = 0; $flag[7] = 0; $flag[6] = 0; $flag[5] = 0;}  

			if  ( 
				($array[6] eq "=")   # the paired reads are all mapped in the same chromsome;
				and ($flag[8] == 0)   # the alignment is primary;
				and (abs($array[3] - $array[7]) > 20) 
				and (abs($array[3] - $array[7]) < 100000) # within the one gene; maximam length of gene;
				)
			{
				print OUT $line;
			}
		}
	}
	close(IN);
	close(OUT);
}

sub mapped_IMR_mating
{
	open (IN, "mapped_IMR_location.sam");
	open (OUT,">IMR_to_reads.txt");

	my @file= <IN>;
	my $line_num = @file;
	my $sam_line_i;
	my $sam_line_j;
	my $line_i;
	my $line_j;
	my @sam_array_i;
	my @sam_array_j;
	my @array_i;
	my @array_j;
	my $i;
	my $j;

	my $read_site;
	my $read_label;
	my $read_gene;

	my $mate_site;
	my $mate_label;
	my $mate_gene;

	for ( $i = 0; $i < $line_num; $i +=1 )
	{
		chomp $file[$i];
		$sam_line_i = $file[$i];
		
		if ( $sam_line_i =~ /XC:Z:([^\t]+)/ )
		{
			$line_i = $1;
			@array_i = split (/\s/,$line_i);
		}
		
		@sam_array_i = split(/\t/,$file[$i]);
		
		$mate_site = $sam_array_i[7]; 

		$read_label = $sam_array_i[0];
		$read_gene = $array_i[6];
		
		if ( $sam_array_i[3] > $mate_site ) {next;} # if mate site is in the downstream of read site, then skip;
		
		for ( $j = $i+1; $j < $line_num; $j +=1 )
		{
			chomp $file[$j];
			$sam_line_j = $file[$j];
		
			if ( $sam_line_j =~ /XC:Z:([^\t]+)/ )
			{
				$line_j = $1;
				@array_j = split (/\s/,$line_j);
			}
			
			@sam_array_j = split(/\t/,$file[$j]);
			
			$read_site = $sam_array_j[7];
			$mate_label = $sam_array_j[0];
			$mate_gene = $array_j[6];
			
			if ( $sam_array_j[2] ne $sam_array_i[2] ) {last;} # if the read and mate are not in the same chromosome, then ending mate seaching loop;
			
			if ( $sam_array_j[3] < $read_site ) {next;} # if mate site is in the uptream of read site, then skip;
			
			#if read and mate are real pair-end seq, then output circRNA junction info;
			if ( ( $sam_array_j[3] eq $mate_site ) and ( $sam_array_i[3] eq $read_site) ) 
			{
				if ( ($mate_gene eq $read_gene) and ($mate_label eq $read_label) )
				{
					#output circRNA_Junctions.txt
					#/Chr0/Start1/End2/Strand3/GeneName4/GenomicLength5/read number6/matched Transcript_id7/
					$array_i[5] =~ s/strand//; 
					$array_i[1] = $array_i[1]-1;
					print OUT join ("\t",$array_i[0],$array_i[1],$array_j[3],$array_i[5],$array_i[6],$array_j[3]-$array_i[1],"1",$array_i[8],);
					print OUT "\n";
					last;
				}
			}
		}
	}
	close(IN);
	close(OUT);
}




sub sorting_samfile
{
	my $count = 0;
	my $line;
	my @sam_array;
	my $out_line;
	my @array;
	my $i;
	my $j;
	my @arr = ();
	
	open(IN,"$_[0]");
	open(OUT,">$_[1]");
	
	$line = <IN>;
	@sam_array = split(/\t/,$line);
	my $element_length = @sam_array;

	push @array,[split(/\t/,$line)];
	$count += 1;
	while ($line = <IN>)
	{
		push @array,[split(/\t/,$line)];
		$count += 1;
	}
	@arr = sort{$a->[2] cmp $b->[2] or $a->[3] <=> $b->[3]} @array;
	for($i=0; $i<$count;$i++)
	{
		$out_line = $arr[$i][0];

		for ( $j=1; $j< $element_length; $j++)
		{
			$out_line = join("\t",$out_line,$arr[$i][$j] );
		}
		print OUT $out_line;
	}
	close(IN);
	close(OUT);
}

sub accepted_reads_processing
{
        `samtools sort -n -@ $opt_p -o accepted_hits.sorted.bam $ARGV[1]`;
        `samtools view accepted_hits.sorted.bam > accepted_hits.sorted.sam`;

        my @time = localtime();
        my $month = $time[4] + 1;
        my $tmpdir = "tmp_$month\_$time[3]\_$time[2]\_$time[1]";
        `mkdir $tmpdir`;

        my @chrid = `cat $opt_gtf | awk '{print \$1;}' | sort | uniq`;
        chomp @chrid;
        `mkdir $tmpdir/tmp_gtf`;
        my %chrid = ();
        foreach (@chrid){
                if (!/_/ && !/^chr.*[^0-9XY]/){
                        `cat $opt_gtf | awk '\$1 ~ /^$_\$/' | awk '\$3 ~ /^exon\$/ {print \$4,\$5,\$7,\$10,\$14,\$16;}'| sort  > $tmpdir/tmp_gtf/$_.gtf`;
                        $chrid{$_} = 1;
                }
        }

        open(IDNEW,">accepted_reads.txt") || die "$!";
        open(ID,"accepted_hits.sorted.sam") || die "$!";
        my $line = '';
        while(<ID>){
                if (!/^@/){
                        $line = $_;
                        last;
                }

        }
        chomp $line;
        my %sam_flag = ('99' => 1, '147' => 1, '83' => 1, '163' => 1);
        while($line){
                if ($line =~ /\s\=\s/){
                        my @array1 = split(' ',$line);
                        if (exists($sam_flag{$array1[1]})){
                                my @array2 = ();
                                while(<ID>){
                                        if (/\s\=\s/){
                                                @array2 = split(' ',$_);
                                                if (exists($sam_flag{$array2[1]})){
                                                        last;
                                                }
                                        }
                                }
                                if ($array1[3] > $array2[3]){
                                        my @array3 = @array2;
                                        @array2 = @array1;
                                        @array1 = @array3;
                                }
                                if ($array1[1] == 147 || $array1[1] == 83){
                                        if ($array1[5] !~ /[^MNID0-9]/  && $array2[5] !~ /[^MNID0-9]/){
                                                my $cigar1 = $array1[5];
                                                my $cigar2 = $array2[5];
                                                $cigar1 =~ s/[0-9]+I//g;
                                                $cigar2 =~ s/[0-9]+I//g;

                                                my @num1 = split(/[MND]/,$cigar1);
                                                my @num2 = split(/[MND]/,$cigar2);

                                                my $t11 = $array1[3];
                                                my $t12 = $t11-1;
                                                my $t21 = $array2[3];
                                                my $t22 = $t21-1;
                                                foreach (@num1){
                                                        $t12 += $_;
                                                }
                                                foreach (@num2){
                                                        $t22 += $_;
                                                }
                                                if ($t12 < $t22){
                                                        if (exists($chrid{$array1[2]})){
                                                                my $t11_l = $t11+2;
                                                                my $t11_r = $t11-2;
                                                                my @exons1 = `cat  $tmpdir/tmp_gtf/$array1[2].gtf | awk '\$1 <= $t11_l' | awk '\$2 >= $t11_r' `;
                                                                chomp @exons1;
                                                                my $t22_l = $t22+2;
                                                                my $t22_r = $t22-2;
                                                                my @exons2 = `cat  $tmpdir/tmp_gtf/$array2[2].gtf | awk '\$1 <= $t22_l' | awk '\$2 >= $t22_r' `;
                                                                chomp @exons2;
                                                                my $exonnum1 = @exons1;
                                                                my $exonnum2 = @exons2;
                                                                if (@exons1 >= 1 && @exons2 >=1){
                                                                        if (@exons1 == 1 && @exons2 == 1) {
                                                                                my @exon_info1 = split(' ',$exons1[0]);
                                                                                my @exon_info2 = split(' ',$exons2[0]);
                                                                                if ($exon_info1[3] eq $exon_info2[3]){
                                                                                        my $inter_len = $t11 - $exon_info1[0] + $exon_info2[1] - $t22;
                                                                                        if ($inter_len < 300){
                                                                                                my @geneid = split(/"/,$exon_info1[3]);
                                                                                                my @transcriptid = split(/"/,$exon_info1[4]);
                                                                                                if($transcriptid[1] !~ /^N/){
                                                                                                        @transcriptid = split(/"/,$exon_info1[5]);
                                                                                                }
                                                                                                my $len = $exon_info2[1] - $exon_info1[0] + 1;
                                                                                                my $strand = $exon_info1[2];
                                                                                                print IDNEW $array1[2],"\t",$exon_info1[0]-1,"\t",$exon_info2[1],"\t",$strand,"\t",$geneid[1],"\t",$len,"\t","1\t",$transcriptid[1],"\n";
                                                                                        }
                                                                                }
                                                                        } elsif (@exons1 > 1 ||  @exons2 > 1){
                                                                                my %exonpos = ();
                                                                                my %exoninfo =();
                                                                                foreach (@exons1){
                                                                                        my @exon_info1 = split(' ',$_);
                                                                                        foreach (@exons2){
                                                                                                my @exon_info2 = split(' ',$_);
                                                                                                if ($exon_info1[3] eq $exon_info2[3]){
                                                                                                        my $inter_len = $t11 - $exon_info1[0] + $exon_info2[1] - $t22;
                                                                                                        if ($inter_len < 300){
                                                                                                                my @geneid = split(/"/,$exon_info1[3]);
                                                                                                                my @transcriptid = split(/"/,$exon_info1[4]);
                                                                                                                if($transcriptid[1] !~ /^N/){
                                                                                                                        @transcriptid = split(/"/,$exon_info1[5]);
                                                                                                                }
                                                                                                                my $len = $exon_info2[1] - $exon_info1[0] + 1;
                                                                                                                my $strand = $exon_info1[2];
                                                                                                                my $left_edge = $exon_info1[0] -1;
                                                                                                                $exonpos{$geneid[1]} = "$left_edge\t$exon_info2[1]\t$strand\t";
                                                                                                                $exoninfo{$geneid[1]} = "$len\t1\t$transcriptid[1]\n";
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                                my @exonids = keys(%exonpos);
                                                                                if (@exonids >= 1){
                                                                                        my @sortids = sort{$a cmp $b}@exonids;
                                                                                        print IDNEW $array1[2],"\t",$exonpos{$sortids[0]},$sortids[0],"\t",$exoninfo{$sortids[0]};
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }
                                        }
                                }
                                $line = readline(ID);
                                if ($line){
                                        chomp $line;
                                }

                        } else {
                                $line = readline(ID);
                                if ($line){
                                        chomp $line;
                                }
                        }

                }else {
                        $line = readline(ID);
                        if ($line){
                                chomp $line;
                        }
                }
        }

        close(ID);
        close(IDNEW);

        if ($opt_temp){
                `rm -rf $tmpdir`;
        }

}



sub cR_statistics
{
	open (IN,"mapped_reads_list.txt");
	open (OUT,">circRNA_list_temp.txt");

	my $line;
	my %hash;
	my $cR_count = 1;
	my @array;
	my $index_string;

	#print OUT join ("\t","NO.","Chromosome","GeneName","JunctionLowerSite","JunctionUperSite","LowerExon","UperExon","ReadStrand","NumberOfMappedReads");
	#print OUT "\n";

	while ($line=<IN>)
	{
		chomp $line;
		@array = split(/\t/, $line);
		
		$index_string = join ("*",$array[1],$array[2],$array[3]);

		if (!exists $hash{$index_string})
		{
			$hash{$index_string} =1; 
		}
		else
		{
			$hash{$index_string} += 1;
		}
	}
	close (IN);

	open (IN,"mapped_reads_list.txt");
	my %list_hash;
	while ( $line=<IN> )
	{
		chomp $line;
		@array = split(/\t/, $line);
		
		$index_string = join ("*",$array[1],$array[2],$array[3]);
		
		if (!exists $list_hash{$index_string})
		{
			$list_hash{$index_string} = join ("\t",$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$hash{$index_string},$array[7]);
		}
	}
	foreach my $key (keys %list_hash)
	{
		print OUT $list_hash{$key},"\n";
	}
	close(IN);
	close(OUT);

	#sorting the circRNA list;
	my $count = 0;
	my @arr;
	my $i;
	@array = ();
	open(IN,"circRNA_list_temp.txt");
	open(OUT,">circRNA_list_temp2.txt");
	$line = <IN>;
	while ($line = <IN>)
	{
		push @array,[split(/\t/,$line)];
		$count += 1;
	}
	@arr = sort{$b->[6] <=> $a->[6]} @array;
	for($i=0; $i<$count;$i++)
	{
		print OUT join("\t",$arr[$i][0],$arr[$i][1],$arr[$i][2],$arr[$i][3],$arr[$i][4],$arr[$i][5],$arr[$i][6],$arr[$i][7]);
	}
	close(IN);
	close(OUT);
#filtering the circRNA list >2 reads supported;
	my @line1;
	open(IN,"circRNA_list_temp2.txt");
	open(OUT,">circRNA_list.txt");
	while($line=<IN>){
	@line1=split(/\t/,$line);
	if($line1[6]>1){
		print OUT $line;
		}
	
	}
close(IN);
close(OUT);
}



#----------------------------------------------------------

sub seq_comp
{
	my $exon_string = $_[0];
	my $read_string = $_[1];
	my @exon_array = split ("",$exon_string);
	my @read_array = split ("",$read_string);
	my $exon_number = @exon_array;
	my $read_number = @read_array;
	
	my $i=0;
	my $mismatch_number = 0;
	
	if ( $exon_number == $read_number ) #checking the length!
	{
		for ( $i = 0; $i < $read_number; $i += 1 ) 
		{
			if ( $exon_array[$i] ne $read_array[$i] )
			{
				$exon_array[$i] =~ tr/atcg/ATCG/;
				if ( $exon_array[$i] ne $read_array[$i] )
				{
					$mismatch_number += 1;
				}
			}
		}
	}
	return $mismatch_number;
}


#fetching exon forward sequence from *.fa file according to the locating site.
sub fetch_Lexon_seq
{
	my $chr_file = join ("/",$opt_fasta,join (".",$_[0],"fa") );#define the directiory of the chromosome.
	
	open (CHR_L, "$chr_file");
	my $pos = $_[1]%50;
	my $line_num = int ( $_[1]/50 );
	my $line;
	my $line1;
	my $line2;
	my $line3;
	my $exon_tmp;
	my $Lower_exon;
	
	my $count=-1;

	while ( $line = <CHR_L> )
	{
		$count += 1;
		if ( $count > $line_num )
		{
			$line1 = <CHR_L>; $line2 = <CHR_L>;$line3 = <CHR_L>;
				
			chomp $line;chomp $line1;chomp $line2;chomp $line3;
				
			$exon_tmp = join ("",$line,$line1,$line2,$line3 );
			if ($pos != 0)
			{
				$Lower_exon = substr($exon_tmp, $pos-1, 123);
			}
			if ($pos == 0)
			{
				$Lower_exon = substr($exon_tmp, $pos, 123);
			}
			last;
		}
	}
	close(CHR_L);
	return($Lower_exon);
}

#fetching exon reverse sequence from *.fa file according to the locating site.
sub fetch_Rexon_seq
{
	my $chr_file = join ("/",$opt_fasta,join (".",$_[0],"fa") );
	open (CHR_U, "$chr_file");
	my $pos = $_[1]%50;
	my $line_num = int ( $_[1]/50 );
	my $line;
	my $line1;
	my $line2;
	my $line3;
	my $exon_tmp;
	my $Uper_exon;
	
	my $count=-1;

	while ( $line = <CHR_U> )
	{
		$count += 1;
		if ( $count > $line_num -3 )
		{
			$line1 = <CHR_U>;$line2 = <CHR_U>;$line3 = <CHR_U>;
		
			chomp $line;chomp $line1;chomp $line2;chomp $line3;
		
			$exon_tmp = join ("",$line,$line1,$line2,$line3 );
			$Uper_exon = reverse ( substr(reverse($exon_tmp), 50-$pos, 123) );
		
			last;
		}
	}
	close(CHR_U);
	return($Uper_exon);
}



