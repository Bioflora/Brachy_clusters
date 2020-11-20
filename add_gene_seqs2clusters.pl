#!/usr/bin/perl
use strict;
use warnings;

# Takes as input 
# 1) precomputed FASTA files with transcript clusters <input cluster dir>
# 2) new sequences to be added to those clusters
#  2.1) a FASTA file per new species to be added (%references)
#  2.2) a folder with pre-computed BLASTN searches matching reference seqs to clusters
#
# In the case of Bmex/Bhyb two alleles/copies are expected for being tetraploid.
# When two sequences are found they are ranked by similarity to defined ancestral and 
# recent species
#
# Produces enlarged clusters with new added/substituted sequences in <output dir>:
# i) when a reference sequence is found in the original cluster, it is replaced
#    with a new one from the reference FASTA file
# ii) if a reference species was not present in the original cluster it is added
#
# MaDecena Decena, Bruno Contreras Moreira EPS-UNIZAR & EEAD-CSIC 2020

if(!$ARGV[1]){
	print "# usage: $0 <input cluster dir> <output dir>\n\n";
	print "# Please make sure you configure the following:\n";
	print "# %blastdirs, %references, \@order, %max_references_seqs\n\n";
	exit(0);
}

my ($transcriptdir, $finalclusterdir) = @ARGV;

# make sure paths can be concatenated
$transcriptdir   .= '/';
$finalclusterdir .= '/';

print "# params:\n";
print "# input cluster dir: $transcriptdir\n";
print "# output dir: $finalclusterdir\n";

my %blastdirs = (
	'Bdis' => 'blast_Bdis/',
	'Bsta' => 'blast_Bsta/',
	'Bsyl' => 'blast_Bsyl/',
	'Barb' => 'blast_Barb/',
	'Bmex' => 'blast_Bmex/',
	'Bpin' => 'blast_Bpin/',
	'Bgla' => 'blast_Bgla/',
	'Bhyb' => 'blast_Bhyb/',
	'Osat' => 'blast_Osat/',
	'Hvul' => 'blast_Hvul/',
);

my %references = (
        'Bdis' => 'referencias/Bdistachyon.unspliced.gene.100flank.fna',
        'Bsta' => 'referencias/Bstacei.unspliced.gene.100flank.fna',
        'Bsyl' => 'referencias/Bsylvaticum.unspliced.gene.100flank.fna',
	'Barb' => 'referencias/Barbuscula_genes.nr.fna',
	'Bmex' => 'referencias/Bmexicanum_genes.nr.fna',
	'Bpin' => 'referencias/Bpin_consensusG_shortname.1L.depth4min1000_nomissing.fna',
	'Bgla'=> 'referencias/Bglau_consensusG_shortname.1L.depth4.min1000.nomissing.fna',
	'Bhyb'=> 'referencias/Bhybridum.unspliced.gene.100flank.fna',
	'Osat' => 'referencias/Osativa.unspliced.gene.100flank.fna',
	'Hvul' => 'referencias/Hvulgare.unspliced.gene.100flank.fna',
);

# any species not in the list will be ignored
# Note: Bmex, Bhyb are not diploids
my @order = ( 'Osat', 'Hvul', 'Bmex', 'Bsta', 'Bhyb', 'Bdis', 
		'Barb', 'Bpin', 'Bsyl', 'Bgla' );

# By default only the 1st, the best hit, is taken.
# However, for polyploids or species absent from original clusters,
# this won't work.
# Two species names are provided that represent ancestral and recent species,
# respectively. These are used to rank multiple alleles
my %max_reference_seqs = (
	'Bmex' => [ 
		2,      # expected copies/alleles 
		'Hvul', # ancestral species, usually an outgroup
		'Bsta', # recent species, expected to descend from ancestral in history
		2,      # max allowed alleles, as assembly does not resolve subgenomes
	],
	'Bhyb' => [
		2,      
		'Bsta',
		'Bdis',
		0	# max not defined, as we expect exactly 1 D & 1 S allele 
	],
	'Bgla' => [
		1,
		'Bsta',
		'Bsta',
		0
	]
);

printf("# blastdirs: %s\n",join(',',sort keys(%blastdirs)));
printf("# references: %s\n",join(',',sort values(%references)));
printf("# order: %s\n",join(',',@order));
printf("# max_reference_seqs: %s\n\n",join(',',sort keys(%max_reference_seqs)));

##########################################################

# make sure missing BLASTN outfiles are regenerated, in case original jobs failed
my $BLASTNCMD = '/home/soft/ncbi-blast-2.10.0+/bin/blastn -outfmt 6 '.	
		'-max_target_seqs 10 -perc_identity 90 ';


my ($sp,$seqid,%ref_sequences);

## read in reference sequences ahead
foreach $sp (keys(%references)){
	open(FASTA,"<",$references{$sp}) || die "# ERROR: cannot read $references{$sp}\n";
        while(my $line = <FASTA>){
                chomp($line);
                if($line =~ m/^>(\S+)/){ $seqid = $1 }
                else{ $ref_sequences{$sp}{$seqid} .= lc($line) }
        }
        close(FASTA);
}

## loop one cluster at a time
opendir(INFOLDER,$transcriptdir) || die "# ERROR: cannot list $transcriptdir\n";
my @infiles = grep {/.fna/} readdir(INFOLDER);
closedir(INFOLDER);

FILE: foreach my $file (@infiles) {

	# debugging
        #next if ($file !~ '99843_c28760_g1_i1_chr3' && $file !~ '93618_c36457_g1_i2_chr2');
	#next if ($file !~ '99843_c28760_g1_i1_chr3');  # 2 Bmex seqs in the original, 4 alleles. skipped 
	#next if ($file !~ '93140_c31422_g1_i5_chr4');  # 1 Bmex seq in the original, 1+dummy at final_cluster	
	#next if ($file !~ '99835_c35434_g2_i3_chr8');  # 2 Bmex seqs in the original, 4 alleles. skipped
	#next if ($file !~ '92929_c38893_g1_i5_chr1');  # 0 Bmex seqs in the original, 2 new at final_cluster
        #next if ($file !~ '93820_c38673_g1_i9_chr5');  # 1 Bmex seq in the original, 1+dummy at final_cluster. 8Oct skipped
	#next if ($file !~ '94018_c35272_g1_i1_chr6');  # 1 Bmex seq in the original, 2 at the final_cluster. 8Oct 1+dummy 
	#next if ($file !~ '96277_c29024_g1_i1_Chr02'); # 1 Bmex seq in the original, 1+dummy at the final cluster. 8Oct only 1 Bmex	
	#next if ($file !~ '104363_c38290_g4_i2_chr2'); # 3 Bmex seq in the original, 2 at the final cluster. 8oct skipped
	#next if ($file !~ '95103_c34081_g1_i3_chr5');  # Same Bit-score, different sequences
        #next if ($file !~ '103577_c35382_g1_i2_chr4'); # 1 Bmex seq in the original, 1 at the final cluster, no blast hits, should be discarded 
	#next if ($file !~ '93618_c36457_g1_i2_chr2');  # 1 Bmex seq in the original, 2 hits, 1 dummy at the final cluster. 8Oct 2Bmex at final cluster.
	#next if ($file !~ '98148_c37636_g1_i2_chr9');  # 8Oct.0 Bmex seqs in the original, 2 sequences at the final cluster. 
	# next if ($file !~ '93023_c37202_g2_i1_chr5'); # Not at 03_blocks_genes 08Oct
	#next if ($file !~ '102250_c37694_g3_i4_chr9');
	next if ($file !~ '108368_c31249_g2_i1_chr1');

	print "$file\n";

	my ($name,$shortsp,%seqs,%genes,@sorted_ids);

	# read transcript sequences
	open(FASTA,"<",$transcriptdir.$file) || 
		die "# ERROR: cannot read $transcriptdir.$file\n";
	while(my $line = <FASTA>){
   		chomp($line);
	
		if($line =~ m/^>(\S+)/){ 
			$seqid = $1; 
			push(@sorted_ids,$seqid);
		}
		else{ $seqs{$seqid} .= lc($line) }
	}
	close(FASTA);

	# read BLAST results for these transcripts (against references)
	# and sort alleles
	my ($qseqid,$sseqid,$pident,$len,$mism,$gap);
	my ($qstart,$qend,$sstart,$send,$eval,$bitscore);
	foreach $sp (sort(keys(%blastdirs))){
		
		# produce BLASTN results if original cannot be found
		if(!-e $blastdirs{$sp}.$file){
			my $cmd = $BLASTNCMD .
				"-out $blastdirs{$sp}$file ".
				"-query $transcriptdir$file ".
				"-db $references{$sp} ";

			system("$cmd");
		} 

		# rank hits for this species ($sp) by parsing 
		# BLASTN results and tallying occurrences and bit-score
      # next if($sp ne 'Bgla');
		my (%stats,%sp_stats,%queries,%bit,%ancestral_bit,@alleles);
		if(open(BLAST,"<",$blastdirs{$sp}.$file)){
			while(my $line = <BLAST>){
				chomp ($line);

				# qseqid -> transcript-name_Bxxx
				# sseqid -> name of reference sequence in blastdirs{$sp}
				($qseqid,$sseqid,$pident,$len,$mism,$gap,
				$qstart,$qend,$sstart,$send,$eval,$bitscore) = 
					split(/\t/,$line);

				$shortsp = '';
				if($qseqid =~ /_(\w{4})$/){ $shortsp = $1 }

				# acumulate scores for reference sequences for $sp
				if($shortsp eq $sp){
					$stats{$sseqid}{'tot'}++;
					$stats{$sseqid}{'bit'}+=$bitscore;
					#$queries{$qseqid}++;
				}
                                
            # save scores to other species as well, 
				# useful to compute distances and ancestry
				$sp_stats{$shortsp}{$sseqid}{'tot'}++;
				$sp_stats{$shortsp}{$sseqid}{'bit'}+=$bitscore;

				# init %bit and %ancestral_bit
				$bit{$sseqid}{'ancestral'} = 0;
				$bit{$sseqid}{'recent'} = 0;
				$ancestral_bit{$sseqid}	= 0;
			}
			close(BLAST);

			if(scalar(keys(%stats)) == 0) {
				print "# no blast hits for $sp\n";

				if($max_reference_seqs{$sp}){
					# take best hits to other species as best hits for $sp
					print "# using hits from previous species\n";
					foreach my $sp2 (keys(%sp_stats)){
						foreach $sseqid (keys(%{ $sp_stats{$sp2} })){
							$stats{$sseqid}{'tot'} += 
								$sp_stats{$sp2}{$sseqid}{'tot'};
                    		$stats{$sseqid}{'bit'} += 
									$sp_stats{$sp2}{$sseqid}{'bit'};
						}
					} 
				}
			}
		} else {
			die "# ERROR: cannot find $blastdirs{$sp}$file, ".
				"use original available transcript\n";
		} 

		# if required, sort alleles as ancestral/recent
		# Note: resulting cumulative bitscore measure similarity to recent,
		# so alleles can sorted from small to large bit scores, in other words,
		# from ancestral to recent
		if($max_reference_seqs{$sp}){

			foreach my $sp2 (keys(%sp_stats)){

				# substract similarity to ancestor [1] species
				if($sp2 eq $max_reference_seqs{$sp}[1]){
					foreach $sseqid (keys(%{$sp_stats{$sp2}})){
						# bit score when aligned to ancestral
						$bit{$sseqid}{'ancestral'} += 
							$sp_stats{$sp2}{$sseqid}{'bit'};
						# in order to produce a single score
						# to sort alleles
						$ancestral_bit{$sseqid} -= 
							$sp_stats{$sp2}{$sseqid}{'bit'};

					}
				} # add similarity to recent [2] species
				elsif($sp2 eq $max_reference_seqs{$sp}[2]){
					foreach $sseqid (keys(%{$sp_stats{$sp2}})){
						# bit score when aligned to recent
						$bit{$sseqid}{'recent'} +=
                                                        $sp_stats{$sp2}{$sseqid}{'bit'};
						# in order to produce a single score
                                                # to sort alleles
                                                $ancestral_bit{$sseqid} +=
                                                        $sp_stats{$sp2}{$sseqid}{'bit'}
                                        }
				}
			}

			#debugging
			#foreach $sseqid (keys(%ancestral_bit)){
                        #               print "$sseqid $ancestral_bit{$sseqid}\n";
                        #}

			## uncoment if needed
			## print number of Bmex alleles for report
			#printf("%s\t%d\t%s\n",
			#	$file,
			#	scalar(keys(%stats)),
			#	join(',',keys(%stats)));
			#next FILE;
				
			# count alleles and stop if too many
			if($max_reference_seqs{$sp}[3] > 0 && #only if defined as in Bmex
				scalar(keys(%stats)) > $max_reference_seqs{$sp}[3]){
				foreach $sseqid (keys(%stats)){
				print "$sseqid $stats{$sseqid}{'tot'} $stats{$sseqid}{'bit'}\n";
				}
				print "# WARNING; skip cluster due to too many $sp alleles\n";
                                next FILE;
			}
		}

   		# choose best hit(s)
   		my ($besthit, $hits) = ('',0);	
   		foreach $sseqid (sort {$stats{$b}{'tot'}<=>$stats{$a}{'tot'} ||
			$stats{$b}{'bit'}<=>$stats{$a}{'bit'} } keys(%stats)){
			
			$besthit = $sseqid;
			$hits++;

			#debugging	
			#print "$sp $besthit $stats{$besthit}{'tot'} $stats{$besthit}{'bit'}\n";
			#print "$ref_sequences{$sp}{$besthit}\n";

			# species where a single copy is expected
			if(!defined($max_reference_seqs{$sp})){
				push(@{$genes{$sp}{'name'}}, (split(/\|/,$besthit))[0]);
				push(@{$genes{$sp}{'sequence'}}, 
					$ref_sequences{$sp}{$besthit});
				last;
			}
			else {  # manage species where multiple copies/hits are expected                 
				# Note: in such cases alleles will be added to genes list
				# from more ancestral to recent order

				# add copies/alleles to temp list
				if($hits <= $max_reference_seqs{$sp}[0]){
					push(@alleles, $besthit); #print "$hits $besthit\n";
				} 

				# when max alleles or last hit is found  
				if($hits == $max_reference_seqs{$sp}[0] || 
					$hits == scalar(keys(%stats)) || 
					$hits == scalar(keys(%{$sp_stats{$sp}}))){

					$name = (split(/\|/,$besthit))[0];

					foreach my $allele (sort {$ancestral_bit{$a} <=> 
						$ancestral_bit{$b}} @alleles){
			
						printf("# ancestral bitscore: ".
							"%s %s %1.1f %1.1f %1.1f\n",
							$allele, $sp,
							$ancestral_bit{$allele},
							$bit{$allele}{'ancestral'},
							$bit{$allele}{'recent'});

						$name = (split(/\|/,$allele))[0];	
						push(@{$genes{$sp}{'name'}},$name);
						push(@{$genes{$sp}{'sequence'}}, 
							$ref_sequences{$sp}{$allele});
					}
				

					# add dummy allele if less than max copies found,
					# considering bit score to ancestral and recent
					# sequences
					if($hits < $max_reference_seqs{$sp}[0]){ 
						
						if($bit{$sseqid}{'ancestral'} >=
							$bit{$sseqid}{'recent'} ) {
							push(@{$genes{$sp}{'name'}}, 
								'dummy');
							push(@{$genes{$sp}{'sequence'}}, 
								'-' x 100);
						} else {
							unshift(@{$genes{$sp}{'name'}}, 
								'dummy');
                     unshift(@{$genes{$sp}{'sequence'}}, 
								'-' x 100);
						}
					}

					last;
				}
			}
		} 
	} 

			 
	# create updated cluster with (complete) sequence for selected species
	my ($keep_transcript,%species_seen,%ref_subs);
	open(GENECLUSTER,">",$finalclusterdir.$file) || 
		die "# cannot write to $finalclusterdir$file\n";

	foreach my $sp (@order) {

	foreach $seqid (@sorted_ids){

		$shortsp = (split(/_/,$seqid))[-1];
		next if($shortsp ne $sp);

		$species_seen{ $shortsp }++;

		$keep_transcript = 1;

		# if there are reference sequences for this species see if they
		# are longer and can thus substitute transcripts from the original cluster
		# Note: this is performed only with the 1st sequence of each species 
		if($genes{$shortsp} && $species_seen{$shortsp} == 1){

			my ($ref_length, $ref_total, $ref_mean) = (0, 0, 0);

			# compute mean length of ref sequence
			foreach my $hit (0 .. scalar(@{$genes{$shortsp}{'name'}})-1){

				$ref_length += length($genes{$shortsp}{'sequence'}[$hit]);
				$ref_total++;
			}
			$ref_mean = sprintf("%1.0f",$ref_length/$ref_total);

			# take always ref sequences for species in %max_reference_seqs,
			# otherwise make sure ref sequences are not shorter than transcripts
			if($max_reference_seqs{$sp} ||
				length($seqs{$seqid}) < $ref_mean){

				$keep_transcript = 0;
				$ref_subs{ $shortsp } = 1;

				# actually add ref sequences to cluster
				foreach my $hit (0 .. scalar(@{$genes{$shortsp}{'name'}})-1){

					if($hit == 0){
                                	print GENECLUSTER ">$shortsp $genes{$shortsp}{'name'}[$hit]\n";
					} else { # number extra alleles
						printf(GENECLUSTER ">%s.%d %s\n",
							$shortsp,
							$hit+1, 
							$genes{$shortsp}{'name'}[$hit]);
					}


					print GENECLUSTER "$genes{$shortsp}{'sequence'}[$hit]\n";
            
					if($genes{$shortsp}{'name'}[$hit] !~ 'dummy'){ 
						printf("%s\t%d > %d\n",$genes{$shortsp}{'name'}[$hit],
                    		length($genes{$shortsp}{'sequence'}[$hit]),
                    		length($seqs{$seqid}));
					} else {
						printf("%s\t%d\n",$genes{$shortsp}{'name'}[$hit],
							length($genes{$shortsp}{'sequence'}[$hit]));
					}
                        	}
		
			} else {
				print "$shortsp ref sequences shorter than transcript\n"
			}		
		} 

		# add original transcript instead
		if($keep_transcript == 1 && !$ref_subs{$shortsp}) { 
			print GENECLUSTER ">$shortsp $seqid\n$seqs{$seqid}\n";
		}
	}

	# add dummy gapped sequence for missing species/alleles
	if(!defined($species_seen{$sp})){

		if(!defined($max_reference_seqs{$sp})){
			print GENECLUSTER ">$sp missing\n";
			print GENECLUSTER '-' x 100;
			print GENECLUSTER "\n";	
		} 
		else {
			# recorre lista de genes de species en %max_reference_seqs y añadirlos al cluster
			  foreach my $hit (0 .. scalar(@{$genes{$sp}{'name'}})-1){
				print GENECLUSTER ">$sp $genes{$sp}{'name'}[$hit]\n";
				print GENECLUSTER "$genes{$sp}{'sequence'}[$hit]\n";
			}
		}
	} 
		
	}
	close(GENECLUSTER);

} # foreach
