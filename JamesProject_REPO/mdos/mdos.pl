#!/usr/bin/perl -w
use lib './util';
use gene_utils;

die "usage: perl $0 input.seq kmer_size single_strand(1:single,0:both) with_extension(1:yes,0:no) num_samples" if ($#ARGV<4);


# use single or both strand
$USE_SINGLE_STRAND = $ARGV[2];
# whether to expand to degenerate code
$WITH_IUB_EXTENSION = $ARGV[3];
# number of sampling iterations
$SAMPLE_ITER_NUM = $ARGV[4];


# size of kmers
$kmer_size = $ARGV[1];
# z-score threshold for expanding
$ZSCORE_TH = 3;
# number of degenerate bases to expand
$NUM_DEG_BASE = 2;
# minimum number of shared genes for testing
$MIN_NUM_SHARED_GENES = 5;



# read in sequence list
@seqList = ();
@all = get_file_data($ARGV[0]); chomp @all;
for ($i=0; $i<@all; $i++) {
    next if (not $all[$i]=~/^>/);
    @a = split(" ", $all[$i]);
    $id = $a[0]; $id=~s/^>//;

    $seqA = $all[++$i];
    $seqB = $all[++$i];
  
    if ($USE_SINGLE_STRAND < 1) {  # if consider both strand
	$rc = reverse_comp($seqA);
	$seqA .= ('x' . $rc);
	$rc = reverse_comp($seqB);
	$seqB .= ('x' . $rc);
    }

    $numB = 0;
    $L = length($seqB)-$kmer_size;
    for ($k=0; $k<=$L; $k++) {
        $tmp = substr($seqB, $k, $kmer_size); 
	next if ($tmp=~/[^ACGT]/);   # skip if anything except ACGT
	$numB++;
    }
    
    push @seqList, {'seqA'=>$seqA, 'seqB'=>$seqB, 'numB'=>$numB, 'id'=>$id};
}




# generate kmer hash table,store the index of sequences, for both sequences
%kmerPosA=(); %kmerPosB=();
for ($i=0; $i<@seqList; $i++) {
    $seqA = $seqList[$i]->{'seqA'};
    $seqB = $seqList[$i]->{'seqB'};

    # for seqA
    %foo = ();   $L = length($seqA)-$kmer_size;
    for ($k=0; $k<=$L; $k++) {
	$tmp = substr($seqA, $k, $kmer_size);
	next if ($tmp=~/[^ACGT]/);   # skip if anything except ACGT
	$foo{$tmp}=1;
    }
    foreach (keys %foo) {
	push @{$kmerPosA{$_}}, $i;
    }
    
    # for seqB
    %foo = ();   $L = length($seqB)-$kmer_size;
    for ($k=0; $k<=$L; $k++) {
	$tmp = substr($seqB, $k, $kmer_size);
	next if ($tmp=~/[^ACGT]/);   # skip if anything except ACGT
	$foo{$tmp}=1;
    }
    # for all kmer in the seq
    foreach (keys %foo) {
	push @{$kmerPosB{$_}}, $i;
    }
}


# get candidate kmer list
# get kmer count stat based on sequence A
%count_conserv_AB =();
foreach $m (keys %kmerPosA) {
    %fooA=(); %fooB=();
    foreach (@{$kmerPosA{$m}}) {
	$fooA{$_}=1;
    }
    if (exists $kmerPosB{$m}) {
	foreach (@{$kmerPosB{$m}}) {
	    $fooB{$_}=1;
	}
    }
    $n = 0;
    foreach (keys %fooA) {
	$n++ if (exists $fooB{$_});
    }
    $count_conserv_AB{$m} = $n;
}

# get conservation rate based on hypergeometric
%motif2HG = ();
$N = scalar(@seqList);
#print "#hyergeometric stat\n";    
foreach $m (keys %count_conserv_AB) {
    next if ( ($m=~tr/[GC]//) < 1);
    next if ($count_conserv_AB{$m} < 5);
    $nA = scalar(@{$kmerPosA{$m}});
    $nB = (exists $kmerPosB{$m} ? scalar(@{$kmerPosB{$m}}) : 0);
    $mu = $nA*$nB/$N;
    $sig = sqrt($nA*$nB/$N*(1-$nB/$N)*($N-$nA)/($N-1));
    next if ($sig <= 0);
    $z = ($count_conserv_AB{$m}-$mu)/$sig;
    $motif2HG{$m} = $z;
    #print "#HG ", join("\t", ($m, $nA, $nB, $count_conserv_AB{$m}, $mu, $sig, $z)),"\n";
}


@t = keys %motif2HG;
@t = sort {$motif2HG{$b}<=>$motif2HG{$a}} @t;


# set neutral rate based on gc content, store rate of occurrence in seqB in seqList
set_neutral_rate_for_kmer();


# to get a list of kmers for testing
my %fooHash=();
foreach $m (keys %count_conserv_AB) {
    next if ( ($m=~tr/[GC]//) < 1);
    next if ($count_conserv_AB{$m} < $MIN_NUM_SHARED_GENES);

    if ($USE_SINGLE_STRAND != 1) { # if use both strand
	$rc = reverse_comp($m);
	next if (exists $fooHash{$rc});
    }
    $fooHash{$m}=1;
}
# used kmer list for testing
@kmerList = keys %fooHash; %fooHash=();
@kmerList = sort {$count_conserv_AB{$b}<=>$count_conserv_AB{$a}} @kmerList;


# used global variable to store probability for each sequence
@probListB=(); $#probListB =$#seqList;

# screen kmer first, identify those above z-score threshold
print "#Kmer conservation z-score\n";
print "#Kmer\tnumA\tnumB\tnumAB\tmean\tstd\tz-score\n";
@screened_kmer_list = ();
foreach $kmer (@kmerList) {
    ($totalNumA, $totalNumB, $numAB, $mu, $sig, $zscore) = sample_given_kmer_list($kmer);

    print join("\t", ($kmer, $totalNumA, $totalNumB, $numAB, $mu, $sig, $zscore)),"\n";
    
    next if ($zscore < $ZSCORE_TH);
    push @used_kmer_list, {'motif'=>$kmer, 'numA'=>$totalNumA, 'numB'=>$totalNumB, 'numAB'=>$numAB, 'mu'=>$mu, 'sig'=>$sig,'z'=>$zscore};
}



if ($WITH_IUB_EXTENSION < 1) {
    #print "#with no IUB extension\n";
    exit;
}



# sort the returned list based on z score
@screened_kmer_list = sort {$b->{'z'}<=>$a->{'z'}} @used_kmer_list;


# expand to degenerate codes
print "\n#start to expand to degenerate code ----\n";

# used degnerate code for testing
%base2IUB=('A'=>'RMWN', 'C'=>'YMSN', 'G'=>'RKSN', 'T'=>'YKWN');


foreach $m (@screened_kmer_list) {
    $kmer = $m->{'motif'};
    ($numA, $numB, $numAB, $mu, $sig, $zscore) = ($m->{'numA'},$m->{'numB'},$m->{'numAB'},$m->{'mu'},$m->{'sig'},$m->{'z'});

    print ">$kmer\t", join("\t",($numA, $numB, $numAB, $mu, $sig, $zscore)),"\n";


    $numDegCode = $NUM_DEG_BASE;  %used_pos_to_IUB=(); 
    @prevKmerList=(); push @prevKmerList, $kmer;

    while ($numDegCode-- > 0) {
	
	@statList = (); 
	push @statList, {'IUB'=>'', 'pos'=>-1,'numA'=>$numA,'numB'=>$numB,'numAB'=>$numAB,'sig'=>$sig,'mu'=>$mu,'zscore'=>$zscore,'kmerList'=>join("\t",@prevKmerList)};
    
	for ($i=0; $i<$kmer_size; $i++) {   # try degenerate code at each position 
	    next if (exists $used_pos_to_IUB{$i});

	    $b = substr($kmer, $i, 1);
	    next if (not exists $base2IUB{$b});
	    @IUB_list = split("", $base2IUB{$b});
	
	    
	    foreach $IUB (@IUB_list) {
		$regexp = IUB_to_regexp($IUB);
		$regexp=~s/[^ACGT]//g;
		@baseList = split("", $regexp);
		
		%fooHash=(); # store kmer to be used
		foreach $prevKmer (@prevKmerList) {
		    foreach (@baseList) {
			$f = $prevKmer;
			substr($f, $i, 1, $_);  # substitute at position i
			next if ( ($f=~tr/[ACGT]//) < 1); # at least 1 CG
			$fooHash{$f}=1;
		    }
		}
		@t = keys %fooHash;
		next if (scalar(@t) <= scalar(@prevKmerList));
		
		($numA, $numB, $numAB, $mu, $sig, $zscore) = sample_given_kmer_list(@t);
		push @statList, {'IUB'=>$IUB, 'pos'=>$i, 'numA'=>$numA,'numB'=>$numB,'numAB'=>$numAB,'sig'=>$sig,'mu'=>$mu,'zscore'=>$zscore, 'kmerList'=>join("\t",@t)};
	    }
	}
	
	@statList = sort {$b->{'zscore'}<=>$a->{'zscore'}} @statList;
	($numA,$numB,$numAB,$mu,$sig,$zscore)=($statList[0]->{'numA'},$statList[0]->{'numB'},$statList[0]->{'numAB'},$statList[0]->{'mu'},$statList[0]->{'sig'},$statList[0]->{'zscore'});
	
	last if ($statList[0]->{'pos'} < 0);
	
	$used_pos_to_IUB{$statList[0]->{'pos'}} = $statList[0]->{'IUB'};
	@prevKmerList = split(" ", $statList[0]->{'kmerList'});


	#print "# ", join("\t", ($numDegCode, $statList[0]->{'pos'}, $statList[0]->{'IUB'})),"\n";
	#print "# ", join("\t", ($numA, $numB, $numAB, $mu, $sig, $zscore)),"\n";
    }
   
    # get new consensus
    $new_motif = $kmer;
    foreach (keys %used_pos_to_IUB) {
	substr($new_motif, $_, 1, $used_pos_to_IUB{$_});
    }

    print "$new_motif\t", join("\t",($numA, $numB, $numAB, $mu, $sig, $zscore)),"\n";
    
    print "\n";
}






sub sample_given_kmer_list {
    my (@kmerList) = @_;
    
    my %fooA=(); my %fooB=(); my ($m, $gc, $k) = ();

    # calculate neutral matching prob for each sequence
    foreach (@probListB) {$_=0;} #probListB: global variable

    foreach $m (@kmerList) {
	$gc = ($m=~tr/[GC]//);
	%fooA=(); # store index of seqA
	if (exists $kmerPosA{$m}) {
	    foreach (@{$kmerPosA{$m}}) {
		$fooA{$_}=1;
	    }
	}
	for ($k=0; $k<@probListB; $k++) {
	    $probListB[$k] += (exists $fooA{$k} ? $seqList[$k]->{'neutral_bg_rate_with_motif_inA'}->[$gc] :  $seqList[$k]->{'neutral_bg_rate'}->[$gc]);
	}
    }
	
    
    # calculate number
    %fooA=(); %fooB=();
    foreach $m (@kmerList) {
	if (exists $kmerPosA{$m}) {
	    foreach (@{$kmerPosA{$m}}) {
		$fooA{$_}=1;
	    }
	}
	if (exists $kmerPosB{$m}) {
	    foreach (@{$kmerPosB{$m}}) {
		$fooB{$_}=1;
	    }
	}
    }
    
    my $totalNumA = scalar(keys %fooA);
    my $totalNumB = scalar(keys %fooB);
    my $numAB = 0;
    foreach (keys %fooA) {
	$numAB++ if (exists $fooB{$_});
    }

    # normalize prob
    my $tmp = sum_given_ref(\@probListB);
    die "Error in probListB: $tmp" if ($tmp <= 0);
    foreach (@probListB) {
	$_ = $_/$tmp;
    }

    die "Error: uneqaul size in probListB" if ($#seqList != $#probListB);
    

    # batch sample-based method
    my $sample_iter = $SAMPLE_ITER_NUM;
    my @sampled_conserv_num=(); my @sampled_list = (); my $n='';
    while ($sample_iter-->0) {
	# generate sampled list
	@sampled_list = generate_sampled_list_batch($totalNumB, \@probListB);
	$n = 0;
	foreach (@sampled_list) {
	    $n++ if (exists $fooA{$_});
	}
	push @sampled_conserv_num, $n;
	#print join("\t", ($totalNumA, scalar(@sampled_list), $numAB, $n)),"\n";
    }
    
    my $mu = mean_given_ref(\@sampled_conserv_num);
    my $sig = std_given_ref(\@sampled_conserv_num);
    my $zscore = ($sig>0 ? ($numAB-$mu)/$sig : 100);
    $mu = sprintf('%.3f', $mu);
    $sig = sprintf('%.3f', $sig);
    $zscore = sprintf('%.3f', $zscore);

    return ($totalNumA, $totalNumB, $numAB, $mu, $sig, $zscore);
}










# generate sampled index from a list of probabilities, without replacement
sub generate_sampled_list_batch {
    my ($n, $probListRef) = @_;
    use strict;
    use warnings;
    
    # get prob:index pair
    my @probIndex=(); $#probIndex=scalar(@{$probListRef})-1;
    for (my $k=0; $k<@probIndex; $k++) {
	$probIndex[$k] = {'prob'=>$probListRef->[$k], 'index'=>$k};
    }

    my @rand_num=(); my @t=(); my @cumProbIndex=();
    my ($prev_index, $sum)=(); my @indexList=();

    while ($n>0) {
	# get cumulative prob
	@cumProbIndex=(); $#cumProbIndex=$#probIndex-1;
	$sum=0;
	for (my $k=0; $k<$#probIndex; $k++) {  # forget the last one
	    $sum += $probIndex[$k]->{'prob'};
	    $cumProbIndex[$k] = {'cumsum'=>$sum, 'index'=>$k+1}; # index of the next one
	}
	
	@rand_num = (); $#rand_num=$n-1; # generate n random number
	for (my $k=0; $k<$n; $k++) {
	    $rand_num[$k] = {'cumsum'=>rand()*$sum, 'index'=>-1};
	}
	
	push @rand_num, @cumProbIndex;  # combine rand number and cumsum prob
	@rand_num = sort {$a->{'cumsum'}<=>$b->{'cumsum'}} @rand_num;
	
	$prev_index = 0; @t=();
	for (my $k=0; $k<@rand_num; $k++) {
	    if ($rand_num[$k]->{'index'} < 0) {  # it is a random number now
		next if ($k>0 and $rand_num[$k-1]->{'index'} < 0);  # if it has alread taken
		push @t, $prev_index;
	    } else {
		$prev_index = $rand_num[$k]->{'index'};
	    }
	}

	# collect sampled index (map to the original ones)
	foreach (@t) {
	    push @indexList, $probIndex[$_]->{'index'};
	}
	$n -= scalar(@t);
	@t = sort {$b<=>$a} @t;  # remove already sampled prob list
	foreach (@t) { 
	    splice @probIndex, $_, 1;
	}
    }
    
    return @indexList;
}



sub get_optimal_conserv_rate {
    my $countA=0;
    my $countB=0;
    my $k='';
    foreach (@count_conserv) {
	$countA += $_;
    }
    foreach (@count_non_conserv) {
	$countB += $_;
    }
    my $pc = ($countA+$countB>0 ? $countA/($countA+$countB) : 0);  # conserv rate without correcting for background rate
    my $gamma = 1-$pc;
    $gamma = 1-0.00001 if ($gamma>=1); # make it robust
    $gamma = 0.00001 if ($gamma<=0); # make it robust

    my $logL = 0;
    my $deriv_logL=0;
    my $second_deriv_logL = 0;
    for (my $j=0; $j<@listZ; $j++) {
	$logL += ($count_conserv[$j]*log(1-$gamma*$listZ[$j]) +  $count_non_conserv[$j]*log($gamma) );
	$deriv_logL += (-$count_conserv[$j]/(1/$listZ[$j]-$gamma) + $count_non_conserv[$j]/$gamma);
	$second_deriv_logL += (-$count_conserv[$j]/((1/$listZ[$j]-$gamma)**2) - $count_non_conserv[$j]/($gamma*$gamma));
    }

    #print "#start iter: $pc\t", $countA, "\t", $countA+$countB,"\n";

    my $iter = 0;
    while (abs($deriv_logL) > 0.0001 and $iter<100) {
	$iter++;
	$gamma -= $deriv_logL/$second_deriv_logL;
	$gamma=0.0000001 if ($gamma<0);
	$gamma=1-0.000001 if ($gamma>1);
	$logL = 0;
	$deriv_logL=0;
	$second_deriv_logL = 0;
	for (my $j=0; $j<@listZ; $j++) {
	    $logL += ($count_conserv[$j]*log(1-$gamma*$listZ[$j]) +  $count_non_conserv[$j]*log($gamma) );
	    $deriv_logL += (-$count_conserv[$j]/(1/$listZ[$j]-$gamma) + $count_non_conserv[$j]/$gamma);
	    $second_deriv_logL += (-$count_conserv[$j]/((1/$listZ[$j]-$gamma)**2) - $count_non_conserv[$j]/($gamma*$gamma));
	}
	#print "#iter $iter: ", 1-$gamma, "\t", $logL, "\t", $deriv_logL, "\n";
    }
    
    return 1-$gamma;
}
    
    

sub set_neutral_rate_for_kmer {

    # get the number of kmers
    my %count_total = ();  # total number in seqB
    my ($s, $seqA, $seqB, $L, $tmp, $m, $k, $j) = ();
    my %fooHash=(); my @tmpA=(); my @tmpB=(); my @listA=();
    my @gcStat=();

    my $total_size = 0;
    foreach (@seqList) {
	$total_size += $_->{'numB'};
	$seqB = $_->{'seqB'};
	$L = length($seqB)-$kmer_size;
	for ($k=0; $k<=$L; $k++) {
	    $tmp = substr($seqB, $k, $kmer_size);
	    next if ($tmp=~/[^ACGT]/);   # skip if anything except ACGT
	    $count_total{$tmp}++;
	}
    }
    # collect all previous stat
    my %kmerStat = ();
    foreach $kmer (keys %count_total) {
	$kmerStat{$kmer} = {'rand_count'=>$count_total{$kmer}, 'rand_size'=>$total_size};
    }
    
    # use array to store GC contenct combined stat, count of GC bases, from 0 to kmer_size
    @gcStat=(); $#gcStat=$kmer_size; 
    foreach $k (keys %kmerStat) {
	$gcNum = ($k=~tr/[GC]//);
	$gcStat[$gcNum]->{'rand_count'} += $kmerStat{$k}->{'rand_count'};
	$gcStat[$gcNum]->{'rand_size'} += $kmerStat{$k}->{'rand_size'};
    }
    # calculate bg rate, with gc motif separate 
    $k = 0;
    foreach (@gcStat) {
	$_->{'rand_rate'} = ($_->{'rand_size'}>0 ? $_->{'rand_count'}/$_->{'rand_size'} : 0);
	#print $k++,"\t", $_->{'rand_rate'}*1000,"\n";
    }



    # filter out kmer to use for conservation stat, remove top highly conserved kmers
    my @fooList = keys %motif2HG;
    
    my %foo_gc_motif = (); my $foo_gc='';
    foreach (@fooList) {
	$foo_gc = ($_=~tr/[GC]//);
	push @{$foo_gc_motif{$foo_gc}}, $_;
    }

    my %motif_not_use_list=();
    foreach (keys %foo_gc_motif) {
	@fooList = sort {$motif2HG{$b}<=>$motif2HG{$a}} @{$foo_gc_motif{$_}};
	my $tmpN = int(scalar(@fooList) * 0.25); 
	$tmpN = 0 if (scalar(@fooList) < 2);
	for (my $k=0; $k<$tmpN; $k++) {
	    $motif_not_use_list{$fooList[$k]}=1 if ($motif2HG{$fooList[$k]} > 1);
	}
    }


    # count number of conservation rate based on GC content, for each sequence
    foreach $s (@seqList) {
	$seqA = $s->{'seqA'};
	$seqB = $s->{'seqB'};
	
	$L = length($seqA)-$kmer_size;
	%fooHash=();
	for ($k=0; $k<=$L; $k++) {
	    $tmp = substr($seqA, $k, $kmer_size);
	    next if ($tmp=~/[^ACGT]/);   # skip if anything except ACGT

	    next if (exists $motif_not_use_list{$tmp});   # not use significant motif for stat

	    $fooHash{$tmp}=0;
	}
	foreach $m (keys %fooHash) {
	    $fooHash{$m}=1 if ($seqB=~/$m/);  # conserved, set to 1
	}
	
	@tmpA =(); $#tmpA=$kmer_size; foreach (@tmpA) {$_=0;}
	@tmpB =(); $#tmpB=$kmer_size; foreach (@tmpB) {$_=0;}
	foreach $m (keys %fooHash) {
	    $gc = ($m=~tr/[GC]//);  # number of GC
	    if ($fooHash{$m}==1) {
		$tmpA[$gc]++;
	    } else {
		$tmpB[$gc]++;
	    }
	}
    
	push @{$s->{'gcMotifCountConserv'}}, @tmpA;
	push @{$s->{'gcMotifCountNonConserv'}}, @tmpB;
    }


    # update random matching rate based on average of similar GC motifs
    for ($k=0; $k<@gcStat; $k++) { # only consider at least one GC nt
	$beta = $gcStat[$k]->{'rand_rate'};
	$logBeta = log(1-$beta);
	
	# set background prob for containing no motif for each sequence
	@listZ = (); $#listZ = $#seqList;
	@count_conserv=(); $#count_conserv=$#seqList;
	@count_non_conserv=(); $#count_non_conserv=$#seqList;
    
	$j=0;
	foreach $s (@seqList) {
	    $numB = $s->{'numB'};
	    $listZ[$j] = exp($logBeta*$numB);  # prob of containing no motif, should be approx: 1-numB*beta/1000;
	    $count_conserv[$j] = $s->{'gcMotifCountConserv'}->[$k];   # count for specific GC 
	    $count_non_conserv[$j] = $s->{'gcMotifCountNonConserv'}->[$k]; # count of non conserved motifs for specifc GC
	    $j++;
	}	
	
	$pc = get_optimal_conserv_rate();

	# estimate overall conservation rate
	foreach $s (@seqList) {
	    $numB = $s->{'numB'};
	    $tmpA = 1-(1-$pc)*exp($logBeta*$numB);  # if seqA contains a motif
	    $tmpB = 1-exp($logBeta*$numB); # if seqA does not contain a motif
	    $s->{'neutral_bg_rate'}->[$k] = $tmpB;
	    $s->{'neutral_bg_rate_with_motif_inA'}->[$k] = $tmpA;
	}
	#print "#", $k, "\t", $pc,"\n\n";
    }
    #print "#done with conservation rate estimation\n";
}

