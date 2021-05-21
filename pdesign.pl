#!/usr/bin/perl
my $former=shift @ARGV;
my $species=shift @ARGV;
my $label=shift @ARGV;
my $thread=shift @ARGV;
########1
open IN,"<spool.txt";
open OUT1,">$species.fasta";
open OUT2,">non_$species.fasta";
my $id;
my $bool;
while(<IN>){
        chomp;
        if(/\>(.*)/){
                $id=$1;
                if(/$former/){
                        if(/$species/){
                                print OUT1 "\n";
                                my @info=split /\s/,$id;
                                print OUT1 ">$info[0]\n";
                                $bool=1;
                        }else{
                                $bool=2;
                                print OUT2 "\n";
                                my @info=split /\s/,$id;
                                print OUT2 ">$info[0]\n";
                        }
                }else{
                        $bool=0;
                }
        }else{
                if($bool==1){
                        print OUT1 "$_";
                }elsif($bool==2){
                        print OUT2 "$_";
                }
        }
}
close IN;
close OUT1;
close OUT2;
`sed -i '1d' $species\.fasta`;
`sed -i '1d' non_$species\.fasta`;

#######2
my $kmin=22;
my $kmax=25;
my $cut=1;
my $file=$species.".fasta";
my $out=$label;
my @kmers;
my $nu=$kmin;
while($nu <= $kmax){
        push @kmers,$nu;
        $nu+=1;
}
foreach my $kmer(@kmers){
	`jellyfish count -m $kmer -s 102410240 -t $thread -o $out\_$kmer -L $cut $file`;
	`jellyfish dump -L $cut $out\_$kmer\_0 > $out\_$kmer.fa`;
}
`cat $out\_*.fa >$out.fa`;
`rm $out\_*.fa`;
`rm $out\_*\_0`;

$file="non_".$species.".fasta";
$out="non_".$label;
foreach my $kmer(@kmers){
        `jellyfish count -m $kmer -s 102410240 -t $thread -o $out\_$kmer -L $cut $file`;
        `jellyfish dump -L $cut $out\_$kmer\_0 > $out\_$kmer.fa`;
}
`cat $out\_*.fa >$out.fa`;
`rm $out\_*.fa`;
`rm $out\_*\_0`;

########3
$file=$label;
my %nc;
open IN,"<non_$file.fa";
while(<IN>){
        chomp;
        if(/\>/){
        }else{
                $nc{$_}=1;
        }
}
close IN;

$nu=0;
open IN,"<$file.fa";
open OUT,">$file\_specific.fa";
while(<IN>){
        chomp;
        if(/\>/){
        }else{

        if($nc{$_}){
        }else{
                print OUT ">$nu\n$_\n";
                $nu+=1;
        }
        }
}
close IN;
close OUT;

##############4
`cd-hit -i $label\.fa -o $label\_clust.fa -c 0.88 -aS 0.8 -d 0`;

##########5
$nu=0;
$file=$label."_clust.fa";
open IN,"<$file";
open OUT,">$label\_clust_used.fa";
while(<IN>){
        chomp;
        if(/\>/){
                print OUT "\>seq$nu\n";
                $nu+=1;

        }else{
                print OUT "$_\n";
        }
}
close IN;
close OUT;

##########6
`patman -P $label\_clust_used.fa -D $species\.fasta -o $label\_cover.out`;
`patman -P $label\_clust_used.fa -D non_$species\.fasta -o non_$label\_cover.out`;

#########7
my $t1=(split /\s+/,`grep -c '>' $species\.fasta`)[0];
my $t2=(split /\s+/,`grep -c '>' non_$species\.fasta`)[0];
my %c1;
my %c;
open IN,"<$label\_cover.out";
while(<IN>){
        chomp;
        my @info=split /\t/;
        $c1{$info[1]}+=1/$t1;
        $c{$info[1]}+=1/$t1;
}
close IN;

my %c2;
open IN,"<non_$label\_cover.out";
while(<IN>){
        chomp;
        my @info=split /\t/;
        $c{$info[1]}-=1/$t2;
        $c2{$info[1]}-=1/$t2;
}
close IN;

my %seq;
my $seq;
open IN,"<$label\_clust_used.fa";
while(<IN>){	
	chomp;
	if(/\>(.*)/){
		$seq=$1;
	}else{
		$seq{$seq}=$_;
	}
}
close IN;

open OUT,">seq.txt";
foreach my $key(sort {$c{$b}<=>$c{$a}} keys %c){
	print OUT "$seq{$key}\n";
}
close OUT;
`RNAfold --jobs=$thread --infile=seq.txt --outfile=structure-log.txt`;

my $seqs;
my %ss;
open IN,"<structure-log.txt";
while(<IN>){
	if(/[AUCG]+/){
		$seqs=$_;
	}elsif(/\((.*)\)/){
		$ss{$seqs}=$1;
	}
}
close IN;	

print "key\tseq\tscore\ttarget\tnon_target\tfree_energy\n";
foreach my $key(sort {$c{$b}<=>$c{$a}} keys %c){
        print "$key\t$seq{$key}\t$c{$key}\t$c1{$key}\t$c2{$key}\t$ss{$seq{$key}}\n";
}


