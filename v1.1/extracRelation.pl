#Programmed by guofeifei 2018.1.11,updated in 2021.6.10
#get drug-target relation and target-pathway/GO relation of different drugs
#usage:perl extracRelation.pl $suffix $randtime $pathdb  
#eg., perl extracRelation.pl test 9 wiki
###Requirement###
#input file: 1. combine_score.tsv  2.pathway data,include goa_human.gaf, 
#OUTPUT:simMatrix-$suffix.tsv#similar matrix between drugs based on pathway fingerprints
########clu-$suffix.pdf#hierarchical cluster of drug based on similar Matrix
use strict;
use List::Util qw(shuffle sum);

my $suffix=shift;#eg., test
my $randtime=shift;#larger than 2
my $pathdb=shift; # please choose one of three from GO,wiki,reactome

#read drug-target relation
open(IN,"combine_score.tsv") or die;
my $head=<IN>;
my $pt_relation;
my $pt_node;
my $pt_tar2num;
my $pt_sim;

while(<IN>){
	chomp(my ($drug,$target,$score,$other)=split/\t/);
	$target=uc($target);
	if($score>=0.5){
		$pt_tar2num->{$drug}->{$target}=$score;
	}

}
close IN;


#extract target-function relation
#read function2gene eg goa 
#Step1:read GOA file and store gene2goterm data to $pt_go	
my $pt_gene2go;
if($pathdb eq "GO"){
	$pt_gene2go=readGO();
}
elsif($pathdb eq "wiki"){
	$pt_gene2go=readWikipathway();
}
elsif($pathdb eq "reactome"){
	$pt_gene2go=readreactome();
}
else{
	print "pathway database error!\n";
}

for my $rand(0..($randtime-1)){

	print "random $rand times";
	foreach my $drug(keys %$pt_tar2num){
		my @target=keys %{$pt_tar2num->{$drug}};
		my @target_rand=shuffle(@target);
		for my $i(0..100){
			$pt_relation->{$target_rand[$i]}->{$drug}=$pt_tar2num->{$drug}->{$target_rand[$i]};
			$pt_node->{'drug'}->{$drug}=1;
			$pt_node->{'target'}->{$target_rand[$i]}=1;
		}
	}
	
	
	close IN;
	foreach my $target(keys %{$pt_node->{'target'}}){
		if(exists $pt_gene2go->{$target}){
			foreach my $term(keys %{$pt_gene2go->{$target}}){
				$pt_relation->{$target}->{$term}=1;
				$pt_node->{'term'}->{$term}=1;
			}
		}
	}
	#print dump $pt_node->{'drug'};
	#<STDIN>;


	printrelation($pt_relation,$pt_node);
	
	foreach my $drug(keys %{$pt_node->{'drug'}}){
		my $drugnum=$pt_node->{'drug'}->{$drug};
		system("python Pathsim.py $drugnum");
		print "python Pathsim.py $drugnum";
		open(IN,"result.txt");
		while(<IN>){
			chomp(my @tmp=split/\t/);
			if(exists $pt_sim->{$drug}->{$tmp[2]}){
				my @sims=@{$pt_sim->{$drug}->{$tmp[2]}};
				push @sims,$tmp[3];
				$pt_sim->{$drug}->{$tmp[2]}=[@sims];
			}
			else{
				my @sims;
				push @sims,$tmp[3];
				$pt_sim->{$drug}->{$tmp[2]}=[@sims];
			}
		}
		close IN;
	}
	unlink("drug.txt","target.txt","result.txt","relation.txt","term.txt");

}


open(OUT,">simMatrix-$suffix.tsv") or die;
print OUT "\t".join("\t",sort keys %{$pt_node->{'drug'}})."\n";
foreach my $drug1(sort keys %{$pt_node->{'drug'}}){
	print OUT "$drug1\t";
	foreach my $drug2(sort keys %{$pt_node->{'drug'}}){
		#print OUT "$pt_sim->{$drug1}->{$drug2}\t";
		my @sims=@{$pt_sim->{$drug1}->{$drug2}};
		print OUT (1-sum(@sims)/($#sims+1))."\t";
	}
	print OUT "\n";
}
close OUT;
print "Rscript hclu.r simMatrix-$suffix.tsv clu-$suffix.pdf\n";
system("Rscript hclu.r simMatrix-$suffix.tsv clu-$suffix.pdf");

sub printrelation($pt_relation,$pt_node){
	my $pt_relation=shift;
	my $pt_node=shift;
	my $pt_nodenum;
	
	my $tnum=1;
	foreach my $type(keys %$pt_node){
		open(OUT,">$type.txt") or die;
		my $i=0;
		foreach my $node(sort keys %{$pt_node->{$type}}){
			$pt_node->{$type}->{$node}=$tnum.'0'.$i;
			$pt_nodenum->{$node}=$tnum.'0'.$i;
			print OUT $tnum.'0'.$i."\t".$node."\n";
			$i++;
		}
		close OUT;
		$tnum++;
	}
	
	open(OUT,">relation.txt") or die;
	foreach my $node1(keys %$pt_relation){
		foreach my $node2(keys %{$pt_relation->{$node1}}){
			#print OUT "$node1\t$node2\t$pt_nodenum->{$node1}\t$pt_nodenum->{$node2}\n";
			print OUT "$pt_nodenum->{$node1}\t$pt_nodenum->{$node2}\t0\n";
		}
	}
	close OUT;
	

}

sub readGO(){
	my $pt_gene2go;
	open(IN,"./PathwayData/goa_human.gaf") or die;
	while(<IN>){
		if($_!~/$\!.*/){
			chomp(my @tmp=split/\t/);
			my $symbol=$tmp[2];
			my $goterm=$tmp[4];
			my $goclass=$tmp[8];
			my $species=$tmp[12];
			if($species=~/9606/){
				$pt_gene2go->{$symbol}->{$goterm}=1;
			}
		}
		
	}
	close IN;
	return $pt_gene2go;
}

sub readWikipathway(){
	my $pt_gene2go;
	open(IN,"./PathwayData/Path2gene_wikipathway.txt") or die;
	while(<IN>){
		chomp(my($path,$des,$genes)=split/\t/);
		#if(exists $pt_imm->{$path}){
			my @genes=split(",",$genes);
			foreach my $gene(@genes){
				$pt_gene2go->{uc($gene)}->{$path}=1;
			}
		#}
	}
	close IN;
	return $pt_gene2go;
}

sub readreactome(){
	my $pt_gene2go;
	open(IN,"./PathwayData/gene2path_reactome.txt") or die;
	while(<IN>){
		chomp(my($id,$path,$des,$genes)=split/\t/);
		#if(exists $pt_imm->{$id}){
			my @genes=split(",",$genes);
			foreach my $gene(@genes){
				$pt_gene2go->{uc($gene)}->{$id}=1;
			}
		#}
	}
	close IN;
	return $pt_gene2go;
}