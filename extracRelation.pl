#Programmed by guofeifei 2018.1.11
#get drug-target relation and target-pathway/GO relation of different drugs
use strict;
use Data::Dump qw(dump);
use List::Util qw(shuffle sum);

#read drug-target relation
open(IN,"combine_score.txt") or die;
my $head=<IN>;
my $pt_relation;
my $pt_node;
my $pt_tar2num;
my $pt_sim;

while(<IN>){
	chomp(my ($drug,$target,$score,$other)=split/\t/);
	$target=uc($target);
	if($score>=0.4){
		$pt_tar2num->{$drug}->{$target}=$score;
	}

}
close IN;

for my $rand(0..99){

	print "random $rand times";
	foreach my $drug(keys %$pt_tar2num){
		my @target=keys %{$pt_tar2num->{$drug}};
		my @target_rand=shuffle(@target);
		for my $i(0..199){
			$pt_relation->{$target_rand[$i]}->{$drug}=$pt_tar2num->{$drug}->{$target_rand[$i]};
			$pt_node->{'drug'}->{$drug}=1;
			$pt_node->{'target'}->{$target_rand[$i]}=1;
		}
	}
	#extract target-function relation
	#read function2gene eg goa 
	#Step1:read GOA file and store gene2goterm data to $pt_go
	#my $pt_gene2go=readGO();
	my $pt_gene2go=readKEGG();
	
	close IN;
	foreach my $target(keys %{$pt_node->{'target'}}){
		if(exists $pt_gene2go->{$target}){
			foreach my $term(keys %{$pt_gene2go->{$target}}){
				$pt_relation->{$target}->{$term}=1;
				$pt_node->{'term'}->{$term}=1;
			}
		}
	}


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

}


open(OUT,">simMatrix.txt") or die;
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
print "Rscript hclu.r test.pdf\n";
system("Rscript hclu.r test.pdf");

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

	#read GO terms related to immune response
	open(IN,"immuneterm.txt") or die;
	my $pt_imm;
	while(<IN>){
		chomp(my($id,$name)=split/\t/);
		$pt_imm->{$id}=$name;
	}
	my $pt_gene2go;
	open(IN,"goa_human.gaf") or die;
	while(<IN>){
		if($_!~/$\!.*/){
			chomp(my @tmp=split/\t/);
			my $symbol=$tmp[2];
			my $goterm=$tmp[4];
			my $goclass=$tmp[8];
			my $species=$tmp[12];
			
			if($species=~/9606/ and exists $pt_imm->{$goterm}){
				$pt_gene2go->{$symbol}->{$goterm}=1;
			}
		}
		
	}
	close IN;
	return $pt_gene2go;
}

sub readKEGG(){
	my $pt_gene2go;
	open(IN,"inflammationpathway.txt") or die;
	while(<IN>){
		chomp(my ($pathway,$symbol)=split/\t/);
		$pt_gene2go->{$symbol}->{$pathway}=1;
	}
	close IN;
	return $pt_gene2go;
}