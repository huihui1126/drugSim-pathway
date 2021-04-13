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




foreach my $drug(keys %$pt_tar2num){
	my @target=keys %{$pt_tar2num->{$drug}};
	my @target_rand=shuffle(@target);
	for my $i(0..$#target_rand){
		$pt_relation->{$target_rand[$i]}->{$drug}=$pt_tar2num->{$drug}->{$target_rand[$i]};
		$pt_node->{'drug'}->{$drug}=1;
		$pt_node->{'target'}->{$target_rand[$i]}=1;
	}
}
#extract target-function relation
#read function2gene eg goa 
#Step1:read GOA file and store gene2goterm data to $pt_go
my $pt_gene2go=readinfla();
#my $pt_gene2go=readKEGG();

close IN;
foreach my $target(keys %{$pt_node->{'target'}}){
	if(exists $pt_gene2go->{$target}){
		foreach my $term(keys %{$pt_gene2go->{$target}}){
			$pt_relation->{$target}->{$term}=1;
			$pt_node->{'term'}->{$term}=1;
		}
	}
}
#foreach my $drug1(keys %{$pt_node->{'drug'}}){
foreach my $drug1("Andrographolide"){
	foreach my $drug2(keys %{$pt_node->{'drug'}}){
		if($drug1 ne $drug2){
			my $pt_target1;
			my $pt_term1;
			my $pt_target2;
			my $pt_term2;
			foreach my $target(keys %{$pt_node->{'target'}}){
				if(exists $pt_relation->{$target}->{$drug1}){
					$pt_target1->{$target}=1;
					foreach my $term(keys %{$pt_relation->{$target}}){
						if(exists $pt_node->{'term'}->{$term}){
							$pt_term1->{$term}->{$target}=1;
						}
					}
					
				}
				if(exists $pt_relation->{$target}->{$drug2}){
					$pt_target2->{$target}=1;
					foreach my $term(keys %{$pt_relation->{$target}}){
						if(exists $pt_node->{'term'}->{$term}){
							$pt_term2->{$term}->{$target}=1;
						}
					}
					
				}
				
			}
			open(OOO,">$drug1.$drug2-metapath.txt") or die;
			open(OOO1,">$drug1.$drug2-metapath_net.txt") or die;
			foreach my $term(keys %$pt_term1){
				if(exists $pt_term2->{$term}){
					print OOO "$term\n";
					foreach my $target1(keys %{$pt_term1->{$term}}){
						print OOO1 "$term\t$target1\n";
						print OOO1 "$target1\t$drug1\n";
						print OOO "\t$target1\t$drug1\n";
					}
					foreach my $target2(keys %{$pt_term2->{$term}}){
						print OOO1 "$term\t$target2\n";
						print OOO1 "$target2\t$drug2\n";
						print OOO "\t$target2\t$drug2\n";
					}
				}
			}
			close OOO;
			close OOO1;
		}
	}
}

sub readinfla(){
	my $pt_gene2go;
	open(IN,"positive-negative of inflammatory response.txt") or die;
	while(<IN>){
		chomp(my ($pathway,$symbol)=split/\t/);
		$pt_gene2go->{$symbol}->{$pathway}=1;
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
