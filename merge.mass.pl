#!/usr/bin/perl
use strict;
use warnings;
use Statistics::KernelEstimation;
use Data::Dumper;
use Storable qw(dclone);
use List::Util qw(sum);

my $output_prefix = 'colon_cancer_desi';

my $input_sample_list = 'sample.list.txt'; 
my $work_dir = work_dir;

my $input_lipid = 'LC_MS_Lipid.txt';
my $input_small_mol = 'LC_MS_molecule.neg.txt';

my $mass_cutoff = 0.02;

######################################################################################################
my ($sample_mass, $mass_sample, $mass) = Parsing_Mass_Table( $input_sample_list, $work_dir );
my $lipid = Parsing_Lipid( $input_lipid );
my $small_mol = Parsing_Small_Molecule( $input_small_mol );

my $output_sample_mass = 'mass_dis_in_samples.txt';
&Print_Mass_Diff_By_Samples( $sample_mass, $output_sample_mass );

my $mass_index_group = Group_Mass ( $mass, $lipid, $small_mol, $mass_cutoff );
my $mass_clustered = Clustering_Mass_by_KDE( $mass_index_group, $lipid, $small_mol, $mass_cutoff );

&Print_Clustered_Mass_By_Sample( $mass_clustered, $mass_sample, $lipid, $small_mol, $output_prefix );


######################################################################################################
sub Print_Clustered_Mass_By_Sample {
    my ($mass_clustered, $mass_sample, $lipid, $small_mol, $output_prefix) = @_;

    # get sample names;
    my ( %samples, @samples);
    for my $m ( keys %{ $mass_sample } ){
	for my $s ( keys %{ $mass_sample->{ $m } } ){
	    $samples{ $s }++;
	}
    }
    @samples = sort{ $a cmp $b } keys %samples;

    my $output_tab = $output_prefix.'.clustered_mass.table.with.anno.txt';
    open( my $fh_out, ">", $output_tab ) or die $!;
    my @header_part_1 = @samples;
    my @header_part_2 = map{ $_.'_mass' } @samples;
    my @header_part_3 = ('anno_lipid', 'anno_small_mol');
    print $fh_out join("\t", ( 'Index', @header_part_1, 'Cluster_type', 'Mass_num',
			       @header_part_2, @header_part_3)), "\n";

    for my $index (
	map { $_->[-1] }
	sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
	map { [split(/_/, $_), $_] } 
	keys %{$mass_clustered} ){

	my %sample_c_mass;
	my %sample_c_mass_flag;
	my @anno_lipid;
	my @anno_small_mol;

#	my $num = $mass_clustered->{ $index }{ 'num' };
	my $tag = $mass_clustered->{ $index }{ 'tag' };
	my @mass = @{ $mass_clustered->{ $index }{ 'mass' } };
	
	for my $m ( @mass ){
	    for my $s ( @samples ){
		if( $mass_sample->{ $m }{ $s } ){
		    $sample_c_mass_flag{ $s }++;
		    push @{ $sample_c_mass{ $s } }, $m; 
		}
	    }
	    
	    if( $lipid->{ $m } ){
		push @anno_lipid, join(";", ($m, 
					     join(",", @{ $lipid->{$m}{'IonFormula'} }),
					     join(",", @{ $lipid->{$m}{'LipidIon'} }) ));
	    }

	    if( $small_mol->{ $m } ){
		push @anno_small_mol, join(";", ($m, 
						 join(",", @{$small_mol->{$m}{'Formula'}}),
						 join(",", $small_mol->{$m}{'Tag'}),
						 join(",", @{$small_mol->{$m}{'Name'}})) );
	    }
	} # for $m
	
	my @output_flag;
	my @output_mass;
	my $output_anno_lipid;
	my $output_anno_small_mol;

	for my $s ( @samples ){
	    if( $sample_c_mass_flag{ $s } ){
		push @output_flag, 1;
	    } else {
		push @output_flag, 0;
	    }

	    if( $sample_c_mass{ $s } ){
		push @output_mass, join(",", @{ $sample_c_mass{ $s } });
	    } else {
		push @output_mass, 'None';
	    }
	}

	if( $anno_lipid[0] ){
	    $output_anno_lipid = join("|", @anno_lipid);
	}
	else {
	    $output_anno_lipid = 'None';
	}

	if( $anno_small_mol[0] ){
	    $output_anno_small_mol = join("|", @anno_small_mol);
	} else {
	    $output_anno_small_mol = 'None';
	}

	my $num = sum( @output_flag );

	print $fh_out join("\t", $index,
			   @output_flag,
			   $tag, $num,
			   @output_mass,
			   $output_anno_lipid,
			   $output_anno_small_mol), "\n";
    }
    return 0;
}






sub Clustering_Mass_by_KDE {
    my $mass_index_group = shift;
    my $lipid            = shift;
    my $small_mol        = shift;
    my $mass_cutoff      = shift;

    my %clustered_mass;

    my $weight_lipid     = 1;
    my $weight_small_mol = 1;
    my $weight_mass      = 1;

    for my $index ( sort {$a <=> $b} keys %{ $mass_index_group } ){
	my @mass = sort { $a <=> $b } @{ $mass_index_group->{$index} };

	my $mass_diff = $mass[-1] - $mass[0];
	my $mass_num  = scalar @mass;
	
	if( $mass_num == 1 ){
	    my $tag;
	    if( $lipid->{ $mass[0] } ){
		$tag = 'Solo_Lipid';
	    }
	    elsif( $small_mol->{ $mass[0] } ){
		$tag = 'Solo_Small_Mol';
	    }
	    else{
		$tag = 'Solo_Mass';
	    }
	    $clustered_mass{ $index }{ 'mass' } = [@mass];
	    $clustered_mass{ $index }{ 'num'  } = $mass_num;
	    $clustered_mass{ $index }{ 'tag'  } = $tag;
	} # good
	elsif ( $mass_num > 1 && $mass_num <=3 && $mass_diff <= $mass_cutoff ){
	    $clustered_mass{ $index }{ 'mass' } = [@mass];
	    $clustered_mass{ $index }{ 'num'  } = $mass_num;
	    $clustered_mass{ $index }{ 'tag'  } = 'Cluster_No_KDE';
	} #good
	elsif ( $mass_num > 1 && $mass_num <= 3 && $mass_diff > $mass_cutoff ){

	    my @mass_split;
	    my $index_split = 0;
	    for (my $i = 0; $i <= ($mass_num-1); $i++){
		if( $i == 0 ){
		    push @mass_split, $mass[$i];
		}
		elsif( $mass[$i]-$mass[$i-1] <= $mass_cutoff ){
		    push @mass_split, $mass[$i];
		}
		elsif( $mass[$i]-$mass[$i-1] > $mass_cutoff ){
		    $clustered_mass{ $index.'_'.$index_split }{ 'mass' } = [@mass_split];
		    $clustered_mass{ $index.'_'.$index_split }{ 'num'  } = scalar @mass_split;
		    $clustered_mass{ $index.'_'.$index_split }{ 'tag'  } = 'Cluster_No_KDE';
		    
		    @mass_split='';
		    push @mass_split, $mass[$i];
		    $index_split++;
		} else {
		    print STDERR 'unknow case:', "\n";
		    print STDERR Dumper( @mass );
		    exit 1;
		}
	    }
	    $clustered_mass{ $index.'_'.$index_split }{ 'mass' } = [@mass_split];
	    $clustered_mass{ $index.'_'.$index_split }{ 'num'  } = scalar @mass_split;
	    $clustered_mass{ $index.'_'.$index_split }{ 'tag'  } = 'Cluster_No_KDE';
	}
	elsif ( $mass_num > 3 && $mass_diff <= $mass_cutoff ){
	    $clustered_mass{ $index }{ 'mass' } = [@mass];
	    $clustered_mass{ $index }{ 'num'  } = $mass_num;
	    $clustered_mass{ $index }{ 'tag'  } = 'Cluster_No_KDE';
	} #good
	elsif ( $mass_num > 3 && $mass_diff > $mass_cutoff ){
	    # return probability density
	    my $pd = Kernel_Density_Estimate ( \@mass, $weight_lipid, $weight_small_mol, $weight_mass );
	    my $pd_bak = dclone( $pd ); # shift of array will empty the array

	    # find probability density clusters
	    my $clusters  = Find_Probability_Density_Regions( $pd_bak );
	    # print table for checking the correction of group the mass
	    &Print_KDE_and_Cluster( $pd, $clusters, $index, \@mass );


	    my @mass_grouped;
	    for my $cluster ( @{ $clusters } ){
		
	    my $mass_start = $cluster->[0];
	    my $mass_end   = $cluster->[1];
	    my @mm;
	    for my $m ( @mass ){
		if( $m >= $mass_start && $m <= $mass_end ){
		    push @mm, $m;
		}
	    }
	    push @mass_grouped, [@mm] if $mm[0];
	    }

	    for ( my $i=0; $i <= $#mass_grouped; $i++ ){
		my $m_grouped = $mass_grouped[$i];
		my $m_diff = $m_grouped->[-1] - $m_grouped->[0];
		my $m_num  = scalar @{ $m_grouped };
		
		my $clustering_tag;

		if( $m_num == 1 ){
		    if( $lipid->{ $m_grouped->[0] } ){
			$clustering_tag = 'KDE_Solo_Lipid';
		    }
		    elsif( $small_mol->{ $m_grouped->[0] } ){
			$clustering_tag = 'KDE_Solo_Small_Mol';
		    }
		    else{
			$clustering_tag = 'KDE_Solo_Mass';
		    }
		}
		elsif( $m_num > 1 && $m_diff <= $mass_cutoff ){
		    $clustering_tag = 'KDE_No_Exceed';
		}
		elsif( $m_num > 1 && $m_diff > $mass_cutoff ){
		    $clustering_tag = 'KDE_Exceed';
		}
		else {
		    print STDERR 'unkown case: KDE clustering!', "\n";
		    print STDERR join("\t", ('mass_num:', $m_num, 'mass_diff', $m_diff)), "\n";
		    exit 1;
		}
		$clustered_mass{ $index.'_'.$i }{ 'mass' } = $m_grouped;
		$clustered_mass{ $index.'_'.$i }{ 'num'  } = $m_num;
		$clustered_mass{ $index.'_'.$i }{ 'tag'  } = $clustering_tag;
	    } # for
	} # elsif ( $mass_num > 1 && $mass_diff > $mass_cutoff )
	else {
	    print STDERR 'unkown case: grouping!', "\n";
	    print STDERR join("\t", ('mass_num:', $mass_num, 'mass_diff', $mass_diff)), "\n";
	    exit 1;
	}
    }
    
    return \%clustered_mass;
}

sub Find_Probability_Density_Regions {
    my $pd = shift;

    # more than 20% of summit ( $pd_max )
    my $peak_change_cutoff = 0.2;

    # fixed bug here and clean the codes
    # up tandency
    if( $pd->[1]->[1] >= $pd->[0]->[1] ){
        my @regions;

        my $pd_min_left  = shift @{ $pd };
        my $pd_min_right = $pd_min_left;
	my $pd_max       = $pd_min_left;
        my $pd_min_right_val;

        while( my $e = shift @{ $pd } ){
            if( $pd->[0] ){
                if( ! $pd_min_left && ! $pd_min_right && ! $pd_max ){
                    $pd_min_left  = $e;
                    $pd_min_right = $e;
                    $pd_max       = $e;
                } # check
                elsif( $e->[1] >= $pd_min_right->[1] && $e->[1] <= $pd->[0]->[1] ){
                    $pd_min_right = $e;
                    $pd_max       = $e;
                } # check
                elsif( $e->[1] >= $pd_min_right->[1] && $e->[1] > $pd->[0]->[1] ){
                    $pd_min_right = $e;
                    $pd_max       = $e;
                } # check
                elsif( $e->[1] <= $pd_min_right->[1] && $e->[1] >= $pd->[0]->[1] ){ # check
                    $pd_min_right = $e;
                }
                elsif ( $e->[1] <= $pd_min_right->[1] && $e->[1] <= $pd->[0]->[1] ){ # check
                    $pd_min_right = $e;
                    push @regions, [$pd_min_left, $pd_min_right, $pd_max];

                    undef $pd_min_left;
                    undef $pd_min_right;
                    undef $pd_max;
                }
                else {
                    print STDERR "unkonwn case: up tandency\n";
                    print STDERR join("\t", (@{$pd_min_left}, ';', @{$pd_max}, ';', @{$pd_min_right})), "\n";
                    print STDERR Dumper( $pd );
                    exit 1;
                }
            }
            else {
                # fixed the bug here
                if( $pd_min_left->[0] && $pd_max->[0] ){
                    $pd_min_right = $e;
                    push @regions, [$pd_min_left, $pd_min_right, $pd_max];
                } elsif( !$pd_min_left->[0] && !$pd_max->[0] ){
                    $pd_min_left  = $e;
                    $pd_min_right = $e;
                    $pd_max       = $e;
		    push @regions, [$pd_min_left, $pd_min_right, $pd_max];
                } else {
                    print STDERR 'up tandency: unknow case',"\n";
                    exit 1;
                }
            }
        }

        my $region_merged = Merge_Regios( \@regions, $peak_change_cutoff );
        return $region_merged;

    } # end if; end up tandency




    
 # down tandency
    elsif( $pd->[1]->[1] < $pd->[0]->[1] ){
	my @regions;

	my $pd_min_left  = shift @{ $pd };
	my $pd_min_right = $pd_min_left;
	my $pd_max       = $pd_min_left;

	while( my $e = shift @{ $pd } ){
	    if( $pd->[0] ){
		if( ! $pd_min_left && ! $pd_min_right && ! $pd_max ){
		    $pd_min_left  = $e;
		    $pd_min_right = $e;
		    $pd_max       = $e;
		}
		elsif( $e->[1] <= $pd_min_right->[1] && $e->[1] >= $pd->[0]->[1] ){
		    $pd_min_right = $e;
		}
		elsif ( $e->[1] <= $pd_min_right->[1] && $e->[1] < $pd->[0]->[1] ){
		    $pd_min_right= $e;
		    push @regions, [$pd_min_left, $pd_min_right, $pd_max];

		    undef $pd_min_left;
		    undef $pd_min_right;
		    undef $pd_max;

		}
		elsif ( $e->[1] >= $pd_min_right->[1] && $e->[1] <= $pd->[0]->[1] ){
		    $pd_min_right = $e;
		    $pd_max       = $e;
		}
		elsif ( $e->[1] >= $pd_min_right->[1] && $e->[1] >= $pd->[0]->[1] ){
		    $pd_min_right = $e;
		    $pd_max       = $e;
		}

#		elsif ( $e->[1] >= $pd_max->[1] && $e->[1] <= $pd->[0]->[1] ){
#		    $pd_min_right = $e;
#		    $pd_max       = $e;
#		}
#		elsif ( $e->[1] >= $pd_max->[1] && $e->[1] > $pd->[0]->[1] ){
#		    $pd_min_right = $e;
#		    $pd_max       = $e;
#		}

		else {
		    print STDERR 'unkown case: down tandency', "\n";
		    print STDERR join("\t", (@{$pd_min_left}, ';', @{$pd_max}, ';', @{$pd_min_right})), "\n";
		    print STDERR Dumper( $pd );
		    exit 1;
		}
	    } else {
		# fixed the bug here
		if( $pd_min_left->[0] && $pd_max->[0] ){
		    $pd_min_right = $e;
		    push @regions, [$pd_min_left, $pd_min_right, $pd_max];
		} elsif( !$pd_min_left->[0] && !$pd_max->[0] ){
		    $pd_min_left  = $e;
		    $pd_min_right = $e;
		    $pd_max       = $e;
		    push @regions, [$pd_min_left, $pd_min_right, $pd_max];
		} else {
		    print STDERR 'down tandency: unknow case', "\n";
		    exit 1;
		}
	    } # while if and else end
	}

	my $region_merged = Merge_Regios( \@regions, $peak_change_cutoff );
	return $region_merged;
    } # end elsif; end down tandency 

# other case
    else {
	print STDERR "unconsidered case:\n";
	print STDERR 'NOT: $pd->[1]->[1] > $pd->[0]->[1]', "\n";
	print STDERR 'print the 1st and 2nd intensity values:', "\n";
	print STDERR join("\t", ($pd->[1]->[1], $pd->[0]->[1])), "\n";
	exit 1;
    }
}


sub Merge_Regios {
    my $regions            = shift;
    my $peak_change_cutoff = shift;

    my @region_merged;
    my $peak_1st           = shift @{$regions};


    my $summit_density     = $peak_1st->[2]->[1];
    my $right_side_density = $peak_1st->[1]->[1];
    
    my $region_merge_left  = $peak_1st->[0]->[0];
    my $region_merge_right = $peak_1st->[1]->[0];
    
    while( my $r = shift @{$regions} ){
	if( ( $summit_density - $right_side_density ) <= $summit_density * $peak_change_cutoff ){
	    
	    $region_merge_right = $r->[1]->[0];
	    $summit_density     = $r->[2]->[1];
	    $right_side_density = $r->[1]->[1];
	}
	else {
	    push @region_merged, [$region_merge_left, $region_merge_right];
	    
	    undef $region_merge_left;
	    undef $region_merge_right;
	    undef $summit_density;
	    undef $right_side_density;
	    
	    $region_merge_left  = $r->[0]->[0];
	    $region_merge_right = $r->[1]->[0];
	    $summit_density     = $r->[2]->[1];
	    $right_side_density = $r->[1]->[1];
	}
	
    }
    push @region_merged, [$region_merge_left, $region_merge_right];

    return \@region_merged;
}



sub Kernel_Density_Estimate {

    my ($mass, $weight_lipid, $weight_small_mol, $weight_mass) = @_;

    my @pdf;
    my $bin_num = 200;

    my $s = Statistics::KernelEstimation->new_epanechnikov();    
#    my $s = Statistics::KernelEstimation->new_gauss();

    for my $m ( @{ $mass } ){
	if( $lipid->{ $m } ){
	    $s->add_data( $m, $weight_lipid );
	}
	elsif( $small_mol->{ $m } ){
	    $s->add_data( $m, $weight_small_mol );
	}
	else {
	    $s->add_data( $m, $weight_mass );
	}
    }

# not work for new_epanechnikov
# $s->optimal_bandwidth();
# Probability density functions

    my $w = $s->default_bandwidth();
    my ( $min, $max ) = $s->extended_range();
    for( my $x=$min; $x<=$max; $x+=($max-$min)/$bin_num ) {
#	print $x, "\t", $s->pdf( $x, $w ), "\t", $s->cdf( $x, $w ), "\n";
	push @pdf, [$x, $s->pdf( $x, $w )];
    }

    return \@pdf;
}


sub Group_Mass {
    my ( $mass, $lipid, $small_mol, $mass_cutoff ) = @_;

    my @mass_all;
    push @mass_all, keys %{ $mass };
    push @mass_all, keys %{ $lipid };
    push @mass_all, keys %{ $small_mol };
    @mass_all = sort{ $a <=> $b } @mass_all;

    my %mass_index_group;
    my @mass_group;
    my $index;
    my $mass_1st;

    for ( my $i = 0; $i <= $#mass_all; $i++ ){
	if( ! $mass_1st && ! @mass_group ){
	    $mass_1st = $mass_all[$i];
	    push @mass_group, $mass_all[$i];
	} else {
	    if( ($mass_all[$i] - $mass_all[$i-1]) <= $mass_cutoff ){
		push @mass_group, $mass_all[$i];
	    } elsif( ($mass_all[$i] - $mass_all[$i-1]) > $mass_cutoff ){
		$index++;
		$mass_index_group{ $index } = [@mass_group];

		undef $mass_1st;
		undef @mass_group;
		
		$mass_1st = $mass_all[$i];
		push @mass_group, $mass_all[$i];
	    } else {
		print STDERR 'unkonwn case: ', join("\t", ($mass_all[$i-1], $mass_all[$i], $mass_all[$i+1])), "\n";
		exit 1;
	    }
	}
    }
    
    $index++;
    $mass_index_group{ $index } = [@mass_group];
    undef $mass_1st;
    undef @mass_group;

    return \%mass_index_group;
}

sub Print_Mass_Diff_By_Samples {
    my $sample_mass = shift;
    my $output      = shift;

    open( my $fh_out, ">", $output ) or die $!;
    print $fh_out join("\t", ('mass_dis', 'sample')), "\n";
    
    for my $sample ( sort { $a cmp $b } keys %{ $sample_mass } ){
	my @mass = sort { $a <=> $b } keys %{ $sample_mass->{ $sample } };
	for ( my $i = 0; $i < $#mass; $i++ ){
	    my $dis = $mass[$i+1] - $mass[$i];
	    print $fh_out join("\t", ($dis, $sample)), "\n";
	}
    }
    close $fh_out;
    return 0;
}

sub Parsing_Small_Molecule {
    my $input = shift;
    my %small_mol;

    my $mass_H    = 1.0078;
    my $mass_HCOO = 44.9977;
    
    open( my $fh_in, $input ) or die $!;
    while( my $line = <$fh_in> ){
	next if $line =~ /^Name/;
	chomp $line;
	my ($name, $formula, $mass) = split(/\t/, $line);

	# fix a bug here
	# add H
	push @{ $small_mol{ $mass - $mass_H }{ 'Name' } }, $name;
	push @{ $small_mol{ $mass - $mass_H }{ 'Formula' } }, $formula;
	$small_mol{ $mass - $mass_H }{ 'Tag' } = 'H';
#	push @{ $small_mol{ $mass + $mass_H }{ 'Name' } }, $name;
#	push @{ $small_mol{ $mass + $mass_H }{ 'Formula' } }, $formula;
#	$small_mol{ $mass + $mass_H }{ 'Tag' } = 'H';

	# add HCOO
	push @{ $small_mol{ $mass + $mass_HCOO }{ 'Name' } }, $name;
	push @{ $small_mol{ $mass + $mass_HCOO }{ 'Formula' } }, $formula;
	$small_mol{ $mass + $mass_HCOO }{ 'Tag' } = 'HCOO';
    }
    close $fh_in;
    return \%small_mol;
}


sub Parsing_Lipid {
    my $input = shift;
    my %lipid;

    open( my $fh_in, $input ) or die $!;
    while( my $line = <$fh_in> ){
	next if $line =~ /^LipidIon/;
	chomp $line;
	my ($lipidIon, $mass, $formula) = split(/\t/, $line);
	push @{ $lipid{ $mass }{ 'LipidIon' } }, $lipidIon;
	push @{ $lipid{ $mass }{ 'IonFormula' } }, $formula;
    }
    close $fh_in;
    return \%lipid;
}

sub  Parsing_Mass_Table {
    my $input_list = shift;
    my $dir_sample = shift;

    my (%sample_mass, %mass_sample, %mass);
    
    open( my $fh_in_1, $input_list ) or die $!;
    while( my $sample_name = <$fh_in_1> ){
	chomp $sample_name;
	my $input_mass = $dir_sample.'/'.$sample_name.'/'.$sample_name.'_3k_signal.lock_mass.txt';
	open( my $fh_in_2, "head -n 1 $input_mass | " ) or die $!;
	while( my $line = <$fh_in_2> ){
	    chomp $line;
	    my @list = split(/\t/, $line);

	    for ( my $i = 3; $i <= $#list; $i++ ){
		$sample_mass{ $sample_name }{ $list[$i] }++;
		$mass_sample{ $list[$i] }{  $sample_name }++;
		$mass{ $list[$i] }++;
	    }
	}
	close $fh_in_2;
    }
    close $fh_in_1;
    return ( \%sample_mass, \%mass_sample, \%mass );
}


sub Print_KDE_and_Cluster {
    my ( $pd, $clusters, $index, $mass ) = @_;

    open( my $fh_out, ">", 'plot.kde.index_'.$index.'.txt' ) or die $!;
    my @kde_mass;
    for my $mmm ( @{ $pd } ){
	push @kde_mass, $mmm->[1];
	print $fh_out join("\t", ($mmm->[0], $mmm->[1], 'KDE')), "\n";
    }
    @kde_mass = sort { $a <=> $b } @kde_mass;
    for my $ccc ( @{ $clusters } ){
	print $fh_out join("\t", ($ccc->[0], $kde_mass[0], 'Cluster')), "\n";
	print $fh_out join("\t", ($ccc->[1], $kde_mass[0], 'Cluster')), "\n";
    }

    for my $m ( @{ $mass } ){
	print $fh_out join("\t", ($m, $kde_mass[0], 'Mass')), "\n";
    }
    close $fh_out;
    return 0;
}
