#! /usr/bin/perl

#################################
# HIV-DRLink_github.pl
# 
# Wei Shao
# 10-16-2019
#
# This script parses JSON files produced by a python client for Stanford HIVDB
# Usage: HIV-DRLink_github.pl <input json file>
# Note, The input sequences with only "other" mutations are excluded from further analysis. Therefore, the number of left sequences may be smaller than the input sequences and the input JSON files. 
# For RT, only major DRM are included. 
###############################

use JSON; # imports encode_json, decode_json, to_json and from_json.
use strict;
use warnings;
use Data::Dumper;

my $filename = shift;

my ($outfile, $mutation_list, $longest);
my (@PI_mutations, @all_PRmutations, %temp_PR,  $all_PR_mut_list, $PR_number);
my (@RT_mutations, @all_RTmutations, %temp_RT,  $all_RT_mut_list, $RT_number, @RT_others, %RT_DRM, %RT_major);
my (@IN_mutations, @all_INmutations, %temp_IN,  $all_IN_mut_list, $IN_number);
my ($input_seq_count, $all_PR_RT_IN_mutations, @all_mutations,  %ID_PR_RT_IN);
my ($ID, $json_text, $json, $data, $element, $number, $json_title, $i, $n, $m);
my ($id, %ID_seq, @ID_line, $mutation_file);
my ($seq, %unique_seq,  %unique_seq2, %essential, %newhash, $cu_base, $u_base);
my ($cu_del, $u_del, $diff, $seq_length, @id_list, $id_list_size, $seq_with_DRM);
my ($found, $mutation_type, @found_mut, %mut_count);
my ($first,  $mut, @DRM_pattern, $DRM_list, $pattern_sum, $var, $each_pattern_percent);
my ($enzyme_line, $total_pattern_percent, $DRM_percent);
#########################################################

$outfile = $filename."_Mutations.fas";
 
%RT_major = (
    M41L => 1,
    K65R => 1,
    D67N => 1,
    T69Insertion => 1,
    K70E => 1,
    K70R => 1,
    L74V => 1,
    L74I => 1,
    L100I => 1,
    K101E => 1,
    K101P => 1,
    K103N => 1,
    K103S => 1,
    V106A => 1,
    V106M => 1,
    Y115F => 1,
    Q151M => 1,
    Y181C => 1,
    Y181I => 1,
    Y181V => 1,
    M184V => 1,
    M184I => 1,
    Y188L => 1,
    Y188C => 1,
    Y188H => 1,
    G190A => 1,
    G190S => 1,
    G190E => 1,
    G190Q => 1,
    L210W => 1,
    T215F => 1,
    T215Y => 1,
    K219Q => 1,
    K219E => 1,
    M230L => 1
    );

$input_seq_count = 0;   
@all_PRmutations = ();
%temp_PR = ();
@all_RTmutations = ();
%temp_RT = ();
@all_INmutations = ();
%temp_IN = ();


$json_text = do {
   open(my $json_fh, "<:encoding(UTF-8)", $filename)
      or die("Can't open \$filename\": $!\n");
   local $/;
   <$json_fh>
};

$json = JSON->new;
$data = $json->decode($json_text);
my $count = 0;
for $element (@$data) {   ## $element is a hash reference. 
    $mutation_list = $element->{'mutations'};
    $number = scalar@$mutation_list;

    foreach $json_title(sort keys %{$element} ) {  
	
	@PI_mutations = ();
	@RT_mutations = ();
	@RT_others = ();
        @IN_mutations = ();
       
	if ($json_title =~ /inputSequence/i) { 
	    $count++;
	    $ID = $$element{$json_title}->{'header'};
            $ID =~ s/\n//g;
            $input_seq_count++;            
        }
        else {  

	   for ($i=0; $i<$number; $i++  ) {

               if ( $$element{$json_title}[$i]->{'gene'}->{'name'} eq "PR") {	  
		  
		   if ($$element{$json_title}[$i]->{'primaryType'}  =~ /Major|Accessory/i) {    #### no "Other"
                      push @PI_mutations, $$element{$json_title}[$i]->{'text'};   #### PI only
		     $temp_PR{$$element{$json_title}[$i]->{'text'} } =1;   
                   }


                }    
	       elsif ( $$element{$json_title}[$i]->{'gene'}->{'name'} eq "RT") {	 

		   if ($$element{$json_title}[$i]->{'primaryType'}  =~ /NNRTI/i) {                

		       if ($RT_major{$$element{$json_title}[$i]->{'text'} } ) {     ### only take RT major 
	
			   push @RT_mutations, $$element{$json_title}[$i]->{'text'};  ## NNTRI first
			   $temp_RT{$$element{$json_title}[$i]->{'text'} } =1;
		       }
		   }

		   if ($$element{$json_title}[$i]->{'primaryType'} =~ /\bNRTI\b/i) {
                  

		       if ($RT_major{$$element{$json_title}[$i]->{'text'} } ) {     ### only take RT major 
		
			   push @RT_mutations, $$element{$json_title}[$i]->{'text'};  ## then NRTI 
			   $temp_RT{$$element{$json_title}[$i]->{'text'} } =1;
		       }
		   }

	      }
              elsif ( $$element{$json_title}[$i]->{'gene'}->{'name'} eq "IN") {	 
                  if ($$element{$json_title}[$i]->{'primaryType'}  =~ /Major/i) {
     
			
			   push @IN_mutations, $$element{$json_title}[$i]->{'text'};  ## NNTRI first
			   $temp_IN{$$element{$json_title}[$i]->{'text'} } =1;
		      
		   }
                   if ($$element{$json_title}[$i]->{'primaryType'}  =~ /Accessory/i) {
             	
			   push @IN_mutations, $$element{$json_title}[$i]->{'text'};  ## NNTRI first
			   $temp_IN{$$element{$json_title}[$i]->{'text'} } =1;
		      
		   }

	      }	     
	  }    
	   @all_mutations = (@PI_mutations, @RT_mutations,@IN_mutations);
	  
	   $ID_PR_RT_IN{$ID} = "@all_mutations";
	  
	}                
    }          
}  

   $longest = 0;
   @all_PRmutations = keys %temp_PR;
   $all_PR_mut_list = "@all_PRmutations";
   $all_PR_mut_list = ORGANIZEMUTS($all_PR_mut_list, \$longest);
   $PR_number = scalar@all_PRmutations;

   @all_RTmutations = keys %temp_RT;
   $all_RT_mut_list = "@all_RTmutations";
   $all_RT_mut_list = ORGANIZEMUTS($all_RT_mut_list, \$longest);  
   $RT_number = scalar@all_RTmutations;

   @all_INmutations = keys %temp_IN;
   $all_IN_mut_list = "@all_INmutations";
$all_IN_mut_list = ORGANIZEMUTS($all_IN_mut_list, \$longest);
   $IN_number = scalar@all_INmutations;

    $all_PR_RT_IN_mutations = " ".$all_PR_mut_list." ".$all_RT_mut_list." ".$all_IN_mut_list;

 
   open (OUT, ">$outfile");

   print OUT "enzyme#\t";
if ( ($PR_number > 0)&&($RT_number >0) && ($IN_number >0) ){    ### all 3 enzymes
    print OUT "PR DRM\t", "\t"x($PR_number-1), "RT DRM", "\t"x$RT_number, "IN DRM\n";   
}
if (($PR_number > 0)&&($RT_number >0) && ($IN_number ==0) ){   ### PR and RT only

    print OUT "PR DRM\t", "\t"x($PR_number-1), "RT DRM\n",
}
if (($PR_number == 0)&&($RT_number >0) && ($IN_number >0) ){   ## RT and IN only
    print OUT "RT DRM", "\t"x$RT_number, "IN DRM\n";
}

if (($PR_number > 0)&&($RT_number == 0) && ($IN_number >0) ){   ## PR and IN only  
    print OUT "PR DRM", "\t"x$PR_number, "IN DRM\n";
}

if (($PR_number > 0)&&($RT_number ==0) && ($IN_number ==0) ){   ## PR only
     print OUT "PR DRM\n";
}
if (($PR_number == 0)&&($RT_number >0) && ($IN_number ==0) ){   ## RT only
     print OUT "RT DRM\n";
}
if (($PR_number == 0)&&($RT_number ==0) && ($IN_number >0) ){   ## IN only
     print OUT "IN DRM\n";
}

FASTA_OUT ($all_PR_RT_IN_mutations, \%ID_PR_RT_IN, $filename, $input_seq_count);  
#################################################################################
$mutation_file = $outfile;

$outfile = $filename."_DRM_pattern_percent_08192019.xls"; 
$seq_with_DRM = 0;

open (F, "$mutation_file") || die "can not open $mutation_file\n";
open (OUT, ">$outfile") || die "cann't open $outfile\n";


while (<F>) {

	chomp;
        next unless ($_ =~ /\w+/);
        next if ($_ =~ /total/);
	if ($_ =~ /enzyme/) {
	  ($first,  $enzyme_line) = split(/#/, $_);
	  
        } 
	elsif ($_ =~ /mutation_types/) {
	    $mutation_type = $_;
        }
	elsif ($_ =~ />/)  {

	    $ID = $_;	    
        }
	else {
	    $ID_seq{$ID} = $_;
	    
        }
}

($first, @found_mut) = split (/\t/, $mutation_type);

print OUT "\t\t$enzyme_line\n";   


print OUT "DRM pattern\t# sequence with same pattern\tDRM pattern% "; 

foreach $mut (@found_mut) {
    print OUT "$mut\t";
    $mut_count{$mut} = 0;
    
}

print OUT "\n";
	
%newhash = reverse %ID_seq;   #get ride of identical sequences. 
%unique_seq = reverse %newhash;
%unique_seq2 = reverse %newhash;

 $n = 0;
 $m = 0;
$found = 0;

foreach $ID (sort keys %unique_seq) {
 
    next unless ($unique_seq{$ID});  
  


    $seq_length = length($unique_seq{$ID});

    
    foreach $id (sort keys %unique_seq2) {

        next if ($ID eq $id);
	next if (length$unique_seq2{$id} != $seq_length); 
    
	$cu_del = 0;
	$u_del = 0;
	$diff = 0;

	for ($i = 0; $i <= ($seq_length -1); $i++) {
	    $cu_base = substr($unique_seq{$ID}, $i, 1);
	    $u_base = substr($unique_seq2{$id}, $i, 1);


	    if ( ($cu_base eq "-") || ($u_base  eq "-") || ($cu_base eq "") || ($u_base  eq "") ){

		$cu_del++ if ( ($cu_base eq "-") || ($cu_base eq "") );
		$u_del++ if ( ($u_base  eq "-") || ($u_base  eq "") );
		next;
	    }
	    elsif ( ($cu_base ne "-") && ($u_base  ne "-") && ($cu_base ne $u_base) ) {

                $diff++;

	    }
	}

   if ($diff == 0) {
       if ($cu_del <= $u_del) {    

	   delete $unique_seq{$id} if $unique_seq{$id} ;  
 
       }
       elsif ($cu_del > $u_del) {

	   delete $unique_seq{$ID} if $unique_seq{$ID};
	   
 
           last;
       }
           
   }
	
	  		   
  }
} 

%essential = reverse%unique_seq;


$i = 0;

$n = 0;
foreach $ID (sort keys %ID_seq) {
  
    $seq_length = length($ID_seq{$ID});
 
    $m = 0;
    $n++;
 

    foreach $seq (keys  %essential) {
  

	next if ($essential{$seq} =~ /$ID/);  
   
	$diff = 0;
        $m++;

        next if (length($seq) != $seq_length); 
	
	for ($i = 0; $i <= ($seq_length -1); $i++) {
	    $cu_base = substr($ID_seq{$ID}, $i, 1);
	    $u_base = substr($seq, $i, 1);

	    if ( ($cu_base ne "-") && ($u_base ne "-") && ($cu_base ne $u_base) ) {
		
                $diff++;           
   
	    }
            
	}

        if ($diff == 0) {
	    if ($essential{$seq} !~ $ID) {
		$ID =~ s/>//;
		$essential{$seq} = $essential{$seq}."+".$ID;
                $found++;               
                last;   
            }
	    
	}
	 
  }
} 


my %new_hash = reverse %essential;

$pattern_sum = 0; 
$var = 0;
$i = 1;

foreach $ID (sort keys %new_hash) {
     
    @id_list = split(/\+/, $ID);   
  


    $id_list_size = scalar@id_list;

    $seq_with_DRM += $id_list_size;  

    $var++;

    $new_hash{$ID} =~ s/(%%)/\t/g; 

    $new_hash{$ID} =~ s/%//g;
    $new_hash{$ID} =~ s/=/\t/g;   
   @DRM_pattern = split(/\t/, $new_hash{$ID});
   $DRM_list = "@DRM_pattern";
    $DRM_list =~ s/^\t+|^\s+//;   
    
    $DRM_list =~ s/\t+|\s+/,/g;
 
    next unless ($DRM_list); 
    $pattern_sum++;
    
    $DRM_list =~ s/(,,)/,/g;
    $DRM_list =~ s/^,//;

    $each_pattern_percent = sprintf ("%2.3g", ($id_list_size/$input_seq_count)*100);   
    
 print OUT "$DRM_list\t$id_list_size\t$each_pattern_percent";  

   foreach $m (@DRM_pattern) {	
	$m =~ s/\s+|\t//g;
	if ($m) {
	    print OUT "+\t";                         
      
        }
	else {
	    print OUT "\t";

        }
       
    }	
 print OUT "\n";

  foreach $mut(@DRM_pattern) { 

      $mut =~ s/\s+|\t//g;
      $mut =~ s/,//g;
      next unless ($mut);
      $mut_count{$mut} += $id_list_size;
 
  }
   
    $i++;


}

$total_pattern_percent  = sprintf("%2.3g", ($seq_with_DRM/$input_seq_count)*100);  
print OUT "total (% based on $input_seq_count total input seq)\t $seq_with_DRM ($total_pattern_percent)\tDRM total (%)";


foreach $mut (@found_mut) {

    next unless ($mut);
    $mut =~ s/\s+|\t//g;
    $DRM_percent = sprintf("%2.3g", ($mut_count{$mut}/$input_seq_count)*100);
    print  OUT "\t$mut_count{$mut} ($DRM_percent)";
}

print "summary (1) # of DRM patterns: $pattern_sum, (2) seq_with_DRM:  $seq_with_DRM\n";

system ("rm -f $mutation_file");


######################################################################
sub ORGANIZEMUTS{

    my ($str, $longest) = @_;
    my ($length, $mut, $above, $here, $higher,$m, @str);

    unless ($str =~ /\S/ ) {
        return "";
    }


    @str = split(" ", $str);
    $length = scalar @str;

   
    if ($length > $$longest) {
        $$longest = $length;
    }


    for ($m=0; $m <$ length-1; $m++) {
        $mut = $str[$m];
        $above = $str[$m+1];
        $here = GETNUM($mut);
        $higher = GETNUM($above);

	if (  ($here > $higher) || (  ($here == $higher) && (ord(substr($mut, length($mut)-1 ,1)) > ord(substr( $above, length($above)-1, 1) )  ) ) ) {
            $str[$m] = $above;
            $str[$m+1] = $mut;
            $m = -1;
        }

    }
    $str = join(" ", @str);
    return $str;
}

######################################################################
#Gets the protein number from the formatted style
#Example of format:  "L36P"
#Returns:            36
######################################################################

sub GETNUM {
    my ($str) = @_;
    $str =~ /^\w(\d+)\w?/;
    return $1;
    my ($length);
    if (length($str) == 5) {
        return substr($str, 1, 3);
    } elsif (length($str) == 4) {
        return substr($str, 1, 2);
    } else {
        die "Error in MutTable::GETNUM()\nInvalid mutation ID\nScript Halted";
    }
}

######################################################################
sub FASTA_OUT {

   my($allmutdatabase1,$DRM1, $input, $seq_count) = @_;

   my ($mut_type, $id, @mutdb, @mutations, %DRM_hash_hash, $mut, $mut_seq, $found); 

   @mutdb = split(/\t|\s+/,$allmutdatabase1);

   print OUT "mutation_types\t";
   foreach $mut_type (@mutdb) {
       print  OUT "$mut_type\t";
   }
   print OUT "\n";

   my %DRM1 = %$DRM1;

   foreach $id (sort keys  %DRM1) {

     
       @mutations = ();
       @mutations = split(/\s+/, $DRM1{$id});
       foreach $mut(@mutations) {

	  $DRM_hash_hash{$id}{$mut} = "$mut%"; 
 
       } 
       
   }

   
   $found = 0;
  foreach $id (sort keys  %DRM_hash_hash) {
  
	  print OUT ">$id\n";
   
	  foreach $mut_type (@mutdb) {
	      $mut_type =~ s/\s+//g;

	      if ($DRM_hash_hash{$id}{$mut_type} ) {
		  

		   if ($found == 0) {
		       print OUT "$DRM_hash_hash{$id}{$mut_type}\%";
                   }
		   if ($found ==1) {
       
                       print OUT ",$DRM_hash_hash{$id}{$mut_type}\%"; 
                   }
		  $mut_seq++;
		   $found = 1;
		  }

	      else {
		 
		  print OUT "=%";  
		  $found = 0;
	      }
	  }

       print OUT "\n";
  }

   print OUT "total_sequences: $seq_count\t mutations_only: $mut_seq\n";

} 
