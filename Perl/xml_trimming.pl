#!/usr/bin/perl
 
# use module
use XML::Simple;

print "This can take a long time.\n";
print "Better go get a cup of coffee.\n";
#####   START - Read Database of UniProt PTMs              #####
my %PTM_type_hash = (); #keys are PTM ID full name and values are PTM target full name
my %PTM_mass_hash = (); #keys are PTM mass shift rounded to three decimal places and values are PTM ID full name (e.g. Phosphoserine)
my $PTM_ID; # PTM ID full name
my $PTM_Target_Name; # PTM target full name
my $PTM_Mass_Shift;
  
#####   START - Read XML Search Database from UniProt            #####
my %observed_uniprot_ptms = (); #key is accession and value is UniProt PTM and location
my %global_sequence_data = (); #key is accession and value detailed sequence data (length, mass, etc.)
 
my @accessions_with_PTMs = ();
my $xml_deref='uniprot_trimmed.xml';#####################              MAIN OUTPUT FILE        ########################
open(MYOUTFILE, ">$xml_deref"); #open for write, overwrite

# create object
$xml = new XML::Simple (ForceArray => 1,KeyAttr=>[]);
 
# read XML file
$data = $xml->XMLin("uniprot-all.xml");#####################          XML DATABASE FILE ########################


print MYOUTFILE '<?xml version=\'1.0\' encoding=\'UTF-8\'?>', "\n";
print MYOUTFILE '<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">',"\n";
foreach $e (@{$data->{entry}})
              {
              my %local_entry_data_hash = ();
              $local_entry_data_hash{dataset}= $e->{dataset}, "\t";
              $local_entry_data_hash{entry_created}= $e->{created}, "\t";
              $local_entry_data_hash{entry_modified}= $e->{modified}, "\t";
              $local_entry_data_hash{entry_version}= $e->{version}, "\n";
              my $first_accession = $e->{accession}->[0]; #could go back later and add the other accession numbers???
              $local_entry_data_hash{name}= $e->{name}->[0];
              $local_entry_data_hash{fullName}= $e->{protein}->[0]->{recommendedName}->[0]->{fullName}->[0];
			  if (defined $local_entry_data_hash{fullName}) 
				  {
				  if(ref($local_entry_data_hash{fullName}) eq "HASH")
					  {
					  $local_entry_data_hash{fullName} = $first_accession;
					  }  
				  }
			  else  
				  {
				  $local_entry_data_hash{fullName} = $first_accession;
				  }
			  $local_entry_data_hash{gene}=$e->{gene}->[0]->{name}->[0]->{content};
              $local_entry_data_hash{organism_name}=$e->{organism}->[0]->{name}->[0]->{content}, "\n";
              $local_entry_data_hash{dbReference_id}=$e->{organism}->[0]->{dbReference}->[0]->{id};
              $local_entry_data_hash{dbReference_type}=$e->{organism}->[0]->{dbReference}->[0]->{type};
              $local_entry_data_hash{proteinExistence}=$e->{proteinExistence}->[0]->{type};
              $local_entry_data_hash{length}=$e->{sequence}->[0]->{length};
              $local_entry_data_hash{mass}=$e->{sequence}->[0]->{mass};
              $local_entry_data_hash{checksum}=$e->{sequence}->[0]->{checksum};
              $local_entry_data_hash{modified}=$e->{sequence}->[0]->{modified};
              $local_entry_data_hash{version}=$e->{sequence}->[0]->{version};
              (my $aminoAcidSequence =$e->{sequence}->[0]->{content}) =~ s/\s//g; #gets sequence with no whitespace
              $local_entry_data_hash{content}=$aminoAcidSequence;
 
              my $feature_array_ref = $e->{feature};
              foreach my $feature (@$feature_array_ref)
                             {
                             if (${$feature}{type}eq'modified residue')
                                           {
                                           push (@{$observed_uniprot_ptms{$first_accession}{${$feature}{location}->[0]->{position}->[0]->{position}}},${$feature}{description});
                                           }
                             }
              $global_sequence_data{$first_accession}=\%local_entry_data_hash;
              }
#####   END - Read XML Search Database from UniProt #####
 
#####   START - Compile Array of Accession Numbers with PTMs             #####
foreach my $AN (sort keys %observed_uniprot_ptms)
              {
              push (@accessions_with_PTMs,$AN);
              }
#####   END - Compile Array of Accession Numbers with PTMs  #####
 
#####   START - Write New XML Database with Novel PTMs        #####
foreach my $AN (sort keys %global_sequence_data)
	{
		if((index(${$global_sequence_data{$AN}}{content},'Z')==-1)&&(index(${$global_sequence_data{$AN}}{content},'B')==-1))
			{
				print MYOUTFILE '<entry dataset="', ${$global_sequence_data{$AN}}{dataset},'" created="', ${$global_sequence_data{$AN}}{entry_created},'" modified="', ${$global_sequence_data{$AN}}{entry_modified},'" version="', ${$global_sequence_data{$AN}}{entry_version},'">',"\n";
				print MYOUTFILE '<accession>',$AN ,'</accession>',"\n";
				print MYOUTFILE '<name>',${$global_sequence_data{$AN}}{name} ,'</name>',"\n";
				print MYOUTFILE '<protein>',"\n";
				print MYOUTFILE '<recommendedName>',"\n";
				print MYOUTFILE '<fullName>',${$global_sequence_data{$AN}}{fullName},'</fullName>',"\n";
				print MYOUTFILE '</recommendedName>',"\n";
				print MYOUTFILE '</protein>',"\n";
				print MYOUTFILE '<gene>',"\n";
				print MYOUTFILE '<name type="primary">',${$global_sequence_data{$AN}}{gene},'</name>',"\n";
				print MYOUTFILE '</gene>',"\n";
				print MYOUTFILE '<organism>',"\n";
				print MYOUTFILE '<name type="scientific">',${$global_sequence_data{$AN}}{organism_name} ,'</name>',"\n";
				print MYOUTFILE '<dbReference type="', ${$global_sequence_data{$AN}}{dbReference_type},'" id="', ${$global_sequence_data{$AN}}{dbReference_id},'"/>',"\n";
				print MYOUTFILE '</organism>',"\n";
				print MYOUTFILE '<proteinExistence type="', ${$global_sequence_data{$AN}}{proteinExistence},'"/>',"\n";
		
				if (grep /$AN/, @accessions_with_PTMs)
								{
								foreach my $base (sort keys %{ $observed_uniprot_ptms{$AN} })
											{
											my @unique_ptms = do { my %seen; grep { !$seen{$_}++ } @{$observed_uniprot_ptms{$AN}{$base}} };
											foreach my $ptm (@unique_ptms)
															{
															print MYOUTFILE '<feature type="modified residue" description="',$ptm,'" evidence="3">', "\n";
															print MYOUTFILE "\t",'<location>', "\n";
															print MYOUTFILE "\t\t",'<position position="',$base,'"/>', "\n";
															print MYOUTFILE "\t",'</location>', "\n";
															print MYOUTFILE '</feature>', "\n";
															}
											}
								}
		
				print MYOUTFILE '<sequence length="', ${$global_sequence_data{$AN}}{length},'" mass="',${$global_sequence_data{$AN}}{mass},'" checksum="',${$global_sequence_data{$AN}}{checksum},'" modified="',${$global_sequence_data{$AN}}{modified},'" version="',${$global_sequence_data{$AN}}{version},'">',"\n";
				print MYOUTFILE ${$global_sequence_data{$AN}}{content};
				print MYOUTFILE "\n";
				print MYOUTFILE '</sequence>',"\n";
				print MYOUTFILE '</entry>',"\n";
			}
		else
			{
				print "Rejecting $AN because the sequence contains ambiguous amino acids.\n";
			}
	}
 
 
print MYOUTFILE '<copyright>', "\n";
print MYOUTFILE 'Copyrighted by the UniProt Consortium, see http://www.uniprot.org/terms Distributed under the Creative Commons Attribution-NoDerivs License', "\n";
print MYOUTFILE '</copyright>', "\n";
print MYOUTFILE '</uniprot>', "\n";
close MYOUTFILE;
#####   END - Write New XML Database with Novel PTMs            #####