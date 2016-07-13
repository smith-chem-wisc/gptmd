#GPTMD

Software to augment a UniProt XML database with PTMs discovered using [Morpheus](github.com/cwenger/Morpheus)

##General Overview

G-PTMD is a tool used to expand the scope of peptide identification to include 
specific post-translational modifications.  Currently, identifying peptides with 
post-translational modifications relies on the variable toggling of modifications 
on all the residues that can be modified or the documentation of a specific 
modification in a database.  The first method is computationally expensive and 
wasteful while databases are often incomplete.  Thus, G-PTMD addresses the weaknesses 
of both methods to improve protein idenfication results.  
  
The purpose of G-PTMD is to build a new proteome reference .xml database by annotating
an existing reference database using a set of peptide spectral matches (PSMs).  
The PSMs are obtained by running a Morpheus search on a .raw or .mzml file of the 
experimentally obtained mass spectrometry data.   The program identifies the peptide
spectral matches that have a mass shift indicative of certain post-translational 
modifications.  For example, if the peptide “PEPTIDE” was identified in the PSMs 
file with a mass shift of 79.966 Da, the corresponding protein entry in the 
reference .xml database would be granted a new feature consisting of a phosphorylated 
threonine at the appropriate position because 79.966 Da is exactly the change in mass 
of a peptide when one residue is phosphorylated.  

Now, when the original .raw or .mzml file is run against this new database in an open 
search, if the phosphorylation on the threonine is present in the sample, the computer 
will be able to recognize and identify it.  Phosphorylation is just one example, and this
technique extends to other post-translational modifications.  The user can add or remove
modifications to the list of post-translational modifications as needed.  This effectively
allows for variable post-translational modifications at targeted positions in the protein 
leading to better search results upon the second pass without incurring huge data costs 
associated with blindly adding variable modifications.


##General Requirements

The following files must be present in the folder with the executable. If not, they are automatically downloaded (to update a file to a newer version, delete it, and the application will download a new version). This is the only network usage by the application. 

* ptmlist.txt: A PTM library
 
  [http://www.uniprot.org/docs/ptmlist.txt](http://www.uniprot.org/docs/ptmlist.txt) 


##Operating System Requirements and Usage

####System
####Usage

###Windows Version
####System
#####

###Python Version
####System
- 8 GB of RAM is recommended
- python v2.7.10 (64 bit)
See https://www.python.org/downloads/ for installation instructions.
This includes the "pip" package manager.
- lxml python packageINstall using the command: pip install lxml
Or see http://lxml.de for installation instructions.
- If you encounter errors installing the package, we recommend 
trying an alternate package manager, such as Canopy, which can be found 
[here](https://www.enthought.com/products/canopy/).

####Usage

Options:

	-h, --help			show this help message and exit
	
	-x REFERENCE_XML, --reference_xml=REFERENCE_XML
						The reference UniProt-XML file.  New PTM features are
						appended to this database to generate the output
						UniProt-XML protein database.
						
	-t PTM_DATABASE, --ptm_database=PTM_DATABASE
						Slightly modified database of UniProt PTMs.  This file
						determines which types of PTMs are included.
						
	-s PSMS, --psms=PSMS  Peptide spectral matches tab-separated.  This file is
						from first-pass open search and contains the mass
						shifts that correspond to PTMs.
						
	-o OUTPUT, --output=OUTPUT
						Output file path.  Outputs a UniProt-XML file.


Example Command Line:

\WorkingDirectory>bpython gptmd.py -x ../uniprot.xml -t ../sub_ptmlist_regular.txt
-s ../PSMs.tsv -o ../test_output.xml

####Notes

* 

###Relevant Manuscripts

* [GPTMD manuscript in revision](http://pubs.acs.org/journal/jprobs)

* [Global Identification of Protein Post-translational Modifications in a Single-Pass Database Search--J. Proteome Res., 2015, 14 (11), pp 4714–4720](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00599)

* [A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra--J. Proteome Res., 2013, 12 (3), pp 1377–1386](http://pubs.acs.org/doi/abs/10.1021/pr301024c)


## License

The software is currently released under the [GNU GPLv3](http://www.gnu.org/licenses/gpl.txt).

Copyright 2016 Lloyd M. Smith Group.