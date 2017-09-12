# GPTMD We've Moved!!!!
# GPTMD is now fully integrated into [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus). All users should MetaMorpheus.

The original G-PTM-D program that augmented a UniProt XML database with PTMs discovered using an open mass search in [Morpheus](http://cwenger.github.io/Morpheus/)

## New and Improved
The manuscript describing this work entitled "Global Post-translational Modification Discovery" can be found [here.](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00034) 
Software provided below will enable the user to reproduce that work and also create custom databases.
However, we have made major improvements in the process since that time. We direct your attention
to [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus) - single integrated software package for bottom-up proteomics that searches, calibrates 
and discovers new PTMs.

New users are STRONGLY encouraged to begin using [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus), rather than use the former combination of 
Morpheus and G-PTMD. 


# You should really use GPTMD in [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus) because it is much faster and much more accurate!


## General Overview

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


## General Requirements

The following files must be present in the folder with the executable. If not, they are automatically downloaded (to update a file to a newer version, delete it, and the application will download a new version). This is the only network usage by the application. 

* uniprot.xml: A UniProt reference database in .xml format

  [http://www.uniprot.org](http://www.uniprot.org)

* ptmlist.txt: A PTM library
 
  [http://www.uniprot.org/docs/ptmlist.txt](http://www.uniprot.org/docs/ptmlist.txt) 


## Operating System Requirements and Usage

### Windows Perl Version
#### System
#### Usage

### Windows C# Version
#### You should really use GPTMD in [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus) because it is much faster and much more accurate!

### Linux Python Version
#### System
- 8 GB of RAM is recommended
- python v2.7.10 (64 bit)
See https://www.python.org/downloads/ for installation instructions.
This includes the "pip" package manager.
- This program uses lxml, a package for interpreting XML databases. lxml can be installed using the command: pip install lxml. See http://lxml.de for additional installation instructions.
- If you encounter errors installing lxml, we recommend 
trying an alternate package manager, such as Canopy, which can be found 
[here](https://www.enthought.com/products/canopy/).

#### Usage

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

	>python gptmd.py -x ../uniprot.xml -t ../sub_ptmlist_regular.txt -s ../PSMs.tsv -o ../test_output.xml

### Relevant Manuscripts

* [Global Post-translational Modification Discovery--J. Proteome Res., 2016 Just Accepted](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00034)

* [Global Identification of Protein Post-translational Modifications in a Single-Pass Database Search--J. Proteome Res., 2015, 14 (11), pp 4714–4720](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00599)

* [A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra--J. Proteome Res., 2013, 12 (3), pp 1377–1386](http://pubs.acs.org/doi/abs/10.1021/pr301024c)


## License

The software is currently released under the [GNU GPLv3](http://www.gnu.org/licenses/gpl.txt).

Copyright 2016 [Lloyd M. Smith Group](http://smith.chem.wisc.edu/).
