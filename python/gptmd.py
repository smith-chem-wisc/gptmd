__author__ = "Anthony J. Cesnik, Stefan Solntsev"
__date__ = "$Oct 29, 2015 1:25:17 PM$"

import sys
import os.path
import optparse
from lxml import etree as et
import utility

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None: HTML_NS, "xsi": XSI_NS}
UP = '{'+HTML_NS+'}'

MOD_MASS_TOLERANCE = 0.02
PTMLIST_HEADER_ABBREVS = ['ID', 'AC', 'FT', 'TG', 'PA', 'PP', 'CF', 'MM', 'MA', 'LC', 'TR', 'KW', 'DR', '//']
MIN_FDR_FIRST_PASS = 100

usedAccessionList = []
unusedAccessionList = []
badAAList = []

count = 0

def condense_xml_entry(entry):
    for element in entry:
        if element.tag not in [UP+'protein', UP+'accession', UP+'name', UP+'gene', UP+'organism', UP+'proteinExistence', UP+'depth', UP+'sequence', UP+'feature', UP+'dbReference']:
            entry.remove(element)
        elif element.get('type') != 'modified residue' and element.tag == UP+'feature': entry.remove(element)
        elif element.get('type') != 'Ensembl' and element.tag == UP+'dbReference': entry.remove(element)
        elif element.tag == UP+'organism':
            for field in element:
                if field.tag != UP+'name': element.remove(field)
        elif element.tag == UP+'protein':
            for name in element:
                if name.tag != UP+'recommendedName': element.remove(name)
        else: continue


def enter_modification(seq_elements, prot_seq, prot_position, ptm_type):
    global count
    new_feature = et.Element(UP+'feature', type="modified residue", description=ptm_type, evidence="3")
    et.SubElement(et.SubElement(new_feature, UP+'location'), UP+'position', position=str(prot_position))
    for element in seq_elements:
        if element.text.replace('\n', '').replace('\r', '') == prot_seq:
            while True:
                element = element.getprevious() #Iterate backwards from the sequence element 
                if element.tag != UP + 'feature': # If suddenly found an element that is not a modified residue, add the new one!
                    element.addnext(new_feature)
		    count += 1
                    break
                else:
                    this_position = int(element.find('.//'+UP+'position').get('position'))
                    if prot_position > this_position: # If found one that is before the one to be added, add the new one!
                        element.addnext(new_feature)
		        count += 1
                        break
                    elif prot_position == this_position and  element.get('description') == ptm_type:  # Prevents duplicate feature entries
                        break
            break # Break out of outer for loop, the one that looks for correct sequence


def equals_within_tolerance(x, value, tolerance):
    return x > value - tolerance and x < value + tolerance


def keep_psm(line, ptm_masses):  # So far, there is no fail-safe for blanks lines or other problems.
    line = line.split('\t')
    is_target = line[26] in ['TRUE', 'True', 'true']
    q_value = float(line[30])
    precursor_mass_error = float("%.3f" % float(line[18]))

    if not is_target or q_value > MIN_FDR_FIRST_PASS: return False
    if equals_within_tolerance(precursor_mass_error, -89.0299, MOD_MASS_TOLERANCE): return True  # For Case 1 N-term
    for dm in ptm_masses:
        if equals_within_tolerance(precursor_mass_error, dm, MOD_MASS_TOLERANCE): return True

    return False


def add_open_search_results(line, sequence_elements, sequences, ptm_types, pp_types, ptm_masses, nterm_acetyls):
    global unusedAccessionList, usedAccessionList, nonmatchingSequences, badAAList

    # Extract necessary information from the PSMs line.
    line = line.split('\t')
    protein_description, base_peptide_sequence, start_residue = line[13], line[12], int(line[14])
    precursor_mass_error = float(line[18])
    protein_description = protein_description.split('|')
    accession = protein_description[1]

    nameList = protein_description[2].split(' ')
    length = len(protein_description[2].split(' '))
    if nameList[length-6] == "(Fragment)": rangeValue = 6
    else: rangeValue = 5
    for i in range(rangeValue):  # Strip name field of OS, GN, PE, SV, and (Fragment)
        nameList.pop(length-i-1)
    if accession not in sequences:  # Ensures that the "sequences" dictionary has an entry for the "accession."
        unusedAccessionList.append(accession)
        return

    if any((AA in set('BXZ')) for AA in sequences[accession]):  # Extra caution with "bad" (ambiguous) residues.
        badAAList.append(accession)
        protein_sequence = sequences[accession]
        possiblePeptidePositions = [i for i in range(len(protein_sequence)) if protein_sequence.startswith(base_peptide_sequence, i)]
        if len(possiblePeptidePositions) < 1:  # In case you cannot find the peptide in the protein.
            print "Accession,", accession, ", removed because one peptide from PSMs list could not be found in this protein's sequence."
            unusedAccessionList.append(accession)
            return
        if len(possiblePeptidePositions) == 1:
            start_residue = possiblePeptidePositions[0] + 1
        if len(possiblePeptidePositions) > 1:  # It becomes very messy to deal with multiple possible peptide possitions, so we'll remove them.
            print "Accession,", accession,", removed because ambiguity with the peptide position in this protein's sequence."
            unusedAccessionList.append(accession)
            return

    protein_sequence = sequences[accession]
    usedAccessionList.append(accession)

    # All modifications handled below.
    possible_precursor_mass_errors = [deltaM for deltaM in ptm_masses if equals_within_tolerance(precursor_mass_error, deltaM, MOD_MASS_TOLERANCE)]
    if len(possible_precursor_mass_errors) > 0:
        for pme in possible_precursor_mass_errors:
            for ptm_type in ptm_masses[pme]:
                mod_aa = ptm_types[ptm_type]
		mod_pp = pp_types[ptm_type]
		if mod_pp == 'Anywhere':
                    position_in_peptide = base_peptide_sequence.find(mod_aa)
                    while position_in_peptide != -1:
                        position_in_protein = start_residue + position_in_peptide
                        enter_modification(sequence_elements, protein_sequence, position_in_protein, ptm_type)
                        position_in_peptide = base_peptide_sequence.find(mod_aa, position_in_peptide + 1)  # Increment the start of the search
                elif mod_pp == 'Any N-terminal':
		    if base_peptide_sequence[0] == mod_aa or mod_aa == "Any":
                        position_in_protein = start_residue 
                        enter_modification(sequence_elements, protein_sequence, position_in_protein, ptm_type)
                elif mod_pp == 'Any C-terminal':
		    if base_peptide_sequence[len(base_peptide_sequence)-1] == mod_aa or mod_aa == "Any":
                        position_in_protein = start_residue + len(base_peptide_sequence)-1
                        enter_modification(sequence_elements, protein_sequence, position_in_protein, ptm_type)
		else:
                    print >> sys.stderr, "Failed: unknown mod_pp " + str(mod_pp)
                    exit(2)
		
			



def __main__():
    # Parse Command Line
    parser = optparse.OptionParser()
    # I/O
    parser.add_option( '-x', '--reference_xml', dest='reference_xml', help='The reference UniProt-XML file.  New PTM features are appended to this database to generate the output UniProt-XML protein database.')
    parser.add_option( '-t', '--ptm_database', dest='ptm_database', help='Slightly modified database of UniProt PTMs.  This file determines which types of PTMs are included.')
    parser.add_option( '-s', '--psms', dest='psms', help='Peptide spectral matches tab-separated.  This file is from first-pass open search and contains the mass shifts that correspond to PTMs.')
    # parser.add_option( '-m', '--threads', dest='threads', help='Number of threads to use for adding annotations.')
    parser.add_option( '-o', '--output', dest='output', help='Output file path.  Outputs a UniProt-XML file.')
    (options, args) = parser.parse_args()

    # OUTPUT: new xml database
    if options.output != None:
        outF = os.path.abspath(options.output)
        outF = open(outF, 'w')
    else:
        print >> sys.stderr, "Failed: no output file specified with -o or --output tag."
        exit(2)


    #Iterative Parse the reference XML
    try:
        reference_xml = os.path.abspath(options.reference_xml)
        refXml = open(reference_xml, 'r')
        xml_iter = et.iterparse(refXml, remove_blank_text=True)
        root = et.Element(UP + 'uniprot', nsmap=NAMESPACE_MAP)
        db = et.ElementTree(root)
        for event, el in xml_iter:
            if el.tag == UP + "entry":
                condense_xml_entry(el)
                root.append(el)
    except Exception, e:
        print >> sys.stderr, "Failed: error occurred while trying to parse the reference XML file. %s" %e
        exit(2)

    # Parse the ptmlist
    try:
        aa_dict = utility.get_aa_dict()
        ptm_types, ptm_masses, pp_types = {}, {}, {}
        id, tg, mm, pp = None, None, None, None

        ptm_database = os.path.abspath(options.ptm_database)
        ptm_database = open(ptm_database, 'r')
        for line in ptm_database:
            line = line.strip()
            if line[:2] not in PTMLIST_HEADER_ABBREVS: continue  # Ignore the header information
            line = ''.join(line).split('   ')  # ptmlist.txt has three spaces
            if line[0] == '//':
                if id != None and mm != None:
                    ptm_types[id] = tg
                    pp_types[id] = pp
                    if mm not in ptm_masses: ptm_masses[mm] = [id]
                    else: ptm_masses[mm].append(id)
                print "Just added " + id + " and " + tg
                id, tg, mm, pp = None, None, None, None
            elif line[0] == 'ID': id = line[1]
            elif line[0] == 'TG': tg = aa_dict[line[1].strip('.')]  # Original ptmlist ends lines with periods
            elif line[0] == 'PP': pp = line[1].strip('.') # Original ptmlist ends lines with periods
            elif line[0] == 'MM': mm = float(line[1])  
        ptm_database.close()
    except Exception, e:
        print >> sys.stderr, "Failed: could not open the PTM database. %s" % e
        exit(2)

    print ptm_types
    print ptm_masses

    # Process the first-pass psms
    # try:
    psms = os.path.abspath(options.psms)
    linect = sum(1 for line in open(psms))
    psms = open(psms, 'r')
    psms_list = []

    # Preprocess psms
    for i, line in enumerate(psms):
        if i % 10000 == 0: print "Preprocessing psms line " + str(i) + " of " + str(linect)
        if line.startswith('Filename'): continue
        if keep_psm(line, ptm_masses):
            psms_list.append(line)

    # Add open search results
    nterm_acetyls = utility.get_nterm_dict()
    sequences = {}
    sequence_elements = root.findall('.//'+UP+'sequence')

    for seq in sequence_elements:  # Extract the name and full name from the reference xml to tag a sequence.
        nameFullNameList = []
        accession = seq.getparent().find(UP+'accession').text
        #  Avoid problem with retrieving text from None type.
        if seq.getparent().find(UP+'protein').find(UP+'recommendedName') != None:
            if seq.getparent().find(UP+'name').text != None and seq.getparent().find(UP+'protein').find(UP+'recommendedName').find(UP+'fullName').text != None:
                nameFullNameList.append(seq.getparent().find(UP+'name').text)
                nameFullNameList.append(seq.getparent().find(UP+'protein').find(UP+'recommendedName').find(UP+'fullName').text)
            elif seq.getparent().find(UP+'name').text != None:
                nameFullNameList.append(seq.getparent().find(UP+'name').text)
            elif seq.getparent().find(UP+'protein').find(UP+'recommendedName').find(UP+'fullName').text != None:
                nameFullNameList.append(seq.getparent().find(UP + 'protein').find(UP + 'recommendedName').find(UP + 'fullName').text)
        elif seq.getparent().find(UP+'name').text != None:
            nameFullNameList.append(seq.getparent().find(UP+'name').text)

        #print "adding "+seq.text.replace('\n','').replace('\r','')+" to sequences!!!"+ " with accession " + accession + " nameAndFullName = " + nameAndFullName 
        sequences[accession] = seq.text.replace('\n','').replace('\r','')
	

    for i, line in enumerate(psms_list):
        if i % 1000 == 0: print "Processing psm " + str(i) + " of " + str(len(psms_list))
        if line.startswith('Filename'): continue
        add_open_search_results(line, sequence_elements, sequences, ptm_types, pp_types, ptm_masses, nterm_acetyls)


    psms.close()

    print "Number of new features added: ", count

    # Write the database
    db.write(outF, pretty_print=True)
    outF.close()
    print "Finished writing new database."
if __name__ == "__main__" : __main__()
