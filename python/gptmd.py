# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "Anthony J. Cesnik"
__date__ = "$Oct 29, 2015 1:25:17 PM$"

import sys
import os.path
import optparse
from lxml import etree as et
import utility

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

MOD_MASS_TOLERANCE = 0.02
PTMLIST_HEADER_ABBREVS = ['ID', 'AC', 'FT', 'TG', 'PA', 'PP', 'CF', 'MM', 'MA', 'LC', 'TR', 'KW', 'DR', '//', 'ZZ']
MIN_FDR_FIRST_PASS = 100

def condense_xml_entry(entry):
    for element in entry:
        if element.tag not in [UP+'protein',UP+'accession',UP+'name',UP+'gene',UP+'organism',UP+'proteinExistence',UP+'depth',UP+'sequence',UP+'feature',UP+'dbReference']:
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
    new_feature = et.Element(UP+'feature', type="modified residue", description=ptm_type, evidence="3")
    et.SubElement(et.SubElement(new_feature, UP+'location'), UP+'position', position=str(prot_position))
    for element in seq_elements:
        if element.text.replace('\n', '').replace('\r', '') == prot_seq:
            while True:
                element = element.getprevious()
                if element.tag != UP+'feature':
                    element.addnext(new_feature)
                    break
                else:
                    this_position = int(element.find('.//'+UP+'position').get('position'))
                    if prot_position > this_position:
                        element.addnext(new_feature)
                        break
                    elif prot_position == this_position:
                        if element.get('description') <= ptm_type: element.addnext(new_feature)
                        else: element.addprevious(new_feature)
                        break
            break

def equals_within_tolerance(x, value, tolerance):
    return x > value - tolerance and x < value + tolerance

def keep_psm(line, ptm_masses):
    line = line.split('\t')
    is_target = line[26] in ['TRUE', 'True', 'true']
    q_value = float(line[30])
    precursor_mass_error = float("%.3f" % float(line[18]))

    if not is_target or q_value > MIN_FDR_FIRST_PASS: return False

    delta_m_of_interest = False
    for dm in ptm_masses:
        if equals_within_tolerance(precursor_mass_error, dm, MOD_MASS_TOLERANCE): return True

    return False

def add_open_search_results(line, sequence_elements, sequences, ptm_types, ptm_masses, nterm_acetyls):
    line = line.split('\t')
    protein_description, base_peptide_sequence, start_residue = line[13], line[12], int(line[14])
    precursor_mass_error = float("%.3f" % float(line[18]))

    protein_description = protein_description.split('|')
    accession = protein_description[1]
    protein_sequence = sequences[accession]

    possible_precursor_mass_errors = [deltaM for deltaM in ptm_masses if equals_within_tolerance(precursor_mass_error, deltaM, MOD_MASS_TOLERANCE)]
    if len(possible_precursor_mass_errors) > 0:
        for pme in possible_precursor_mass_errors:
            for ptm_type in ptm_masses[pme]:
                for mod_aa in ptm_types[ptm_type]:
                    position_in_peptide = base_peptide_sequence.find(mod_aa)
                    while position_in_peptide != -1:
                        # print base_peptide_sequence + ' ' + str(position_in_peptide)
                        position_in_protein = start_residue + position_in_peptide # 1-indexed
                        position_is_mod_aa = base_peptide_sequence[position_in_peptide:position_in_peptide + 1] == mod_aa
                        is_n_term_acetyl_site = position_in_protein == 2 and protein_sequence[position_in_peptide:position_in_peptide + 1] == mod_aa and ptm_type.find('N-acetyl') >= 0

                        if is_n_term_acetyl_site: enter_modification(sequence_elements, protein_sequence, position_in_protein, ptm_type)
                        elif position_is_mod_aa:
                            if base_peptide_sequence[0] != "X":
                                enter_modification(sequence_elements, protein_sequence, position_in_protein + 1, ptm_type)
                            else:
                                enter_modification(sequence_elements, protein_sequence, position_in_protein, ptm_type)
                        position_in_peptide = base_peptide_sequence.find(mod_aa, position_in_peptide + 1) # Increment the start of the search

    #Special case: N-terminal acetylation
    is_n_term_acetyl_site = start_residue == 1 and equals_within_tolerance(precursor_mass_error, -89.029921, MOD_MASS_TOLERANCE)
    if is_n_term_acetyl_site:
        second_aa = base_peptide_sequence[1]
        if second_aa in nterm_acetyls and protein_sequence[1] == second_aa:
            enter_modification(sequence_elements, protein_sequence, 2, nterm_acetyls[second_aa])
            "entering modification " + nterm_acetyls[second_aa]
            
def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
      #I/O
    parser.add_option( '-x', '--reference_xml', dest='reference_xml', help='The reference UniProt-XML file.' )
    parser.add_option( '-t', '--ptm_database', dest='ptm_database', help='Database of UniProt PTMs (slightly modified)' )
    parser.add_option( '-s', '--psms', dest='psms', help='Peptide spectral matches tab-separated file from open search mode first-pass' )
    parser.add_option( '-o', '--output', dest='output', help='The output UniProt-XML' )
    (options, args) = parser.parse_args()

    ####OUTPUT: new xml database
    if options.output != None:
        outF = os.path.abspath(options.output)
        outF = open(outF, 'w')
    else:
        print >> sys.stderr, "failed: no ouput file specified with -o or --output tag"
        exit(2)

    #Parse the reference XML
    try:
        reference_xml = os.path.abspath(options.reference_xml)
        refXml = open(reference_xml, 'r')
        p = et.XMLParser(remove_blank_text=True)  # required for pretty additions
        db = et.parse(refXml, p)
        root = db.getroot()
        for entry in root: condense_xml_entry(entry)
    except Exception, e:
        print >> sys.stderr, "failed: no ouput file specified with -o or --output tag. %s" % e
        exit(2)

    #Parse the ptmlist
    try:
        aa_dict = utility.get_aa_dict()
        ptm_types, ptm_masses = {}, {}
        id, tg, mm = None, None, None

        ptm_database = os.path.abspath(options.ptm_database)
        ptm_database = open(ptm_database, 'r')
        for line in ptm_database:
            line = line.strip()
            if line[:2] not in PTMLIST_HEADER_ABBREVS: continue # Ignore the header information
            line = line.split('\t')
            if len(line) == 1: line = ''.join(line).split('   ') #original ptmlist.txt has spaces
            if line[0] in ['//', 'ZZ']:
                if id != None:
                    if id not in ptm_types: ptm_types[id] = [tg]
                    else: ptm_types[id].append(tg)
                if mm != None:
                    if mm not in ptm_masses: ptm_masses[mm] = [id]
                    else: ptm_masses[mm].append(id)
                print "just added " + id + " and " + tg
                id, tg, mm = None, None, None
            elif line[0] == 'ID': id = line[1]
            elif line[0] == 'TG': tg = aa_dict[line[1].strip('.')] # original ptmlist ends lines with periods
            elif line[0] == 'MM': mm = float("%.3f" % float(line[1]))  # round to 3rd place after decimal
        ptm_database.close()
    except Exception, e:
        print >> sys.stderr, "Failed to open the UniProt PTM database. %s" % e
        exit(2)

    print ptm_types
    print ptm_masses

    #Process the first-pass psms
    try:
        psms = os.path.abspath(options.psms)
        linect = sum(1 for line in open(psms))
        psms = open(psms, 'r')
        psms_list = []

        #Preprocess psms
        for i, line in enumerate(psms):
            if i % 10000 == 0: print "psms preprocessing line " + str(i) + " of " + str(linect)
            if line.startswith('Filename'): continue
            if keep_psm(line, ptm_masses): psms_list.append(line)

        #Add open search results
        nterm_acetyls = utility.get_nterm_dict()
        sequences = {}
        sequence_elements = root.findall('.//'+UP+'sequence')
        for seq in sequence_elements:
            accession = seq.getparent().find(UP+'accession').text
            sequences[accession] = seq.text.replace('\n','').replace('\r','')
        for i, line in enumerate(psms_list):
            if i % 100 == 0: print "psm " + str(i) + " of " + str(len(psms_list))
            if line.startswith('Filename'): continue
            add_open_search_results(line, sequence_elements, sequences, ptm_types, ptm_masses, nterm_acetyls)
        psms.close()
    except Exception, e:
        print >> sys.stderr, "Failed to open the PSMS.tsv list. %s" % e
        exit(2)

    #Write the database
    db.write(outF, pretty_print=True)
    outF.close()
    
if __name__ == "__main__" : __main__()