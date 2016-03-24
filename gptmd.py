# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "anthony"
__date__ = "$Oct 29, 2015 1:25:17 PM$"

import sys
import os.path
import optparse
from lxml import etree as et

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

def read_and_condense_xml(xml_file_name):
    if xml_file_name != None:
        try:
            reference_xml = os.path.abspath(xml_file_name)
            refXml = open(reference_xml, 'r')
            p = et.XMLParser(remove_blank_text=True) #required for pretty additions
            db = et.parse(refXml, p)
            root = db.getroot()
            for entry in root:
                for element in entry:
                    if element.tag not in [UP+'protein',UP+'accession',UP+'name',UP+'gene',UP+'organism',UP+'proteinExistence',UP+'depth',UP+'sequence',UP+'feature',UP+'dbReference']:
                        entry.remove(element)
                    elif element.get('type') != 'modified residue' and element.tag == UP+'feature': entry.remove(element)
                    elif element.get('type') != 'Ensembl' and element.tag == UP+'dbReference': entry.remove(element)
                    elif element.tag == UP+'organism':
                        for field in element:
                            if field.tag != UP+'name':
                                element.remove(field)
                    elif element.tag == UP+'protein':
                        for name in element:
                            if name.tag != UP+'recommendedName':
                                element.remove(name)
                    else:
                      continue
            return db
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(2)

def read_ptm_db(ptm_database):
    if ptm_database != None:
        try:
            ptm_database = os.path.abspath(ptm_database)
            ptm_database = open(ptm_database, 'r')
            id, tg, mm, ptm_types, ptm_masses = None, None, None, {}, {}
            for line in ptm_database:
                line = line.split()
                if line[0] != '//':
                    if id != None and id not in ptm_types: ptm_types[id] = tg
                    if mm != None and mm not in ptm_masses: ptm_masses[mm] = id
                    id, tg, mm = None, None, None
                elif line[0] == 'ID': id = line[1]
                elif line[0] == 'TG': tg = line[1]
                elif line[0] == 'MM': mm = float("%.3f" % line[1]) # round to 3rd place after decimal
            ptm_database.close()
            return ptm_types, ptm_masses
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(2)

def enter_modification(seq_elements, prot_seq, prot_position, ptm_type):
    new_feature = et.Element(UP +'feature', type="modified residue", description = ptm_type, evidence = "3")
    et.SubElement(et.SubElement(new_feature, UP + 'location'), UP+'position', position = str(prot_position))
    for element in seq_elements:
        if element.text.replace('\n','').replace('\r','') == prot_seq:
            while True:
                element = element.getprevious()
                if element.tag != UP + 'feature': 
                    element.addnext(new_feature)
                    break
                else:
                    thisPosition = int(element.find('.//'+UP+'position').get('position'))
                    if prot_position > thisPosition: 
                        element.addnext(new_feature)
                        break
                    elif prot_position == thisPosition:
                        if element.get('description') <= ptm_type: element.addnext(new_feature)
                        else: element.addprevious(new_feature)
                        break
            break
    
def add_open_search_results(db, psms, ptm_types, ptm_masses):
    root = db.getroot()
    mod_mass_tolerance = 0.02
    sequences = {}
    sequence_elements = root.findall('.//' + UP + 'sequence')
    for seq in sequence_elements:
        accession = seq.getparent().find(UP+'accession').text
        sequences[accession] = seq.text.replace('\n','').replace('\r','')
    
    nterm_acetyls = dict()
    nterm_acetyls['A']='N-acetylalanine'
    nterm_acetyls['C']='N-acetylcysteine'
    nterm_acetyls['E']='N-acetylglutamate'
    nterm_acetyls['G']='N-acetylglycine'
    nterm_acetyls['M']='N-acetylmethionine'
    nterm_acetyls['S']='N-acetylserine'
    nterm_acetyls['T']='N-acetylthreonine'
    nterm_acetyls['V']='N-acetylvaline'
    nterm_acetyls['D']='N-acetylaspartate'
    
    if psms != None:
        #try:
        psms = os.path.abspath(psms)
        psms = open(psms, 'r')
        for line in psms:
            if line.startswith('Filename'): continue

            line = line.split('\t')
            isTarget = line[26] in ['TRUE','True','true']
            q_value = float(line[30])
            protein_description, base_peptide_sequence, start_residue,  = line[13], line[12], int(line[14])
            precursor_mass_error = float("%.3f" % float(line[18]))

            if isTarget and q_value <= 100 and protein_description.find('DECOY') == -1:
                protein_description = protein_description.split('|')
                accession = protein_description[1]
                protein_sequence = sequences[accession]

                possible_precursor_mass_errors = [m for m in ptm_masses if precursor_mass_error > m - mod_mass_tolerance and precursor_mass_error < m + mod_mass_tolerance]
                if len(possible_precursor_mass_errors) > 0:
                    for pme in possible_precursor_mass_errors:
                        for ptm_type in ptm_masses[pme]:
                            for mod_aa in ptm_types[ptm_type]:
                                position_in_peptide = base_peptide_sequence.find(mod_aa)
                                while position_in_peptide != -1:
                                    position_in_protein = start_residue + position_in_peptide # 1-indexed
                                    #Special case is for N-terminal acetylated peptides
                                    if ptm_type.find('N-acetyl') >= 0 and position_in_protein == 2 and protein_sequence[position_in_peptide:position_in_peptide + 1] == mod_aa: 
                                        enter_modification(sequence_elements, protein_sequence, position_in_protein, ptm_type)
                                    elif base_peptide_sequence[position_in_peptide:position_in_peptide + 1] == mod_aa:
                                        if sequence[0] != "X": enter_modification(sequence_elements, protein_sequence, position_in_protein + 1, ptm_type)
                                        else: enter_modification(sequence_elements, protein_sequence, position_in_protein, ptm_type)
                                    position_in_peptide = base_peptide_sequence.find(mod_aa, position_in_peptide)

                #Special case: N-terminal acetylation
                if precursor_mass_error > -89.029921 - mod_mass_tolerance and precursor_mass_error < -89.029921 + mod_mass_tolerance and start_residue == 1:
                    second_aa = base_peptide_sequence[1]
                    if second_aa in nterm_acetyls and protein_sequence[1] == second_aa: 
                        enter_modification(sequence_elements, protein_sequence, 2, nterm_acetyls[second_aa])
        psms.close()
        #except Exception, e:
         #   print >> sys.stderr, "failed: %s" % e
          #  exit(2)
            
def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
      #I/O
    parser.add_option( '-x', '--reference_xml', dest='reference_xml', help='The reference UniProt-XML file.' )
    parser.add_option( '-t', '--ptm_database', dest='ptm_database', help='Database of UniProt PTMs (slightly modified)' )
    parser.add_option( '-s', '--psms', dest='psms', help='Peptide spectral matches tab-separated file from open search mode first-pass' )
    parser.add_option( '-o', '--output', dest='output', help='The output UniProt-XML' )
    (options, args) = parser.parse_args()
    
    if options.output != None:
        outF = os.path.abspath(options.output)
        outF = open(outF, 'w')
    else:
        print >> sys.stderr, "failed: no ouput file specified with -o or --output tag"
        exit(2)
    
    db = read_and_condense_xml(options.reference_xml)
    ptm_types, ptm_masses = read_ptm_db(options.ptm_database)
    add_open_search_results(db, options.psms, ptm_types, ptm_masses)
    db.write(outF, pretty_print=True)
    
if __name__ == "__main__" : __main__()
