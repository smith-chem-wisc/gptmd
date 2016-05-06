# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "anthony"
__date__ = "$Nov 3, 2015 8:29:23 PM$"

from lxml import etree as et

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

def compare1(db1, db2):
    db1 = et.parse(db1)
    db2 = et.parse(db2)
    for entry1 in db1.getroot():
        for entry2 in db2.getroot():
            if entry1.tag != UP+'copyright' and entry2.tag != UP+'copyright' and entry2.find(UP+'accession').text == entry1.find(UP+'accession').text:
                mods1, mods2 = entry1.findall(UP+'feature'), entry2.findall(UP+'feature')
                if len(mods1) != len(mods2):
                    print "FAILED"
                    break
                break
                
def compare2(db1, db2):
    db1 = et.parse(db1)
    db2 = et.parse(db2)
    for entry1 in db1.getroot():
        for entry2 in db2.getroot():
            if entry1.tag != UP+'copyright' and entry2.tag != UP+'copyright' and entry2.find(UP+'accession').text == entry1.find(UP+'accession').text:
                mods1 = [(x.get('description'),x.find('.//'+UP+'position').get('position')) for x in entry1.findall(UP+'feature')]
                mods2 = [(x.get('description'),x.find('.//'+UP+'position').get('position')) for x in entry2.findall(UP+'feature')]
                if len(mods1) != len(mods2):
                    print "FAILED"
                    break
                for mod in mods1:
                    if mod not in mods2:
                        print mod if mod != None else "None"
                        print "FAILED"
                        break
                        break
                for mod in mods2:
                    if mod not in mods1:
                        print mod if mod != None else "None"
                        print "FAILED"
                        break
                        break
                break