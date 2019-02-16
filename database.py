from Bio import ExPASy
from Bio import SwissProt

#accesses Swissprot database and returns polypeptide sequence
def access_sequence(accession):
    handle = ExPASy.get_sprot_raw(accession)
    try:
        record = SwissProt.read(handle)
    except ValueException:
        print("WARNING: Accession %s not found" % accession)
    return record.sequence

#accesses Swissprot database and returns polypeptide entry_name
def access_name(accession):
    handle = ExPASy.get_sprot_raw(accession)
    try:
        record = SwissProt.read(handle)
    except ValueException:
        print("WARNING: Accession %s not found" % accession)
    return record.entry_name