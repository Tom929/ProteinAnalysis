from Bio import ExPASy
from Bio import SwissProt

accessions = ["O23729", "O23730", "O23731"]
records = []

for accession in accessions:
     handle = ExPASy.get_sprot_raw(accession)
     try:
         record = SwissProt.read(handle)
     except ValueException:
         print("WARNING: Accession %s not found" % accession)
     records.append(record)

