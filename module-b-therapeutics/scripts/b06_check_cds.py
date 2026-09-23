from Bio import Entrez, SeqIO

Entrez.email = "shahzadmuneeb15@gmail.com"
handle = Entrez.efetch(db="nucleotide", id="NM_000533.5", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")

for feature in record.features:
    if feature.type == "CDS":
        cds_len = len(feature)
        print(f"CDS location: {feature.location}")
        print(f"CDS length: {cds_len} bp ({cds_len // 3 - 1} amino acids)")
        print(f"Product: {feature.qualifiers.get('product')}")