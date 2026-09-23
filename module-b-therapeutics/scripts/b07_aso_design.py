from Bio import Entrez, SeqIO
import RNA  # this is the import name for the ViennaRNA package
import pandas as pd
import re

Entrez.email = "shahzadmuneeb15@gmail.com"

# --- Step 1: pull the real PLP1 mRNA sequence (NM_000533.5) ---
print("Fetching PLP1 mRNA (NM_000533.5) from NCBI...")
handle = Entrez.efetch(db="nucleotide", id="NM_000533.5", rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
seq = str(record.seq).upper().replace("T", "U")  # DNA -> RNA alphabet
print(f"Fetched {record.description}, length {len(seq)} nt")

# --- Step 2: fold the whole transcript to get base-pair probabilities ---
print("Computing partition function (this takes a few seconds)...")
fc = RNA.fold_compound(seq)
(mfe_struct, mfe) = fc.mfe()
fc.exp_params_rescale(mfe)
fc.pf()
bpp = fc.bpp()  # base-pair probability matrix (1-indexed internally)

n = len(seq)
unpaired_prob = [1.0] * (n + 1)
for i in range(1, n + 1):
    paired = sum(bpp[i][j] for j in range(1, n + 1) if j != i and bpp[i][j] > 0)
    unpaired_prob[i] = max(0.0, 1.0 - paired)

# --- Step 3: tile candidate 18-20mer ASO target sites ---
WINDOW = 18
candidates = []
for start in range(0, n - WINDOW):
    site = seq[start:start + WINDOW]
    dna_site = site.replace("U", "T")  # ASOs are reported as DNA sequence

    gc_content = (site.count("G") + site.count("C")) / WINDOW

    # avoid long homopolymer runs (poly-A/T/G/C reduces specificity and synthesis quality)
    has_homopolymer = bool(re.search(r"(.)\1{3,}", dna_site))

    # mean accessibility across the window (higher = less structured = better target)
    accessibility = sum(unpaired_prob[start + 1 : start + WINDOW + 1]) / WINDOW

    candidates.append({
        "position": start + 1,
        "sequence_dna": dna_site,
        "gc_content": round(gc_content, 2),
        "accessibility": round(accessibility, 3),
        "has_homopolymer": has_homopolymer,
    })

df = pd.DataFrame(candidates)

# --- Step 4: filter to reasonable ASO design rules ---
filtered = df[
    (df["gc_content"] >= 0.40) & (df["gc_content"] <= 0.60) &
    (~df["has_homopolymer"])
].copy()

filtered = filtered.sort_values("accessibility", ascending=False)

filtered.to_csv("../results/aso_targets.csv", index=False)

print(f"\n{len(df)} total windows scanned, {len(filtered)} passed GC/homopolymer filters")
print("\nTop 15 candidate ASO target sites (ranked by predicted accessibility):")
print(filtered.head(15).to_string(index=False))