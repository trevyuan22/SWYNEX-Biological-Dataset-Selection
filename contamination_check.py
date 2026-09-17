from Bio import SeqIO

suspicious = []
for rec in SeqIO.parse("primary.fna", "fasta"):
    L = len(rec.seq)
    if L == 0: continue
    gc = (rec.seq.count("G") + rec.seq.count("C")) / L * 100
    # Flag short contigs with anomalous GC content
    if L < 10000 and (gc < 25 or gc > 45):
        suspicious.append((rec.id, L, round(gc, 2)))

print(f"Flagged {len(suspicious)} suspicious contigs in the primary assembly.")
for sid, L, gc in suspicious[:20]:
    print(f"  {sid}  length={L}  GC={gc}%")