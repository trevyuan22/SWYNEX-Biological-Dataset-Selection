import json

def extract(filename, label):
    print(f"--- {label} ---")
    try:
        with open(filename, encoding="utf-8") as f:
            for line in f:
                d = json.loads(line)
                stats = d.get("assemblyStats", {})
                asm = d.get("assemblyInfo", {})
                
                # FIX: Look for 'accession' at the top level first, then fall back to 'assemblyAccession'
                accession = d.get("accession") or asm.get("assemblyAccession") or "Unknown"
                
                print(f"Accession: {accession}")
                print(f"Total bp : {stats.get('totalSequenceLength')}")
                print(f"GC %     : {stats.get('gcPercent')}")
                print(f"Scaffold N50: {stats.get('scaffoldN50')}")
                print(f"Contig N50  : {stats.get('contigN50')}")
                print(f"Scaffolds   : {stats.get('numberOfScaffolds')}")
                print(f"Contigs     : {stats.get('numberOfContigs')}\n")
    except FileNotFoundError:
        print(f"Error: Could not find {filename}\n")

extract("primary_assembly_data_report.jsonl", "PRIMARY ASSEMBLY")
extract("alternate_assembly_data_report.jsonl", "ALTERNATE HAPLOTYPE")