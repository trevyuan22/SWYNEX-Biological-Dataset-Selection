**Data Cleaning \& Basic Analysis of *Leptodirus hochenwartii* Genome**



Organism - *Leptodirus hochenwartii*

Primary Assembly - GCA\_947310635.1 (primary.fna)

Alternate Assembly - GCA\_947310685.1 (alternate.fna)

primary\_assembly\_data\_report.jsonl and alternate\_assembly\_data\_report.jsonl — metadata files provided by NCBI

primary\_md5sum.txt and alternate\_md5sum.txt — checksum files for verifying download integrity



Required Python Packages - biopython, pandas, numpy, and matplotlib



NOTE: All assembly data files downloaded from NCBI must be in the same folder





1. Run 'assembly\_stats.py'

   * Computes contiguity and composition statistics
   * Outputs 'stats\_primary.txt' and 'stats\_alternate.txt'
2. Run 'python assembly\_stats.py primary.fna "Primary Assembly" > stats\_primary.txt' and 'python assembly\_stats.py alternate.fna "Alternate Haplotype" > stats\_alternate.txt'
3. Run 'assembly\_plots.py'

   * Generates 4-panel diagnostic figures
4. Run 'python assembly\_plots.py primary.fna "Primary"' and python assembly\_plots.py alternate.fna "Alternate"'
5. Run 'contamination\_check.py'

   * Screens for contaminant contigs
   * Outputs 'contamination\_report.txt'
6. Run 'extract\_metadata.py'

   * Parses NCBI's official assembly metadata
   * Outputs 'ncbi\_official\_stats.txt'
7. Run 'python extract\_metadata.py > ncbi\_official\_stats.txt'
8. Run 'generate\_report.py'

   * Computes assembly statistics for primary + alternate FASTA
   * Runs a contamination screen
   * Parses NCBI official metadata (JSONL)
   * Generates 4 diagnostic plots per assembly
   * Produces a complete markdown report with all tables filled in





**Final Results**



The Leptodirus hochenwartii genome assembly (GCA\_947310635.1) was successfully analyzed using Python covering statistics, visualization, contamination screening, and metadata extraction.



Primary assembly: 492,378,077 bp across 75 sequences, with a scaffold N50 of 37.27 Mb and L50 of 6 was found and a GC content was 31.79% with only 0.01% N bases.



Alternate haplotype: 434,099,467 bp across 1,130 sequences (scaffold N50 = 1.23 Mb), representing heterozygous regions.



All computed values matched NCBI's official metadata within 0.004%, confirming data integrity.



Zero suspicious contigs flagged during contamination checks. The unimodal GC distribution and clean length-vs-GC scatter confirm the assembly is uncontaminated.



The assembly exceeds every Earth Biogenome Project benchmark, including contig N50 (3.37 Mb), scaffold N50 (37.27 Mb), QV (65.3), k-mer completeness (96.74%), and BUSCO (97.4%).







