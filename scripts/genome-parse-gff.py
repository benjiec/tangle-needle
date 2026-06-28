import sys
import gffutils
from Bio import SeqIO
from defaults import Defaults


def extract_intron_sequences(gff_path, fna_file, target_protein_id):
    """
    Parses a GFF3 and FNA file to isolate and slice intron sequences 
    flanked by the exons of a targeted protein product ID.
    """
    # 1. Initialize an on-the-fly relational database from the GFF text file
    print(f"[*] Parsing and indexing GFF structure...", file=sys.stderr)
    db = gffutils.create_db(
        gff_path, 
        dbfn=':memory:',  # Avoids disk overhead for standalone execution
        force=True, 
        keep_order=True, 
        merge_strategy='merge', 
        sort_attribute_values=True
    )
    
    # 2. Load the genomic nucleotide assembly into memory
    print(f"[*] Loading FNA genome reference...", file=sys.stderr)
    genome_dict = SeqIO.to_dict(SeqIO.parse(fna_file, "fasta"))
    
    # 3. Discover the parent mRNA identifier tied to the requested protein_id
    parent_transcript_id = None
    for cds in db.features_of_type('CDS'):
        if 'protein_id' in cds.attributes:
            if any(target_protein_id in pid for pid in cds.attributes['protein_id']):
                parent_transcript_id = cds.attributes['Parent'][0]
                break
                
    if not parent_transcript_id:
        raise ValueError(f"Error: Target protein ID '{target_protein_id}' not resolved within GFF CDS features.")
        
    print(f"[+] Found parent transcript lineage: {parent_transcript_id}", file=sys.stderr)
    
    # 4. Extract and sort child exons matching that specific mRNA parent
    exons = list(db.children(parent_transcript_id, featuretype='exon', order_by='start'))
    
    if len(exons) <= 1:
        print(f"[-] Warning: Transcript {parent_transcript_id} contains {len(exons)} exons. No introns exist.", file=sys.stderr)
        return []
        
    strand = exons[0].strand
    chrom = exons[0].seqid
    
    if chrom not in genome_dict:
        raise KeyError(f"Error: Sequence ID '{chrom}' from GFF not found in FNA file records.")
        
    chromosome_sequence = genome_dict[chrom].seq
    intron_records = []
    
    # 5. Inter-feature coordinate logic mapping to isolate coordinate spans
    # GFF coordinate space: 1-based, fully inclusive [start, end]
    for i in range(len(exons) - 1):
        # Intron boundary begins 1 base past current exon end, terminates 1 base before next exon start
        intron_start = exons[i].end + 1
        intron_end = exons[i+1].start - 1
        
        if intron_start > intron_end:
            continue  # Filters out overlapping features or annotation discrepancies
            
        # Transform 1-based inclusive GFF bounds to Python 0-based half-open slice notation
        intron_seq = chromosome_sequence[intron_start - 1 : intron_end]
        
        # Enforce biological strandedness reverse-complementation
        if strand == '-':
            intron_seq = intron_seq.reverse_complement()
            
        intron_records.append({
            "genomic_coord": f"{chrom}:{intron_start}-{intron_end}",
            "strand": strand,
            "sequence": str(intron_seq)
        })
        
    # 6. Apply biological numbering relative to transcription direction
    if strand == '-':
        for idx, record in enumerate(reversed(intron_records)):
            record["biological_index"] = idx + 1
    else:
        for idx, record in enumerate(intron_records):
            record["biological_index"] = idx + 1
            
    return sorted(intron_records, key=lambda x: x["biological_index"])


if __name__ == "__main__":

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("genome_accession")
    ap.add_argument("protein_id")
    ap.add_argument("--exon", type=int)
    ap.add_argument("--intron", type=int)
    args = ap.parse_args()

    gff_file = Defaults.ncbi_genome_gff(args.genome_accession)
    fna_file = Defaults.ncbi_genome_fna(args.genome_accession)
    print(gff_file, fna_file)

    try:
        introns = extract_intron_sequences(gff_file, fna_file, args.protein_id)
        for intron in introns:
            print(f">Intron_{intron['biological_index']} | {intron['genomic_coord']} | Strand: {intron['strand']}\n{intron['sequence']}")
    except Exception as error:
        print(f"[!] Execution Failure: {error}", file=sys.stderr)
        sys.exit(1)
