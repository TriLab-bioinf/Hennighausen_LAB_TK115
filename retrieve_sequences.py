from Bio import SeqIO
import sys

def retrieve_protein_cds_sequences(gbk_file, output_file):
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(gbk_file, 'genbank'):
            for feature in record.features:
                if feature.type == 'CDS' and 'translation' in feature.qualifiers:
                    gene_name = feature.qualifiers.get('gene', [''])[0]
                    protein_sequence = feature.qualifiers['translation'][0]
                    out_f.write(f">{gene_name}\n{protein_sequence}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python retrieve_sequences.py gbk_file output_file")
        sys.exit(1)

    gbk_file = sys.argv[1]
    output_file = sys.argv[2]

    retrieve_protein_cds_sequences(gbk_file, output_file)
