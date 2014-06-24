__author__ = 'janos'

import csv
import json
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import StringIO
import os


def main(genes_fkpm, fasta_file, transcripts_gtf_file, load_sequence_dict=True, load_gene_dict=True,
         length_of_sequence_to_match=72):

    gene_dict = {}
    scaffold_dict = {}

    if load_gene_dict and os.path.exists("gene_dict.json"): #Load a preexisiting "genes_dict.json"
        with open("gene_dict.json", "r") as fj:
            gene_dict = json.load(fj)

    else:

        with open(genes_fkpm, 'r') as f:
            genes = csv.DictReader(f, delimiter="\t")
            i = 0
            for gene in genes:

                gene_id = gene["gene_id"]
                locus = gene["locus"]

                scaffold, raw_position = tuple(locus.split(":"))
                start_position, end_position = tuple(raw_position.split("-"))

                start_position, end_position = int(start_position), int(end_position)

                if scaffold in scaffold_dict:
                    scaffold_dict[scaffold] += [(gene_id, start_position, end_position)]
                else:
                    scaffold_dict[scaffold] = [(gene_id, start_position, end_position)]

                fpkm_numeric = float(gene["FPKM"])
                gene["FPKM_numeric"] = fpkm_numeric

                gene["locus_parsed"] = {"scaffold": scaffold, "start_position": start_position,
                                        "end_position": end_position, "length": end_position - start_position}

                gene_dict[gene_id] = gene

                i += 1

            for gene_id in gene_dict.keys():
                if gene_dict[gene_id]["FKPM_numeric"] == 0.0:
                    gene_dict.pop(gene_id)

    if load_sequence_dict:
        with open("gene_sequence_dict.json", "r") as fj:
            gene_sequence_dict = json.load(fj)
    else:
        gene_sequence_dict = find_sequences_in_genome(fasta_file, scaffold_dict)

        with open("gene_sequence_dict.json", "w") as fwj:
            json.dump(gene_sequence_dict, fwj, indent=4, sort_keys=True)

    exon_boundaries_dict = find_exon_boundaries(transcripts_gtf_file, gene_dict)
    for gene_id in exon_boundaries_dict.keys():
        gene_dict[gene_id]["exon_boundaries"] = exon_boundaries_dict[gene_id]

    for gene_id in gene_sequence_dict.keys():
        sequence = gene_sequence_dict[gene_id]
        if gene_id in gene_dict:
            gene_dict[gene_id]["locus_parsed"]["sequence"] = sequence

    with open("gene_dict.json", "w") as fjw:
                json.dump(gene_dict, fjw, indent=4, sort_keys=True)

    for gene_id in gene_dict.keys():
        sequence = gene_sequence_dict[gene_id]
        first_start_position, first_end_position = tuple(gene_dict[gene_id]["exon_boundaries"][0])

        if first_end_position < length_of_sequence_to_match:
            len_seq = first_end_position
        else:
            len_seq = length_of_sequence_to_match

        if "blast_matches" not in gene_dict[gene_id]["locus_parsed"]:
            match_list = blast_sequence_get_names_of_top_matches(sequence[0:len_seq])

            gene_dict[gene_id]["locus_parsed"]["blast_matches"] = match_list

            with open("gene_dict.json", "w") as fjw:
                json.dump(gene_dict, fjw, indent=4, sort_keys=True)

        with open("gene_dict.json", "w") as fjw:
            json.dump(gene_dict, fjw, indent=4, sort_keys=True)


def find_sequences_in_genome(fasta_file_path, scaffold_dict):
    from Bio.SeqIO import FastaIO
    sequence_dict = {}
    with open(fasta_file_path, "r") as fg:
        genome = FastaIO.FastaIterator(fg)

        j = 0
        for sequence in genome:
            if j % 100000 == 0:
                print("Read %s scaffolds" % j)

            sequence_id = sequence.id

            if sequence_id in scaffold_dict:
                gene_matches = scaffold_dict[sequence_id]
                sequence_as_a_string = sequence.seq
                for gene_match in gene_matches:
                    sequence_dict[gene_match[0]] = str(sequence_as_a_string[gene_match[1]: gene_match[2] + 1])
            j += 1

    return sequence_dict


def find_exon_boundaries(gtf_file_name, gene_dict):

    with open(gtf_file_name, "r") as fg:
        exon_boundaries_dict = {}
        gtf_file = csv.reader(fg, delimiter="\t")
        for gtf in gtf_file:

            start_position = int(gtf[3])
            end_position = int(gtf[4])

            type_of_row = gtf[2]

            gene_id = gtf[-1].split("; ")[0].split()[1][1:-1]

            if type_of_row == "exon":
                if gene_id in gene_dict:
                    if gene_id in exon_boundaries_dict:
                        exon_boundaries_dict[gene_id] += [(start_position - 1, end_position - 1)]
                    else:
                        exon_boundaries_dict[gene_id] = [(start_position - 1, end_position - 1)]

        for gene_id in exon_boundaries_dict.keys():
            exon_boundaries_dict[gene_id].sort(key=lambda x: x[0])

    return exon_boundaries_dict


def blast_sequence_get_names_of_top_matches(sequence, top_n_matches=5):
    print("Blasting sequence: '%s'" % sequence)
    blast_results = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_results_read = blast_results.read()
    blast_results_fp = StringIO.StringIO(blast_results_read)

    blast_results_obj = NCBIXML.parse(blast_results_fp)

    alignment_list = []

    for searches in blast_results_obj:
        i = 0
        for alignment in searches.alignments:
            if i == top_n_matches:
                break

            alignment_dict = {"title": alignment.title, "expect": alignment.hsps[0].expect}
            alignment_list += [alignment_dict]

            i += 1

    return alignment_list


if __name__ == "__main__":
    main("/Users/janos/rna/genes.fpkm_tracking",
         "/Users/janos/rna/rna_zip/data/sater/to_align/LAEVIS_7.1.repeatMasked.fa",
         "/Users/janos/rna/transcripts.gtf"
         )
