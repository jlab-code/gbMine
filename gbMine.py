from dataclasses import dataclass
from intervaltree import Interval, IntervalTree
import argparse
from scipy import stats
import multiprocessing
import FilterTable
import re
import pandas as pd
import statsmodels.stats.multitest as smt
import os
import traceback


# Definition of Gene using dataclass decorator
@dataclass
class Gene:
    id: str
    chr: str
    strand: str
    start: int
    end: int
    exons: IntervalTree
    exon_num: int
    length: int
    tot_cg: int = 0
    tot_chh: int = 0
    tot_chg: int = 0
    meth_cg: int = 0
    meth_chg: int = 0
    meth_chh: int = 0


@dataclass()
class Exon:
    start: int
    end: int


# Adds gene with start/end to Interval tree
def add_gene_to_tree(gene: Gene, tree: IntervalTree) -> IntervalTree:
    interval = Interval(gene.start, gene.end, gene)
    tree.add(interval)
    return tree


# Write header of file
def write_header(file):
    file.write(
        "Chromosome\tStrand\tGeneID\tCG p-Value\tCorrected_CG\tCHG p-Value\tCorrected_CHG\tCHH p-Value\tCorrected_CHH\t"
        "non-CG p-Value\tCorrected_non_CG\tGene length\tnumber of exons\tnumber of cytosines\n")


# Read gff-File and save genes with characteristics in dictionary:
# key = chromosome + strand, value = Interval Tree of included genes
def read_gff(gff_path: str, exons: bool, contigs: list):
    including_introns = not exons
    gene_trees = {}

    try:
        # open gff file
        with (open(gff_path, 'r')) as gff_file:
            chr = ""
            geneid = ""
            strand = ""
            pos = []
            exons = IntervalTree()
            exon_num = 0
            for line in gff_file:
                if not line.startswith("#"):
                    if not line.rstrip():
                        continue
                    data = line.split("\t")
                    if data[2] == "gene":
                        if geneid != "":
                            start = min(pos)
                            end = max(pos)
                            if chr in contigs or not contigs:
                                # create gene object
                                g = Gene(geneid, chr, strand, start, end, IntervalTree(exons), exon_num, end - start)
                                key = chr + strand

                                # add gene to Intervaltree
                                if key not in gene_trees:
                                    gene_trees[key] = IntervalTree()

                                add_gene_to_tree(g, gene_trees.get(key))

                        chr = data[0]
                        strand = data[6]
                        pos.clear()
                        exons.clear()
                        exon_num = 0
                        pos.append(int(data[3]))
                        pos.append(int(data[4]))
                        attr = data[8].split(';')
                        for attribute in attr:
                            if "gene_id" in attribute or "ID" in attribute:
                                geneid = attribute.replace("gene_id=", "").replace("ID=", "")
                                break

                    # if only exonic regions are relevant, add exon to Intervaltree
                    elif not including_introns and data[2] == "exon":
                        # if "exon_id=" + geneid in data[8] or "Parent=" + geneid in data[8]:
                        exon = Exon(int(data[3]), int(data[4]) + 1)

                        exons.addi(exon.start, exon.end, )
                        exon_num += 1

                    elif including_introns and data[2] == "exon":
                        exon_num += 1

            start = min(pos)
            end = max(pos)
            if chr in contigs or not contigs:
                g = Gene(geneid, chr, strand, start, end, exons, exon_num, end - start)
                key = chr + strand
                if key not in gene_trees:
                    gene_trees[key] = IntervalTree()

                add_gene_to_tree(g, gene_trees.get(key))
            return gene_trees

    # if gff-file path cannot be found
    except FileNotFoundError:
        print(f"The file '{gff_path}' does not exist.")


# increase total counts and gene count for given gene and given status. Methylation counts are increased
# too if cytosine is methylated
def increase_gene_count(context: str, gene: Gene, methylated: str, cg_total: int, cg_methyl: int, chg_total: int,
                        chg_methyl: int, chh_total: int, chh_methyl: int):
    # Depending on the file, status does not necessarily be CG, CHH and CHG, but could be more specific, e.g. CTG
    # CG & CGH status
    if context.startswith("CG"):
        gene.tot_cg += 1
        cg_total += 1
        if methylated == "M" or methylated == "1":
            gene.meth_cg += 1
            cg_methyl += 1
    # CHH status
    elif re.match(r'^C[^G]*$', context):
        chh_total += 1
        gene.tot_chh += 1
        if methylated == "M" or methylated == "1":
            gene.meth_chh += 1
            chh_methyl += 1
    # CHG status
    elif re.match(r'^C[^G]*G$', context):
        gene.tot_chg += 1
        chg_total += 1
        if methylated == "M" or methylated == "1":
            gene.meth_chg += 1
            chg_methyl += 1

    return cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl


def increase_total_count(context: str, methylated: str, cg_total: int, cg_methyl: int, chg_total: int,
                         chg_methyl: int, chh_total: int, chh_methyl: int):

    if context.startswith("CG"):
        cg_total += 1
        if methylated == "M" or methylated == "1":
            cg_methyl += 1
    # CHH status
    elif re.match(r'^C[^G]*$', context):
        chh_total += 1
        if methylated == "M" or methylated == "1":
            chh_methyl += 1
    # CHG status
    elif re.match(r'^C[^G]*G$', context):
        chg_total += 1
        if methylated == "M" or methylated == "1":
            chg_methyl += 1

    return cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl


# calculates p-value for every gene, for every methylation type
def look_through_tree(intervaltree: IntervalTree, cg_prop: float, chg_prop: float, chh_prop: float, non_cg_prop: float):
    results = []
    for interval in intervaltree:
        gene = interval.data
        result = (gene.chr, gene.strand, gene.start, gene.end, gene.id,
                  str(calc_p_value(gene.meth_cg, gene.tot_cg, cg_prop)),
                  str(calc_p_value(gene.meth_chg, gene.tot_chg, chg_prop)),
                  str(calc_p_value(gene.meth_chh, gene.tot_chh, chh_prop)),
                  str(calc_p_value((gene.meth_chh + gene.meth_chg), (gene.tot_chg + gene.tot_chh), non_cg_prop)),
                  (gene.end - gene.start + 1), gene.exon_num,
                  str(gene.tot_chg + gene.tot_chh + gene.tot_cg))
        results.append(result)

    return results


# calculates p-value through binomial test
def calc_p_value(methyl: int, tot_num: int, prop: float) -> float:
    if methyl == 0 or tot_num == 0:
        return 1.0
    sum = stats.binomtest(methyl, tot_num, prop, "greater")
    return sum.pvalue


# read methylome-file and count methylated cytosines for every gene
def map_methyls_on_genes(methylome: str, exons: bool, genomic_wide: bool, gene_dict,
                         mindepth: int, out: str, overlap: bool, cutoff: float, chr_dict, filter: bool):
    included_introns = not exons
    cg_total = 0
    cg_methyl = 0
    chg_methyl = 0
    chg_total = 0
    chh_total = 0
    chh_methyl = 0
    gene_pos = None
    gene_neg = None
    not_found_in_gff = []

    #map cytosines
    def mapper(methyl_file, overlap, chr_dict):
        nonlocal cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl, gene_pos, gene_neg, not_found_in_gff
        header = methyl_file.readline()
        for line in methyl_file:
            data = line.split('\t')

            if line == header:
                continue

            chr = data[0]
            pos = int(data[1])
            strand = data[2]
            context = data[3]
            depth = int(data[5])

            if len(data) >= 9:
                methylated = data[7].strip()

            else:
                methylated = data[6].strip()

            if chr in chr_dict:
                chr = chr_dict.get(chr)

            if depth >= mindepth:
                if not overlap and gene_pos is not None and gene_pos.chr == chr and strand == "+" \
                        and gene_pos.start <= pos <= gene_pos.end:
                    if included_introns:
                        cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl = \
                            increase_gene_count(context, gene_pos, methylated, cg_total, cg_methyl,
                                                chg_total, chg_methyl, chh_total, chh_methyl)
                    else:
                        found_exons = gene_pos.exons[pos]
                        if len(found_exons) > 0:
                            cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl = increase_gene_count(
                                context, gene_pos, methylated, cg_total,
                                cg_methyl, chg_total, chg_methyl, chh_total,
                                chh_methyl)
                        else:
                            if genomic_wide:
                                cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl = \
                                    increase_total_count(context, methylated, cg_total, cg_methyl, chg_total,
                                                         chg_methyl, chh_total, chh_methyl)

                elif not overlap and gene_neg is not None and gene_neg.chr == chr and strand == "-" \
                        and gene_neg.start <= pos <= gene_neg.end:
                    if included_introns:
                        cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl = \
                            increase_gene_count(context, gene_neg, methylated, cg_total, cg_methyl, chg_total,
                                                chg_methyl, chh_total, chh_methyl)

                    else:
                        found_exons = gene_neg.exons[pos]
                        if len(found_exons) > 0:
                            cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl = increase_gene_count(
                                context, gene_neg, methylated, cg_total,
                                cg_methyl, chg_total, chg_methyl, chh_total,
                                chh_methyl)
                        else:
                            if genomic_wide:
                                cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl = (
                                    increase_total_count(context, methylated, cg_total, cg_methyl,
                                                         chg_total, chg_methyl, chh_total, chh_methyl))

                else:
                    key = chr + strand
                    chromosome_tree = gene_dict.get(key)
                    try:
                        found_genes = chromosome_tree[pos]
                    except:
                        if chr not in not_found_in_gff:
                            not_found_in_gff.append(chr)
                        found_genes = []
                    if len(found_genes) > 0:
                        for interval in found_genes:
                            gene = interval.data
                            if strand == "+":
                                gene_pos = gene
                            else:
                                gene_neg = gene

                            if included_introns:
                                cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl = increase_gene_count(
                                    context, gene, methylated, cg_total, cg_methyl, chg_total, chg_methyl, chh_total,
                                    chh_methyl)
                            else:
                                found_exons = gene.exons[pos]
                                if len(found_exons) > 0:
                                    cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl = \
                                        increase_gene_count(context, gene, methylated, cg_total, cg_methyl, chg_total,
                                                            chg_methyl, chh_total, chh_methyl)
                                else:
                                    if genomic_wide:
                                        cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl = \
                                            increase_total_count(context, methylated, cg_total, cg_methyl, chg_total,
                                                                 chg_methyl, chh_total, chh_methyl)

                    else:
                        if genomic_wide:
                            cg_total, cg_methyl, chg_total, chg_methyl, chh_total, chh_methyl = (
                                increase_total_count(context, methylated, cg_total, cg_methyl, chg_total, chg_methyl,
                                                     chh_total, chh_methyl))

    try:
        with (open(methylome, 'r')) as methyl_file:
            mapper(methyl_file, overlap, chr_dict)
            print("Methylome-file processed, calculating p-values...")

    except Exception as e:
        print(str(e))
        return None

    cg_prop = 0
    chg_prop = 0
    chh_prop = 0
    non_gc_prop = 0
    # calculate proportion of methylated cytosine residues per status of whole genome
    if cg_total != 0:
        cg_prop = cg_methyl / cg_total
    if chg_total != 0:
        chg_prop = chg_methyl / chg_total
    if chh_total != 0:
        chh_prop = chh_methyl / chh_total
    if chg_total != 0 and chh_total != 0:
        non_gc_prop = (chh_methyl + chg_methyl)/(chh_total + chg_total)
    print(f' CG background distribution: {cg_prop}')
    print(f' CHG background distribution: {chg_prop}')
    print(f' CHH background distribution: {chh_prop}')
    print(f' non-CG background distribution: {non_gc_prop}')
    results = []
    # loop parallel over genes and calculate p-values
    for key, value in gene_dict.items():
        print(f' Processing Chromosome: {key}...')
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        results += pool.apply(look_through_tree, args=(value, cg_prop, chg_prop, chh_prop, non_gc_prop))

        pool.close()
        pool.join()
    try:
        methyl_filename = os.path.basename(methylome)
        data = pd.DataFrame(results)
        column_names = ["Chromosome", "Strand", "Start", "End", "GeneID", "CG p-Value", "CHG p-Value", "CHH p-Value",
                        "non-CG p-Value", "Gene length", "number of exons", "number of cytosines"]
        data.columns = column_names
        data["CG p-Value"] = data["CG p-Value"].astype(float)
        data["CHG p-Value"] = data["CHG p-Value"].astype(float)
        data["CHH p-Value"] = data["CHH p-Value"].astype(float)
        data["non-CG p-Value"] = data["non-CG p-Value"].astype(float)

        # correct p-values
        data['Corrected_CG'] = smt.multipletests(data["CG p-Value"], method='fdr_bh')[1]
        data['Corrected_CHG'] = smt.multipletests(data["CHG p-Value"], method='fdr_bh')[1]
        data['Corrected_CHH'] = smt.multipletests(data["CHH p-Value"], method='fdr_bh')[1]
        data['Corrected_non_CG'] = smt.multipletests(data["non-CG p-Value"], method='fdr_bh')[1]

        new_order = ["Chromosome", "Strand", "Start", "End", "GeneID", "CG p-Value", "Corrected_CG", "CHG p-Value",
                     "Corrected_CHG", "CHH p-Value", "Corrected_CHH", "non-CG p-Value", "Corrected_non_CG",
                     "Gene length", "number of exons", "number of cytosines"]

        data = data[new_order]
        all_path = out + "all_genes_" + methyl_filename + ".txt"
        data.to_csv(all_path, index=False, sep='\t')

        if filter:
            FilterTable.filter(all_path, out + "filtered_" + methyl_filename + ".txt")

        # gbM filter
        if cutoff is not None:
            filtered_data = data[(data['Corrected_CG'] <= cutoff) & (data['Corrected_CHG'] >= cutoff) &
                                 (data['Corrected_CHH'] >= cutoff)]
            filtered_data.to_csv(out + "gbm_genes_" + methyl_filename + ".txt", index=False, sep='\t')
            return filtered_data
        return None

    except FileNotFoundError:
        print("Could not find path to ", FileNotFoundError)

# Read contig file
def read_non_contigs(path: str):
    contigs = []
    with (open(path, 'r')) as contig_file:
        for line in contig_file:
            contigs.append(line.rstrip('\n'))

    return contigs

# read Conversion file
def read_conversion(path: str):
    chr_dict = {}
    with (open(path, 'r')) as conv_file:
        for line in conv_file:
            content = line.strip().split('\t')
            if len(content) == 2:
                key, value = content
                chr_dict[key] = value
    return chr_dict


if __name__ == '__main__':

    # parse arguments
    parser = argparse.ArgumentParser(description="Command-line parser of Gene Methylation Classification", )
    parser.add_argument("--gff", help="Required Path to gff3-File", required=True)
    parser.add_argument("--methyl",
                        help="Required Path to methylome-File: can be .tsv, .allc (methylpy output format),"
                             " .CX_report.txt.gz (optional genome-wide cytosine report output of Bimarck)",
                        required=True)
    parser.add_argument("--out", help="Path to output-directory", required=True)
    parser.add_argument("--mincov", help="Minimum Coverage of analyzed methylation sites (default = 0)",
                        required=False, default=0)
    parser.add_argument("--exons", help="If set only exonic parts of the gene are analyzed", action="store_true")
    parser.add_argument("--genomicBackground", help="If set, the background distribution is calculated from the entire"
                                                    " genome; otherwise, it is calculated only from gene regions",
                        action="store_true", required=False)
    parser.add_argument("--overlap", help="Genome contains overlapping genes", action="store_true", default=False,
                        required=False)
    parser.add_argument("--gbm", help="Filters automatically gbM genes by given cut-off value",
                        required=False)
    parser.add_argument("--filter", help="Filter output directly (saved separately)", action="store_true")
    parser.add_argument("--coreset", help="Determine core set of gbM genes found in all files", action="store_true")
    parser.add_argument("--draft",
                        help="If genome input is a draft assembly, an additional file with only considered IDs of "
                             "chromosomes is given", required=False)
    parser.add_argument("--convert",
                        help="If GFF3 file and ALLC-file do not have equal chromosome label, a file can be added"
                             "to convert the labels", required=False)

    args = parser.parse_args()
    if args.coreset:
        if args.gbm is None:
            print("Parameter --gbm is mandatory, when flag --coreset is set\n"
                  "Continuing with default value 0.05 for cutoff")
            args.gbm = 0.05

    out = args.out
    if not out.endswith("/"):
        out += "/"

    non_contigs = []
    conv_dict = {}
    cutoff = None

    try:
        if args.draft is not None:
            non_contigs = read_non_contigs(args.draft)

        if args.convert is not None:
            conv_dict = read_conversion(args.convert)

        if args.gbm is not None:
            cutoff = float(args.gbm)

        mincov = int(args.mincov)
        genetrees = read_gff(args.gff, args.exons, non_contigs)
        gff_filename = os.path.basename(args.gff)

        print("Reading gff-file done")
        if os.path.isfile(args.methyl):

            core_df = map_methyls_on_genes(args.methyl, args.exons, args.genomicBackground, genetrees, mincov, out,
                                           args.overlap, cutoff, conv_dict, args.filter)
            if args.coreset:
                core_df.to_csv(out + "coreset_gbm_" + gff_filename + ".txt", index=False, sep='\t')

        else:
            directory = os.fsencode(args.methyl)
            gbm_dfs = []
            for file in os.listdir(directory):
                filename = os.fsdecode(file)
                if (filename.endswith(".txt") or filename.endswith(".tsv") or filename.endswith(".csv")
                        or filename.endswith(".allc")):
                    print(f'\nProcessing {filename}...')
                    path = os.path.join(directory.decode('utf-8'), filename)
                    df_gbm = map_methyls_on_genes(path, args.exons, args.genomicBackground, genetrees,
                                                  mincov, out, args.overlap, cutoff, conv_dict, args.filter)
                    if df_gbm is not None:
                        gbm_dfs.append(df_gbm)

            if args.coreset:
                common_gbm = set(gbm_dfs[0]["GeneID"])

                for df in gbm_dfs[1:]:
                    common_gbm = common_gbm.intersection(set(df["GeneID"]))

                coreset_df = gbm_dfs[0][gbm_dfs[0]["GeneID"].isin(common_gbm)]
                coreset_df = coreset_df[["Chromosome", "Strand", "Start", "End", "GeneID", "Gene length",
                                         "number of exons"]]
                coreset_df.to_csv(out + "coreset_gbm_" + gff_filename + ".txt", index=False, sep='\t')

    except Exception as e:
        traceback.print_exc()