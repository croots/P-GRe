import argparse
import os
import subprocess
import sys
import shutil
from datetime import datetime
from pathlib import Path


def run(cmd, check=True):
    result = subprocess.run(cmd, shell=True, check=check)
    return result


def main():
    print("P-GRe pipeline v1.0\n\nP-GRe is a pipeline dedicated to the automatic detection and annotation of pseudogenes. If you are using P-GRe, please cite:\n")
    print("Cabanac et al. P-GRe: An efficient pipeline for pseudogenes annotation. 2026, Genomics, Volume 118, Issue 2.")

    parser = argparse.ArgumentParser(
        prog="pgre",
        usage="pgre -f fasta -g gff -p proteins [-o output directory] [-u other proteins] [-A] [-Q] [-t threads]"
    )
    parser.add_argument("-f", dest="genome", required=True, metavar="fasta",
                        help="path to the genome to annotate, in FASTA format")
    parser.add_argument("-g", dest="gff", required=True, metavar="gff",
                        help="path to the genome annotation, in GFF format")
    parser.add_argument("-p", dest="proteome", required=True, metavar="proteins",
                        help="path to the organism protein sequences, in FASTA format")
    parser.add_argument("-o", dest="outdir", default=None, metavar="output directory",
                        help="output directory. If none is given results will be written in the directory containing the genome to annotate [fasta]")
    parser.add_argument("-u", dest="proteome2", default=None, metavar="other proteins",
                        help="path to protein sequences from other organisms, used for unitary pseudogene predictions")
    parser.add_argument("-A", dest="check_polya", action="store_true", default=False,
                        help="if activated, P-GRe will look for poly(A) tail near uncategorized pseudogene to reclassify them as retropseudogenes")
    parser.add_argument("-Q", dest="hq_filter", action="store_true", default=False,
                        help="if activated, P-GRe will filter out its predictions if they have no homology with the protein sequences")
    parser.add_argument("-t", dest="thread", type=int, default=1, metavar="threads",
                        help="number of threads to use for multi-threading [1]")

    args = parser.parse_args()

    genome = str(Path(args.genome).resolve())
    if not os.path.isfile(genome):
        print(f"ERROR with -f option: {genome} file does not exist.")
        sys.exit(1)

    gff = str(Path(args.gff).resolve())
    if not os.path.isfile(gff):
        print(f"ERROR with -g option: {gff} file does not exist.")
        sys.exit(1)

    proteome = str(Path(args.proteome).resolve())
    if not os.path.isfile(proteome):
        print(f"ERROR with -p option: {proteome} file does not exist.")
        sys.exit(1)

    proteome2 = None
    if args.proteome2 is not None:
        proteome2 = str(Path(args.proteome2).resolve())
        if not os.path.isfile(proteome2):
            print(f"ERROR with -u option: {proteome2} file does not exist.")
            sys.exit(1)

    if args.outdir is None:
        timestamp = datetime.now().strftime("%d_%m_%Y_%Hh_%Mm_%Ss")
        outdir = os.path.join(os.path.dirname(os.getcwd()), f"PGRe_{timestamp}")
        os.makedirs(outdir)
    else:
        outdir = args.outdir
        os.makedirs(outdir, exist_ok=True)

    script_path = os.path.dirname(os.path.realpath(__file__))

    command_parts = ["pgre", "-f", genome, "-g", gff, "-p", proteome]
    if proteome2:
        command_parts += ["-u", proteome2]
    if args.outdir:
        command_parts += ["-o", outdir]
    if args.check_polya:
        command_parts.append("-A")
    if args.hq_filter:
        command_parts.append("-Q")
    command_parts += ["-t", str(args.thread)]
    print(f"Running P-GRe. Command used: {' '.join(command_parts)}")

    ##########################################################################
    # WORKING FILES
    ##########################################################################

    print("\nCreating working files...")

    mrna_bed = os.path.join(outdir, "mrna_coordinates.bed")
    masked_genome = os.path.join(outdir, "masked_genome.fasta")
    in_organism_prot_id = os.path.join(outdir, "in_organism_prot.id")

    awk_cmd = (
        f"awk 'BEGIN{{OFS=\"\\t\"}} /\\tmRNA\\t/{{if ($4<$5){{print $1,$4,$5}} else {{print $1,$5,$4}}}}' "
        f"{gff} > {mrna_bed}"
    )
    run(awk_cmd)

    run(f"bedtools maskfasta -fi {genome} -bed {mrna_bed} -fo {masked_genome}")

    run(f"awk '/^>/{{print substr($1,2,length($1))}}' {proteome} > {in_organism_prot_id}")

    proteome_full = proteome
    all_prot = os.path.join(outdir, "all_prot.fa")
    if proteome2:
        run(f"cat {proteome} {proteome2} > {all_prot}")
        proteome_full = all_prot

    print("Done.")

    ##########################################################################
    # MAIN
    ##########################################################################

    miniprot_res = os.path.join(outdir, "miniprot_res.gff")
    print("\nRunning miniprot... This may take a few minutes.")
    run(
        f"miniprot -L 15 -B 5 -J 25 -F 23 -I -t {args.thread} -p 0.6 --outs 0.7 --outc 0.1 --gff --aln "
        f"{masked_genome} {proteome_full} > {miniprot_res}"
    )

    print("Done.\n\nFiltering results... This may take a few minutes.")
    run(f"python3 {script_path}/overlap_filter.py {miniprot_res} {mrna_bed} {outdir} {in_organism_prot_id}")

    miniprot_filtered = os.path.join(outdir, "miniprot_res_filtered.gff")
    check_unitarity = os.path.join(outdir, "check_unitarity.list")
    run(
        f"sed -rn '/pseudogene/s/.+Target=([^ ]+).+/\\1/gp' {miniprot_filtered} | "
        f"grep -Fxv -f {in_organism_prot_id} - > {check_unitarity}",
        check=False
    )

    diamond_res = os.path.join(outdir, "diamond.res")
    run(f"touch {diamond_res}")

    unitary_check_fasta = os.path.join(outdir, "unitary_check.fasta")
    diamond_db = os.path.join(outdir, "diamond.db")
    if proteome2:
        print("Done.\n\nLooking for unitary pseudogenes with DIAMOND... This may take a few minutes.")
        run(f"python3 {script_path}/extract_unitary.py {check_unitarity} {proteome2} > {unitary_check_fasta}")
        run(f"diamond makedb --in {proteome} -d {diamond_db} -p {args.thread}")
        run(f"diamond blastp -p {args.thread} -k 1 -f 6 qseqid sseqid -d {diamond_db} --query {unitary_check_fasta} --out {diamond_res}")

    in_op_full_id = os.path.join(outdir, "in_op_full.id")
    run(
        f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $1,$1}}' {in_organism_prot_id} | "
        f"cat - {diamond_res} > {in_op_full_id}"
    )

    exp_struc = os.path.join(outdir, "exp_struc.tsv")
    print("Done\n\nChecking pseudogene's parent structures to determine their types...")
    run(f"python3 {script_path}/compute_expected_structure.py {gff} > {exp_struc}")

    ##########################################################################
    # OUTPUTS
    ##########################################################################

    pgre_unsorted = os.path.join(outdir, "PGRe.unsorted.res")
    pgre_gff = os.path.join(outdir, "PGRe.gff")
    pseudogene_protein = os.path.join(outdir, "pseudogene_protein.fasta")

    print("Done.\n\nInferring pseudogene types and writing output...")
    run(f"python3 {script_path}/write_output.py {in_op_full_id} {exp_struc} {miniprot_filtered} {outdir} > {pgre_unsorted}")
    run(f"sort --version-sort -k1,1 -k4,4 {pgre_unsorted} > {pgre_gff}")
    run(f"python3 {script_path}/get_seq.py {pgre_gff} {miniprot_res} > {pseudogene_protein}")

    ##########################################################################
    # OPTIONAL OUTPUTS
    ##########################################################################

    check_quality_db = os.path.join(outdir, "check_quality.db")
    check_quality_res = os.path.join(outdir, "check_quality.res")

    if args.hq_filter:
        print("Done.\n\nKeeping high-quality predictions only...")
        run(f"diamond makedb --in {proteome_full} -d {check_quality_db} -p {args.thread}")
        run(f"diamond blastp -p {args.thread} -k 1 -f 6 qseqid -d {check_quality_db} --query {pseudogene_protein} --out {check_quality_res}")
        run(f"python3 {script_path}/filter_low_qual.py {check_quality_res} {pgre_gff} {outdir} {pseudogene_protein}")
        os.rename(pgre_gff, os.path.join(outdir, "PGRe_w_low_qual.gff"))
        os.rename(pseudogene_protein, os.path.join(outdir, "pseudogene_protein_w_low_qual.fasta"))
        os.rename(os.path.join(outdir, "PGRe_filtered.gff"), pgre_gff)
        os.rename(os.path.join(outdir, "pseudogene_protein_filtered.fasta"), pseudogene_protein)

    scaffolds_length = os.path.join(outdir, "scaffolds_length.tsv")
    three_prime_bed = os.path.join(outdir, "three_prime.bed")
    three_prime_fasta = os.path.join(outdir, "three_prime.fasta")
    found_polya = os.path.join(outdir, "found_polyA_tail.tsv")
    pgre_not_reclassified = os.path.join(outdir, "PGRe_not_reclassified.gff")

    if args.check_polya:
        os.rename(pgre_gff, pgre_not_reclassified)
        print("Done.\n\nChecking for poly(A)-tails to classify uncategorized pseudogenes...")
        run(
            f"awk '/^>/{{if(l)print n\"\\t\"l;n=substr($0,2);l=0;next}}{{l+=length}}END{{print n\"\\t\"l}}' "
            f"{genome} > {scaffolds_length}"
        )
        run(f"python3 {script_path}/rePolarizeGFF.py {pgre_not_reclassified} {scaffolds_length} > {three_prime_bed}")
        run(f"bedtools getfasta -fi {genome} -fo {three_prime_fasta} -bed {three_prime_bed} -s -name")
        run(f"python3 {script_path}/getPolyA.py {three_prime_fasta} > {found_polya}")

        polya_count = int(subprocess.check_output(f"wc -l {found_polya} | awk '{{print $1}}'", shell=True).decode().strip())
        print(f"Found {polya_count} pseudogene associated poly(A)-tails.")

        run(f"python3 {script_path}/reclassifyGFF.py {pgre_not_reclassified} {found_polya} > {pgre_gff}")

    ##########################################################################
    # CLEANUP
    ##########################################################################

    print("Done.\n\nCleaning working directory...")
    tmp_dir = os.path.join(outdir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    for item in os.listdir(outdir):
        item_path = os.path.join(outdir, item)
        if item_path != tmp_dir:
            try:
                shutil.move(item_path, tmp_dir)
            except Exception:
                pass

    for fname in ["PGRe.gff", "pseudogene_protein.fasta"]:
        src = os.path.join(tmp_dir, fname)
        dst = os.path.join(outdir, fname)
        if os.path.isfile(src):
            shutil.move(src, dst)

    print("Done.\n\nThank you for using P-GRe.")


if __name__ == "__main__":
    main()
