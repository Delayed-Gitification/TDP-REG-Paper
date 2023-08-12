import dnaio
import argparse
from rapidfuzz.fuzz import partial_ratio
from statistics import mode


def rev_c(seq):
    """
    simple function that reverse complements a given sequence
    """
    tab = str.maketrans("ACTGN", "TGACN")
    # first reverse the sequence
    seq = seq[::-1]
    # and then complement
    seq = seq.translate(tab)
    return seq


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--consensus_frac", type=float, default=0.85, help="Fraction of reads with consensus NT to be considered consensus")
    parser.add_argument("--fuzz_score", type=float, default=93, help="score, out of 100, for fuzzy matching")
    parser.add_argument("--splicing_csv", "-s", type=str, required=True, help="filename of csv with splicing information")
    parser.add_argument("--consensus_csv", type=str, required=True, help="filename of csv with consensus sequences")

    args = parser.parse_args()


    # Defined for forward strand
    end_exon_1 = "agatacaccggctggggcag".upper()
    start_exon_2 = "gtAtccggccagggcgat".upper()

    # Define expected sequences for a no_CE read
    no_CE_jun = end_exon_1 + start_exon_2
    no_CE_jun_rc = rev_c(no_CE_jun)

    # Calculate rev_c of these
    end_exon_1_rc = rev_c(end_exon_1)
    start_exon_2_rc = rev_c(start_exon_2)

    # defined on reverse strand, seqs upstream and downstream (read2-wise) of upi 
    upi_up = "GGCCTGCGGATCC"
    upi_down = "CAGAAGTTATCCTTCA"

    # Define some forward strand sequences present in intron 1
    intron1 = ["ATGCACATCACT", "AATGACACTCAGTGCC", "ATATCTACACTTTAAAA"]

    # Define some reverse strand sequences present in intron 2
    intron2 = ["GAGTAAAGATAAGTCCAG".upper(), "cCTGACAAAGGAGTAAA".upper(), "TCCAGTTACAGCCCCTGA"]

    # ce_down_rc = rev_c(downstream)
    # ce_up_rc = rev_c(upstream[12:])

    # The sequence of the CE
    ce_seq = "RCTNTCNCGNAARCTNATHAAYGGNATHCGNGAYAARCARTCNGGNAARACNATHCTNGAYTTYCTNAARTCNGAYGGNTTYGCNAAYCGNAAYTTYATGCARCTNATHCAYGAYGAYTCNCTNACNTTYAARGARGAYATHCARAARGCNCAR"

    # Check a few positions in read 1 to make sure that they are as expected
    pos_f = [1, 2, 4, 5, 7, 8, 10, 11, 13, 14]
    ce_detect = ''.join([ce_seq[i] for i in pos_f])

    ### ANALYSE SPLICING ###
    
    splicing_d = {}
    categories = ["intron1_IR", "intron2_IR", "both_IR", "no_CE", "with_CE"]
    ce_seq_d = {}

    with dnaio.open(args.r1) as file1, dnaio.open(args.r2) as file2, open(args.splicing_csv, 'w') as splicing_csv:
        splicing_csv.write("upi,umi,intron1_IR,intron2_IR,both_IR,no_CE,with_CE,CE_seq\n")

        for r1, r2 in zip(file1, file2):
            
            # Define a dictionary to store details about each read
            this_read_d = {a: False for a in categories}

            # Check that the end of exon 1 is in the sequence
            if end_exon_1 not in r1.sequence:
                continue

            # Check that start of exon 2 is in read 2
            if start_exon_2_rc not in r2.sequence:
                continue

            # Determine whether there is IR
            if sum([1 for a in intron1 if a in r1.sequence]) > 0:
                this_read_d["intron1_IR"] = True

            if sum([1 for a in intron2 if a in r2.sequence]) > 0:
                this_read_d["intron2_IR"] = True

            # Determine whether it's a "no_CE" read
            if partial_ratio(str(r1.sequence), no_CE_jun) > args.fuzz_score:
                this_read_d["no_CE"] = True

            elif partial_ratio(str(r2.sequence), no_CE_jun_rc) > args.fuzz_score:
                this_read_d["no_CE"] = True

            # Determine whether it's a "with_CE" read
            ce_start = r1.sequence.find(end_exon_1) + len(end_exon_1)

            if ''.join([r1.sequence[a + ce_start] for a in pos_f]) == ce_detect:
                if ce_start > 10:
                    this_read_d["with_CE"] = True

            # Some sanity checks 
            if this_read_d["intron1_IR"] or this_read_d["intron2_IR"]:
                if this_read_d["with_CE"] or this_read_d["no_CE"]:
                    continue  # can't have intron retention and splicing

            if this_read_d["with_CE"] and this_read_d["no_CE"]:
                continue  # can't have and not have the CE

            # Some house keeping
            if this_read_d["intron1_IR"] and this_read_d["intron2_IR"]:
                this_read_d["both_IR"] = True
                this_read_d["intron1_IR"] = False
                this_read_d["intron2_IR"] = False

            
            # Let's extract the UPI
            upi_start = r2.sequence.find(upi_up)
            upi_end = r2.sequence.find(upi_down)
            upi = r2.sequence[upi_start + len(upi_up):upi_end]

            if len(upi) != 12:
                continue

            # And also the UMI
            umi = r1.sequence[16:21]

            if this_read_d["with_CE"]:

                ce_start2 = r2.sequence.find(start_exon_2_rc) + len(start_exon_2_rc)
                ce_end2 = r2.sequence.find(end_exon_1_rc) - 9
                ce_seq = r2.sequence[ce_start2:ce_end2]
                if len(ce_seq) == 154 and ''.join([rev_c(ce_seq)[a] for a in pos_f]) == ce_detect:
                    ce_seq = rev_c(ce_seq)
                else:
                    ce_seq = ""
            else:
                ce_seq = ""

            if upi in ce_seq_d.keys():
                ce_seq_d[upi].append(ce_seq)
            else:
                ce_seq_d[upi] = [ce_seq]

            splicing_csv.write(','.join([upi, umi, ','.join([str(int(this_read_d[a])) for a in categories]), ce_seq]) + "\n")


    #### Find consensus CE sequences ####
    with open(args.consensus_csv, "w") as file:
        file.write("upi,n,seq\n")

        for upi, seqs in ce_seq_d.items():
            consensus_seq = ""
            seqs_correct_length = [a for a in seqs if len(a) == 154]
            
            if len(seqs_correct_length) > 0:
                for i in range(154):
                    # Extract nucleotide i from all seqs
                    nts_this_pos = [a[i] for a in seqs_correct_length if a != ""]

                    # Find most common element
                    try:
                        max_nt = mode(nts_this_pos)
                        nt_frac = nts_this_pos.count(max_nt)/len(nts_this_pos)
                        if nt_frac >= args.consensus_frac:
                            consensus_seq += max_nt
                        else:
                            consensus_seq += "N"
                    except:
                        consensus_seq += "N"

            file.write(','.join([upi, str(len(seqs)), consensus_seq]) + "\n")



if __name__ == "__main__":
    main()

