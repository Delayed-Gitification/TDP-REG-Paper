from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import argparse
import numpy as np

paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(resource_filename('spliceai', x)) for x in paths]


def get_probs(input_sequence):
    context = 10000
    x = one_hot_encode('N' * (context // 2) + input_sequence + 'N' * (context // 2))[None, :]
    y = np.mean([models[m].predict(x) for m in range(5)], axis=0)

    acceptor_prob = y[0, :, 1]
    donor_prob = y[0, :, 2]

    return acceptor_prob, donor_prob



def make_dict(filename):
    d = {}
    with open(filename) as file:
        for i, line in enumerate(file):
            if i == 0:
                continue

            split = line.rstrip().split(",")

            upi = split[0]
            n = int(split[1])
            seq = split[2]

            d[upi] = [n, seq]

    return d


# random mutagenesis

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nt", type=str, required=True, help="consensus csv for nt cells")
    parser.add_argument("--dox", type=str, required=True, help="consensus csv for dox cells")
    parser.add_argument("--min_reads", type=int, default=4, help="Minimum number of reads for each cryptic")
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()

    nt_d = make_dict(args.nt)
    dox_d = make_dict(args.dox)

    # Filter the dictionary for UPIs which we obtained a CE sequence for (rather than just '""'')
    only_actual_seqs_nt_d = {a:b for a, b in nt_d.items() if len(b[1]) == 154}
    only_actual_seqs_dox_d = {a:b for a, b in dox_d.items() if len(b[1]) == 154}

    # Combine the two files and check they are consistent
    filtered_d = {}
    all_upis = list(set(only_actual_seqs_nt_d.keys()).union(set(only_actual_seqs_dox_d.keys())))

    for upi in all_upis:
        if upi in only_actual_seqs_nt_d.keys() and upi in only_actual_seqs_dox_d.keys():
            total_reads = nt_d[upi][0] + dox_d[upi][0]

            if total_reads >= args.min_reads:
                if only_actual_seqs_nt_d[upi][1] == only_actual_seqs_dox_d[upi][1]:
                    filtered_d[upi] = only_actual_seqs_nt_d[upi][1]
        
        elif upi in only_actual_seqs_nt_d.keys():
            if only_actual_seqs_nt_d[upi][0] >= args.min_reads:
                filtered_d[upi] = only_actual_seqs_nt_d[upi][1]

        elif upi in only_actual_seqs_dox_d.keys():
            if only_actual_seqs_dox_d[upi][0] >= args.min_reads:
                filtered_d[upi] = only_actual_seqs_dox_d[upi][1]

    print(len(filtered_d.keys()))

    # Make sure that each UPI is detected in both conditions and doesn't contain any Ns

    filtered_d2 = {a:b for a, b in filtered_d.items() if a in nt_d.keys() and a in dox_d.keys() and "N" not in b}
    print(len(filtered_d2))


    upstream = "GATCCTACCATCCACTCGACACACCCGCCAGCGGCCGCTTCTTGGTGCCAGCTTATCAtagcgctaccggtcgccaccatggCgagaACCATGGTAGCCATGGAGaccATGgggctcATGACAACAGATCTGGCAAAATTTGGGagatacaccggctggggcagGTAAGAATGCACATCACTTCTTGAGAGTATGGAGGAGTGAAATGACACTCAGTGCCAGAGTTACTGTATATCTACACTTTAAAAGTGTAGCTTTTAAAAGATAAGCAAGCACAATCTTTTGTGTGTGTGTGTGTGAATGTGTGTGTGTGTGTGTGTCACCCAG"
    downstream = "GTATGCATCACCCCCCCAGCTAATTTTTTTTTGTATTTTTTACCGAGTCGGGGTTTCGCAATGTTGCCCAGGCTGGTCTCAGAGTCTCGCTCTGTTGTCTACGCTGGAGTGCAGTAACATGAGCCACTGTGCCCGGCCAATCCTAAGAATTTCTTTTGCGGTGGTTGCAAGTCTGGGCAGAACTCTTGTCAGGGGCTGTAACTGGACTTATCTTTACTCCTTTGTCAGgtAtccggccagggcgatagcctgcAATCCTCCCAGGCAACATGAAGGATAACTTCTGNNNNNNNNNNNNGGATCCGCAGGCCTCTGCTAGCTTGACTGACTGAGATACAGCGTACCTTCAGCTCACAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGCTATTGCTTTATTTGTAACCATTATAAGCTGCAATAAACAAGTTAACAACAACAATTGCA"

    with open(args.output, 'w') as output:

        output.write("upi,acc,don\n")

        for upi, ce_seq in filtered_d2.items():
            seq = upstream.upper() + ce_seq + downstream.upper()
            
            acceptor_prob, donor_prob = get_probs(seq)

            ce_acc = acceptor_prob[327]
            ce_don = donor_prob[480]

            output.write(upi + "," + str(ce_acc) + "," + str(ce_don) + "\n")


if __name__ == "__main__":
    main()
