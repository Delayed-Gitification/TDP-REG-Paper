import pysam
import sys

# This function checks that the risk SNP is present in the read, which is position 545 of the reference fasta

def main():

	snp_pos = 545

	bam_file = sys.argv[1]

	with pysam.AlignmentFile(bam_file, 'rb') as bam:
		
		out_bam = pysam.AlignmentFile(sys.argv[2], "wb", template=bam)
		
		for record in bam:
			positions = record.get_reference_positions(full_length=True)

			try:
				snp_pos_in_read = positions.index(snp_pos)
			except ValueError:
				continue

			try:

				nt = record.query_sequence[snp_pos_in_read]

				if nt == "C":
					out_bam.write(record)
			except:
				bad = 1


if __name__ == "__main__":
	main()


