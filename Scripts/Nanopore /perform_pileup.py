import pysam
import argparse


def add_to_dict(d, key):
	if key in d.keys():
		d[key] += 1
	else:
		d[key] = 1
	return d


def make_key(rname, position, value):
	return ','.join([str(rname), str(position), str(value)])


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input")
	parser.add_argument("-o", "--output")
	parser.add_argument("-l", "--min_aligned_length", default=400, type=int)
	parser.add_argument("--primary", action="store_true", default=False)

	args = parser.parse_args()

	d = {}
	n_reads = 0
	with pysam.AlignmentFile(args.input, 'rb') as file:
		for record in file:
			if record.is_secondary and args.primary:
				continue

			if record.is_unmapped:
				continue

			positions = record.get_reference_positions(full_length=True)
			seq = record.query_sequence
			rname = record.reference_id

			first_pos = min([i for i in positions if i is not None])
			last_pos = max([i for i in positions if i is not None])

			if last_pos - first_pos < args.min_aligned_length:
				continue

			n_reads += 1

			started = False
			insertion_counter = 0

			for seq_pos, p in enumerate(positions):
				if p == first_pos:
					started = True
					prev_pos = first_pos

				if not started:
					continue

				if p is not None:  # therefore is digit
					insertion_counter = 0

					if p - prev_pos > 1:  # deletion in read
						for i in range(prev_pos+1, p):
							add_to_dict(d, make_key(rname, i, "del"))

					add_to_dict(d, make_key(rname, p, seq[seq_pos]))

					prev_pos = p

				if p is None:  # insertion in read
					insertion_counter += 1
					this_pos = prev_pos + 0.0001 * insertion_counter
					add_to_dict(d, make_key(rname, this_pos, seq[seq_pos]))

				if p == last_pos:
					break

	with open(args.output, 'w') as file:
		file.write("rname,position,nt,n,total_reads\n")
		for key, value in d.items():
			file.write(key + "," + str(value) + "," + str(n_reads) + "\n")






if __name__ == "__main__":
	main()
