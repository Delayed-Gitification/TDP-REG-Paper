import argparse
import pysam

def write_out(file, to_write):
	file.write(to_write)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-b", "--bam", required=True)
	parser.add_argument("--min_intron_length", default=50, type=int)
	parser.add_argument("--chunk_size", default=1_000_000, type=int)
	parser.add_argument("-o", "--output", required=True)
	parser.add_argument("--early_stop", default=-1, type=int)
	args = parser.parse_args()

	output_d = {}
	skipped = 0

	with pysam.AlignmentFile(args.bam) as bam:
		record_number = 0
		for record in bam:

			try:
				if record.is_unmapped:
					to_write = ','.join(["NA", "NA", "NA", "NA"])

					to_write += "," + ";".join([])

				else:

					if record.flag == 256:
						secondary = True
					else:
						secondary = False

					record_number += 1

					if record_number % 10_000 == 0:
						print(record_number)

					positions = record.get_reference_positions()

					first_pos = min(positions)
					last_pos = max(positions)

					mapping_quality = record.mapping_quality

					ref = record.reference_name
					junctions = []

					if record.is_reverse:
					 	strand = "-"
					else:
					 	strand = "+"

					for i in range(len(positions) - 1):
						distance = positions[i + 1] - positions[i]
						if distance >= args.min_intron_length:
							junctions.append(str(positions[i]) + "-" + str(positions[i + 1]))

					to_write = ','.join([ref, str(mapping_quality), str(secondary), strand, str(first_pos), str(last_pos)])

					to_write += "," + ";".join(junctions)

				if to_write in output_d.keys():
					output_d[to_write] += 1
				else:
					output_d[to_write] = 1


				if record_number > args.early_stop > 0:
					break
			except:
				skipped+=1

		with open(args.output, 'w') as out:
			out.write("reference,mapping_quality,secondary,strand,first_pos,last_pos,junctions,number_of_reads\n")
			for key, value in output_d.items():
				out.write(key + "," + str(value) + "\n")

	print(str(skipped) + " skipped")

if __name__ == '__main__':
	main()