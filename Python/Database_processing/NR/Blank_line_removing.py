input_file_path = 'nr_proteinID2annotation.tsv'
output_file_path = 'output.txt'

with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
    lines = [line.strip() for line in input_file if line.strip()]
    output_file.write('\n'.join(lines))
