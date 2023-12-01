with open('nr_proteinID2annotation.txt', 'r') as input_file:
    with open('output.tsv', 'w') as output_file:
        for line in input_file:
            # Looking the first space appearing in each row
            first_space_index = line.find(' ')
            
            # Replace the space to tab
            if first_space_index != -1:
                modified_line = line[:first_space_index] + '\t' + line[first_space_index + 1:]
            else:
                modified_line = line

            output_file.write(modified_line)
