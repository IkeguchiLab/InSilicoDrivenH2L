import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
insert_file = sys.argv[3]

def remove_and_insert_pdb(input_file:str, output_file:str, insert_file:str) -> None:
    with open(input_file, 'r') as f:
        lines = f.readlines()

    index_list = []
    for index, line in enumerate(lines):
        if 'MOL' in line:
            index_list.append(index)

    #print(lines[:index_list[0]])
    #print(lines[index_list[-1]:])

    with open(insert_file, 'r') as f:
        insert_lines = f.readlines()

    new_lines = lines[:index_list[0]] + insert_lines + lines[index_list[-1] + 1:] 
    with open(output_file, 'w') as f:
        f.writelines(new_lines)

remove_and_insert_pdb(input_file, output_file, insert_file)
