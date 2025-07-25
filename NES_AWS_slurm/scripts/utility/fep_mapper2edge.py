import re

input_filename = 'core4_R1_cycle100.edge'
output_filename = 'edgelist.txt'
 
pattern = re.compile(r'# (\d+)-(\d+) -> (\d+)-(\d+)')

converted_lines = []

with open(input_filename, 'r') as file:
    for line in file:
        match = pattern.search(line)
        if match:
            int1, int2, int3, int4 = match.groups()
            converted_line = f"lig_{int1}-{int2},lig_{int3}-{int4}"
            converted_lines.append(converted_line)

with open(output_filename, 'w') as file:
    for converted_line in converted_lines:
        file.write(converted_line + '\n')

print(f"save '{output_filename}'")
