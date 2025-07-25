import re

def remove_section(content, section_name):
    pattern = fr'\[ {section_name} \][^\[]*(?=\[|$)'
    matches = re.finditer(pattern, content)
    for match in matches:
        start = match.start()
        end = match.end()
        content = content[:start] + content[end:]
    return content

def extract_and_remove_section(content, section_name):
    pattern = fr'\[ {section_name} \]([^\[]*)(?=\[|$)'
    match = re.search(pattern, content)
    if match:
        section_content = match.group(1)
        content = re.sub(pattern, '', content)
        return content, section_content
    else:
        return content, None

def process_ini_file(file_path, output_file):
    with open(file_path, 'r') as file:
        content = file.read()
    # Remove [defaults] section
    content = remove_section(content, 'defaults')

    # Extract and remove [atomtypes] section
    content, atomtypes_section = extract_and_remove_section(content, 'atomtypes')

    with open(output_file, 'w') as output:
        output.write(content)

    if atomtypes_section:
        with open('ffMOL.itp', 'w') as atomtypes_file:
            atomtypes_file.write("[ atomtypes ]\n")
            atomtypes_file.write(atomtypes_section)

 

if __name__ == "__main__":
    input_file = 'gromacs_sed.top'  
    output_file = 'MOL.itp'  
    try:
        process_ini_file(input_file, output_file)
    except Exception as e:
        print(f"Error processing the INI file: {e}")
