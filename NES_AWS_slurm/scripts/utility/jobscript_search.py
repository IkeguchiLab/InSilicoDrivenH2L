import os
import re
import argparse
import fnmatch

def search_files(directory, pattern, file_pattern="jobscript_*"):
    matching_files = []
    regex = re.compile(pattern)

    for root, dirs, files in os.walk(directory):
        for file in files:
            if fnmatch.fnmatch(file, file_pattern):
                file_path = os.path.join(root, file)
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                        if regex.search(content):
                            matching_files.append(file_path)
                except (UnicodeDecodeError, FileNotFoundError, PermissionError):
                    continue

    return matching_files


parser = argparse.ArgumentParser(description="seiki")
parser.add_argument("regex", type=str)
args = parser.parse_args()

directory_to_search = './'
search_pattern = args.regex

matching_files = sorted(search_files(directory_to_search, search_pattern))

if matching_files:
    print("contains the given pattern:")
    for file in matching_files:
        print(file)
else:
    print("cannot find any files containing the given pattern")
