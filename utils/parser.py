# import list of vcf files (paths)

import os


def get_files_names(path: str, ext: str) -> list:
    list_files = []
    for (dirname, _, files) in os.walk(path):
        for filename in files:
            if filename.endswith(ext):
                list_files.append(os.path.abspath(os.path.join(dirname, filename)))
    return list_files


