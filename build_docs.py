import os
import webgestaltpy

"""
Auto-generation of the Reference section for the documentation website.

Needs to be run before website is generated. Currently, this script is run automatically in the github action.

Configure the order of the functions by changing the PRIORITY variable in the config section. 

To skip a function or variable, add the variable to SKIP.
"""


ADJUSTED_HEADERS = ["#" * i for i in range(1, 7)]

######################
####### CONFIG #######
######################

OUTPUT_DIR = "docs/reference"
PRIORITY = ["ora", "meta_ora", "gsea", "meta_gsea", "nta", "NTAMethod"]
SKIP = ["webgestaltpy"]


def sanitize_file_name(name) -> str:
    keepcharacters = (" ", ".", "_")
    return "".join(c for c in name if c.isalnum() or c in keepcharacters).rstrip()


def process_method(name):
    func = getattr(webgestaltpy, name)
    doc_string = func.__doc__
    if doc_string is None:
        return ""
    new_doc_string = ""
    for line in doc_string.split("\n"):
        splits = line.split(" ")
        for i in range(len(splits)):
            if splits[i].strip() in ADJUSTED_HEADERS:
                splits[i] += "#"
        new_doc_string += " ".join(splits) + "\n"
    partial_name = sanitize_file_name(name)
    filename = partial_name + ".md"
    with open(OUTPUT_DIR + "/" + filename, "w") as w:
        w.write("# `" + name + "`\n" + new_doc_string + "\n")
    print("Processed " + name + ".")
    return partial_name


def start():
    all_functions = dir(webgestaltpy)
    functions = [f for f in all_functions if "__" not in f and f not in SKIP]
    new_content = "# Reference\nDocumentation for each function in webgestaltpy.\n\n## Functions\n\n"
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    for p in PRIORITY:
        if p in functions:
            path = process_method(p)
            if path != "":
                new_content += f"- [{p}](./{path}.md)\n"
    for f in functions:
        if f not in PRIORITY:
            path = process_method(f)
            if path != "":
                new_content += f"- [{f}](./{path}.md)\n"
    with open(OUTPUT_DIR + "/index.md", "w") as w:
        w.write(new_content)


if __name__ == "__main__":
    start()
