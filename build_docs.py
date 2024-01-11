import os
import webgestaltpy

OUTPUT_DIR = "docs/reference"
ADJUSTED_HEADERS = ["#" * i for i in range(1, 7)]
print(ADJUSTED_HEADERS)
PRIORITY = ["webgestaltpy", "ora_from_files", "gsea_from_files"]


def sanitize_file_name(name) -> str:
    keepcharacters = (" ", ".", "_")
    return "".join(c for c in name if c.isalnum() or c in keepcharacters).rstrip()


def process_method(name):
    func = getattr(webgestaltpy, name)
    doc_string = func.__doc__
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
    functions = [f for f in all_functions if "__" not in f]
    new_content = "# Reference\nDocumentation for each function in webgestaltpy.\n\n## Functions\n\n"
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    for p in PRIORITY:
        if p in functions:
            path = process_method(p)
            new_content += f"- [{p}](./{path}.md)\n"
    for f in functions:
        if f not in PRIORITY:
            path = process_method(f)
            new_content += f"- [{f}](./{path}.md)\n"
    with open(OUTPUT_DIR + "/index.md", "w") as w:
        w.write(new_content)


if __name__ == "__main__":
    start()
