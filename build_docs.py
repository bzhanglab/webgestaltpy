import os
import glob

KEYWORDS = ["#[pyfunction]"]


class PythonFunction:
    name: str
    comments: str

    def __init__(self, name: str, comments: str) -> None:
        self.name = name
        self.comments = comments


def process_file(path):
    print(f"Processing {path}.")
    with open(path, "r") as r:
        lines = r.read().split("\n")
    lines = [l.strip() for l in lines if len(l.strip()) != 0]
    for i in range(len(lines)):
        line = lines[i]
        if line in KEYWORDS:
            print(f"Found {lines[i+1]}.")


def start():
    print("Building auto-generated docs")
    rust_files = glob.glob("src/*.rs")
    print(f"Found {len(rust_files)} Rust files.")
    for file in rust_files:
        process_file(file)


if __name__ == "__main__":
    start()
