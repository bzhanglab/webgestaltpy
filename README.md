# WebGestaltPy

A python package containing bindings to the underlying Rust library [webgestalt_rust](https://www.github.com/bzhanglab/webgestalt_rust).

## Features

- Compute GSEA and ORA
- Run a meta-analysis combining multiple lists
- Combine multiple lists into a single analysis type

The output of the python package is the values. This does not generate any HTML reports. For reports, please use the [R package](https://github.com/bzhanglab/webgestaltr).

## Installation

```
pip install webgestaltpy
```

## Development

WebGestaltPy uses [maturin](https://www.maturin.rs) and [rye](https://rye-up.com) to build the full project. To build WebGestaltPy, run the following commands

```bash
git clone https://github.com/bzhanglab/webgestaltpy.git
cd webgestaltpy
rye sync # download correct python version and install dependencies
rye tools install maturin
rye run python -m ensurepip --default-pip # makes sure maturin can run correctly
rye shell # initialize shell so the python command executes the local version
maturin develop # add -r to build and install the full release version of the rust library.
python test.py # run the test script
```
