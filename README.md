# WebGestaltPy

 A python package containing bindings to the underlying Rust library [webgestalt_rust](https://www.github.com/bzhanglab/webgestalt_rust).

## Features

* Compute GSEA, ORA, and NTA
* Run a meta-analysis combining multiple lists
* Combine multiple lists into a single analysis type

The output of the python package is the values. This does not generate any HTML reports. For reports, please use the [R package](https://github.com/bzhanglab/webgestaltr).

## Installation

WebGestaltPy uses [maturin](https://www.maturin.rs) and [rye](https://rye-up.com) to build the full project. To build WebGestaltPy, run the following commands

```bash
rye sync
rye tools install maturin
rye run python -m ensurepip --default-pip # makes sure maturin can run correctly
rye shell # initialize shell so the python command executes the local version
maturin develop # -r # add -r to build the full release version of the rust library.
python test.py # run the test script
```


