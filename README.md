# NOTE

**This is a development branch that breaks many features.** To install
the the previous release, use the [1.1 release], the [main branch], or
`pip install -U alhambra`.

# Introduction

Alhambra, formerly tilesetdesigner, is a package for designing DNA tile
systems.

This is the 2.0 prerelease, which is being built to support both SST and
DX tiles, and the ability to extensibly add other motifs.

For DX tiles, it uses stickydesign to create sticky end sequences, and
peppercompiler with spuriousSSM to create core sequences. For SSTs, it
will use DSD. It uses an extensible system for tileset design, and can
flexibly take inputs of YAML or similar formats.

# Installation

Alhambra is designed to be installed as a Python package.

To install the prerelease version, use

    pip install -U --pre alhambra

To install the latest git version (which may be broken), use

    pip install -U 'alhambra @ git+https://github.com/DNA-and-Natural-Algorithms-Group/alhambra@v2'

Alhambra 2 is designed to work with Python 3.9 or later.

To install development versions, you can check out this github
repository, and use `pip -e` or some other method for installing Python
packages. Multiple versions can be handled through `virtualenv`.

All Alhambra requirements should be handled through setuptools
dependencies, but note that xgrow and peppercompiler both rely on C code
being compiled, which may fail.

# Usage

[Documentation is available online on readthedocs.io]. In particular,
see [the tutorial]. It is also available in the docs/ folder.

Most user-facing functions are on the TileSet class.

# Questions

Please send any questions to Constantine Evans, at cevans@evanslabs.org
or cge@dna.caltech.edu.

  [1.1 release]: https://github.com/DNA-and-Natural-Algorithms-Group/alhambra/releases/tag/v1.1.0
  [main branch]: https://github.com/DNA-and-Natural-Algorithms-Group/alhambra/tree/main
  [Documentation is available online on readthedocs.io]: https://alhambra.readthedocs.io/en/latest/
  [the tutorial]: https://alhambra.readthedocs.io/en/latest/tutorial.html
