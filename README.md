# Big Bang for Space Settler #

[![Build Status](https://travis-ci.org/Razican/Space-Settler-BigBang.svg)](https://travis-ci.org/Razican/Space-Settler-BigBang)
[![Coverage Status](https://coveralls.io/repos/Razican/Space-Settler-BigBang/badge.svg)](https://coveralls.io/r/Razican/Space-Settler-BigBang)

This is the Rust implementation for generating a new galaxy in Space Settler. It will generate a
random galaxy taking into account real proportions of stars, planets and satellites to the best of
our knowledge. It will try to do as far as possible a true habitability check in rocky bodies, but
the main objective is to create as much habitable bodies as possible so that a future game could use
this database.

## License ##

This software is licensed under the 3-clause BSD license. You may copy, modify and redistribute it
provided that the conditions shown in the [LICENSE](LICENSE) file are met.

## Documentation ##

The documentation of this project is available in [GitHub pages](http://razican.github.io/Space-Settler-BigBang/). There you can find module by module documentation with some big explanations of how does the big bang work.

## References ##

**Stellar classification**:
 - https://en.wikipedia.org/wiki/Stellar_classification

**Chandrasekhar limit**:
 - https://en.wikipedia.org/wiki/Chandrasekhar_limit

**White Dwarf mass - radius relation**:
 - http://iopscience.iop.org/0004-637X/474/2/774/pdf/0004-637X_474_2_774.pdf

**O type star**:
 - https://en.wikipedia.org/wiki/O-type_main-sequence_star
