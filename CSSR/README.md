# CSSR #
CSSR (Coarse-grained Secondary Structure of RNAs) is an algorithm for
assignment of canonical secondary structure to both full-atomic and
coarse-grained 3D structures of RNAs.

## Installation ##
We provide pre-compiled binaries for 64bit Linux, 64bit MacOS and 64bit
Windows at ``exe/CSSR-Linux64``, ``exe/CSSR-Mac64`` and
``exe/CSSR-Win64.exe``, respectively. No compilation is needed.

For operating systems other than those mentioned above, you can compile
the program by:
```bash
make exe/CSSR
```

## Usage ##

```bash
exe/CSSR input.pdb output.dssr
```

Use the ``exe/CSSR -h`` option to get detail document.
Here, ``exe/CSSR`` should be changed to the binaries for your operating
systems if you use pre-compiled binaries.

## License ##

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

See ``src/README.md`` for dependencies on other libraries.

## Reference ##
Chengxin Zhang, Anna Marie Pyle (2022)
"CSSR: assignment of secondary structure to coarse-grained RNA tertiary structures."
Acta Crystallogr D. 78, 466-471.
