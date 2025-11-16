# C source skeleton

## Prerequisites

To use all build targets.

- make version 4.0 or newer
- Some version of ctags
- clang-format for formatting

## Configuration

Configure the toolchain and additional options in `config.mk`:

```make
# Supported: GCC, CLANG, ICC
TOOLCHAIN ?= GCC
ENABLE_OPENMP ?= false

OPTIONS +=  -DARRAY_ALIGNMENT=64
#OPTIONS +=  -DVERBOSE_AFFINITY
#OPTIONS +=  -DVERBOSE_DATASIZE
#OPTIONS +=  -DVERBOSE_TIMER
```

The verbosity options enable detailed output about affinity settings, allocation
sizes and timer resolution.

## Build targets

- Build with:

```sh
make
```

You can build multiple toolchains in the same directory, but note that the
Makefile is only acting on the one currently set. Intermediate build results are
located in the `./build/<TOOLCHAIN>` directory.
This will also trigger the generation of a `.clangd` configuration file including
all includes and defines for correct functioning of clang language server.

To see the output of all executed commands use:

```sh
make Q=
```

- Clean up with:

```sh
make clean
```

to clean intermediate build results of the currently active toolchain. This will
not remove the binary as well as intermediate build results from other
toolchains. This is useful to force a rebuild for the currently active
toolchain.

```sh
make distclean
```

to clean **all** intermediate build results and the binary.

- Generate assembly output:

```sh
make asm
```

The assembler files will also be located in the `./build/<TOOLCHAIN>` directory.

- Output compiler version information:

```sh
make info
```

Useful in automated benchmarking scripts to output compiler version.

- Autoformat all source files according to configuration in `.clang-format`:

```sh
make format
```
