# Build Instructions

## Prerequisites

- GCC/G++ compiler (version 7 or later)
- Make
- Fortran compiler (gfortran) for OpenBLAS
- CBLAS and LAPACK development libraries

## Building on Linux

### Install dependencies

```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran libopenblas-dev liblapack-dev
```

### Build OpenBLAS

```bash
cd OpenBLAS
make -j4
make PREFIX=/usr/local install
cd ..
```

### Build and run tests

```bash
cd test_cblas_l3
make clean
make run_all
```

### Build benchmarks

```bash
cd custom_func
g++ -O3 -march=native -fopenmp -o benchmark benchmark.cpp -lstdc++ -lm -lopenblas
./benchmark
```

## Building on Windows (with MSYS2)

```bash
# Install MSYS2 packages
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-lapack make

# Build OpenBLAS
cd OpenBLAS
make -j4
cd ..

# Build tests
cd test_cblas_l3
make run_all
```

## CI/CD Integration

This project uses GitHub Actions for continuous integration. Every push to the main branch:

1. Installs all dependencies
2. Builds the project
3. Runs interface tests
4. Generates a build report

See `.github/workflows/ci.yml` for details.
