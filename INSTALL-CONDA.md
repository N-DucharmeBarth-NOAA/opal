# Installing `opal` in a Conda R Environment (Windows)

This guide provides a robust installation pathway for users running `opal` within a Conda R environment (e.g., Miniforge or Anaconda) that uses optimised math libraries like OpenBLAS.

## 1. Why Use This Method?

- **Performance:** Accelerates dense matrix operations (likelihood components, simulations) up to 30x via OpenBLAS.
- **Isolation:** Keeps your fisheries modelling environment separate from other system R installations.
- **Stability:** Resolves common "DLL Hell" and "Procedure Not Found" errors caused by mixing Conda and CRAN binaries.

## 2. Prerequisites

1. **Conda Environment:** An active environment with `r-base=4.5` or higher (must match the `renv.lock` R version).
2. **Rtools 4.5:** Standalone Rtools 4.5 installed to the default `C:\rtools45`.
   - Download from: https://cran.r-project.org/bin/windows/Rtools/
   - **Note:** Avoid installing Rtools via `conda install` as it frequently causes version conflicts with the R-base headers.

## 3. Advantages & Disadvantages

| Feature | Advantages | Disadvantages |
|---|---|---|
| **Speed** | Superior performance for dense linear algebra components. | Minimal impact on sparse Cholesky factorisation (TMB bottleneck). |
| **Management** | Clean environment isolation and easy environment cloning. | Requires manual "bridging" to see system compilers (Rtools). |
| **Compatibility** | Access to high-performance math libraries (AVX2/AVX-512). | Sensitive to "DLL locking" from background IDE services (VS Code). |

## 4. Installation Steps

Perform these steps in an Administrator terminal (Anaconda/Miniforge Prompt).

### Step 1: Create the Conda Environment

Create a new environment with R 4.5, OpenBLAS, and required system libraries:

```bash
conda create -n r-openblas-45 -c conda-forge \
  r-base=4.5 \
  libblas=*=*openblas* \
  r-later r-promises r-processx \
  fftw
conda activate r-openblas-45
```

The `libblas=*=*openblas*` spec explicitly selects OpenBLAS-backed BLAS. `fftw` is required to compile RTMB.

> **PowerShell note:** `R` is an alias for `Invoke-History` in PowerShell. Always use `R.exe` to launch R sessions in PowerShell.

### Step 2: Bridge R to Rtools

Launch R from the terminal (`R.exe --vanilla`) and point it to your system compilers:

```r
# Add Rtools to the R session's PATH
Sys.setenv(PATH = paste("C:/rtools45/usr/bin", Sys.getenv("PATH"), sep = ";"))

# Verify Rtools is recognised (should return TRUE)
pkgbuild::has_build_tools(debug = TRUE)
```

### Step 3: Restore Packages via renv

This project uses `renv` for package management. The lockfile was generated with R 4.5.2.

Recommended packages (MASS, mgcv, etc.) ship with R/Conda and must be excluded from the renv restore to avoid version conflicts and DLL errors:

```r
# Exclude recommended packages (Conda manages these)
recommended_base <- c("MASS", "Matrix", "boot", "class", "cluster", "codetools",
                      "foreign", "KernSmooth", "lattice", "mgcv", "nlme",
                      "nnet", "rpart", "spatial", "survival")
rec_pkgs <- rownames(installed.packages(priority = "recommended"))
all_exclude <- union(rec_pkgs, recommended_base)

renv::restore(exclude = all_exclude)
```

After the restore, install mgcv as a binary (source compilation causes DLL relocation errors in Conda):

```r
install.packages("mgcv", type = "binary")
```

### Step 4: Install `opal` Locally

Build and install the local `opal` package:

```r
install.packages(".", repos = NULL, type = "source")
```

### Step 5: Sync the lockfile

Update the lockfile to record all installed packages:

```r
renv::snapshot()
```

## 5. Troubleshooting

### "The specified procedure could not be found" (DLL Error)

- **Cause:** A binary package (e.g., `plyr.dll`) was installed via CRAN, but your environment expects OpenBLAS-linked function signatures.
- **Solution:** Reinstall the failing package from source: `install.packages("package_name", type = "source")`.

### MASS version not found during `renv::restore()`

- **Cause:** The lockfile may record a Conda-patched MASS version (e.g., `7.3.60.0.1`) that does not exist on CRAN.
- **Solution:** Install MASS from CRAN and update the lockfile:
  ```r
  install.packages("MASS")
  renv::snapshot(packages = "MASS")
  ```

### RTMB fails to compile: `fftw3.h: No such file or directory`

- **Cause:** The FFTW library is not installed in the Conda environment.
- **Solution:** Install it from conda-forge:
  ```bash
  conda install -c conda-forge fftw
  ```

### mgcv fails with "32 bit pseudo relocation out of range"

- **Cause:** Source compilation of mgcv produces DLL address conflicts in the Conda environment.
- **Solution:** Install mgcv as a binary instead:
  ```r
  install.packages("mgcv", type = "binary")
  ```

### "Permission Denied" during installation

- **Cause:** A background service (like the VS Code R Language Service) is holding a package DLL in memory.
- **Solution:** Close VS Code and all R sessions. Open a fresh terminal as Administrator and run the install command before loading any libraries.

### R cannot find Rtools (`make` not found)

- **Cause:** Conda ignores system-wide compilers by default.
- **Solution:** Manually set the PATH inside R as shown in Step 2. You can make this permanent by adding the `Sys.setenv` line to your `.Rprofile`.

### `conda install r-base=4.5` fails with `libblas` conflict

- **Cause:** An existing environment has conflicting BLAS versions between the `defaults` and `conda-forge` channels.
- **Solution:** Create a fresh environment rather than upgrading an existing one (see Step 1), or use `--override-channels`:
  ```bash
  conda install -c conda-forge r-base=4.5 --override-channels
  ```

## 6. Verification

After installation, verify the build and OpenBLAS performance:

```r
library(opal)
library(SparseNUTS)

# Confirm OpenBLAS is active
extSoftVersion()[c("BLAS", "LAPACK")]

# Benchmark matrix operations (OpenBLAS: ~50ms, reference BLAS: ~30-60s)
n <- 2000
A <- matrix(rnorm(n * n), n, n)
B <- matrix(rnorm(n * n), n, n)
system.time(A %*% B)
```

Expected results with a correctly configured OpenBLAS environment (24-core machine):
- 2000×2000 matrix multiply: **< 100ms**
- 1000×1000 eigendecomposition: **< 0.5s**
