# Installing `opal` in a Conda R Environment (Windows)

This guide provides a robust installation pathway for users running `opal` within a Conda R environment (e.g., Miniforge or Anaconda) that uses optimised math libraries like OpenBLAS.

## 1. Why Use This Method?

- **Performance:** Accelerates dense matrix operations (likelihood components, simulations) up to 30x via OpenBLAS.
- **Isolation:** Keeps your fisheries modelling environment separate from other system R installations.
- **Stability:** Resolves common "DLL Hell" and "Procedure Not Found" errors caused by mixing Conda and CRAN binaries.

## 2. Prerequisites

1. **Conda Environment:** An active environment (e.g., `r-openblas`) with `r-base=4.3` or higher.
2. **Rtools 4.3:** Standalone Rtools 4.3 installed to the default `C:\rtools43`.
   - **Note:** Avoid installing Rtools via `conda install` as it frequently causes version conflicts with the R-base headers.

## 3. Advantages & Disadvantages

| Feature | Advantages | Disadvantages |
|---|---|---|
| **Speed** | Superior performance for dense linear algebra components. | Minimal impact on sparse Cholesky factorisation (TMB bottleneck). |
| **Management** | Clean environment isolation and easy environment cloning. | Requires manual "bridging" to see system compilers (Rtools). |
| **Compatibility** | Access to high-performance math libraries (AVX2/AVX-512). | Sensitive to "DLL locking" from background IDE services (VS Code). |

## 4. Installation Steps

Perform these steps in an Administrator terminal (Anaconda/Miniforge Prompt).

### Step 1: Align Conda Foundations

Ensure the asynchronous and foundational stack is managed by Conda to prevent session locks:

```bash
conda activate your_env_name
conda install -c conda-forge r-later r-promises r-processx r-lifecycle=1.0.5 --force-reinstall
```

### Step 2: Bridge R to Rtools

Launch R from the terminal (`R --vanilla`) and point it to your system compilers:

```r
# Add Rtools to the R session's PATH
Sys.setenv(PATH = paste("C:/rtools43/usr/bin", Sys.getenv("PATH"), sep = ";"))

# Verify Rtools is recognised (should return TRUE)
pkgbuild::has_build_tools(debug = TRUE)
```

### Step 3: Source Build & Package Installation

To ensure all DLLs are correctly mapped to your OpenBLAS pointers, build dependencies from source.

```r
# 1. Update core dependencies from source to fix DLL links
install.packages(c("rlang", "lifecycle", "plyr", "reshape2"), type = "source")

# 2. Install StanEstimators from the R-Universe
install.packages('StanEstimators', 
                 repos = c('https://andrjohns.r-universe.dev', 'https://cloud.r-project.org'),
                 type = "source")

# 3. Install the NOAA stack (use upgrade = "never" to protect your Conda links)
remotes::install_github("noaa-afsc/SparseNUTS", upgrade = "never")
remotes::install_github("N-DucharmeBarth-NOAA/opal", upgrade = "never")
```

## 5. Troubleshooting

### "The specified procedure could not be found" (DLL Error)

- **Cause:** A binary package (e.g., `plyr.dll`) was installed via CRAN, but your environment expects OpenBLAS-linked function signatures.
- **Solution:** Reinstall the failing package from source: `install.packages("package_name", type = "source")`.

### "Permission Denied" during installation

- **Cause:** A background service (like the VS Code R Language Service) is holding a package DLL in memory.
- **Solution:** Close VS Code and all R sessions. Open a fresh terminal as Administrator and run the install command before loading any libraries.

### R cannot find Rtools (`make` not found)

- **Cause:** Conda ignores system-wide compilers by default.
- **Solution:** Manually set the PATH inside R as shown in Step 2. You can make this permanent by adding the `Sys.setenv` line to your `.Rprofile`.

## 6. Verification

After installation, verify the build with:

```r
library(opal)
library(SparseNUTS)
```
