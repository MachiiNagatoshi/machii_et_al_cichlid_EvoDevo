# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 250526 Run PyWGCNA with time series transcriptome data
#
# reference
# - https://github.com/mortazavilab/PyWGCNA
# - https://mortazavilab.github.io/PyWGCNA/html/installation.html#install-from-pypi-recommended
# - https://anndata.readthedocs.io/en/latest/

# %% [markdown]
# ## install

# %% [markdown]
#
# ```bash
# (base) ~ ❯❯❯ pip install PyWGCNA
# Collecting PyWGCNA
#   Downloading pywgcna-2.2.1.tar.gz (54 kB)
#      ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 54.6/54.6 kB 4.7 MB/s eta 0:00:00
#   Preparing metadata (setup.py) ... done
# Requirement already satisfied: pandas>=2.1.0 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (2.2.3)
# Collecting numpy>=2.1.0 (from PyWGCNA)
#   Downloading numpy-2.2.6-cp310-cp310-macosx_10_9_x86_64.whl.metadata (62 kB)
#      ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 62.0/62.0 kB 9.0 MB/s eta 0:00:00
# Requirement already satisfied: scipy>=1.9.1 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (1.13.1)
# Requirement already satisfied: scikit-learn>=1.2.2 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (1.5.1)
# Requirement already satisfied: statsmodels>=0.14.0 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (0.14.2)
# Collecting matplotlib>=3.9.1 (from PyWGCNA)
#   Downloading matplotlib-3.10.3-cp310-cp310-macosx_10_12_x86_64.whl.metadata (11 kB)
# Requirement already satisfied: seaborn>=0.11.2 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (0.13.2)
# Collecting biomart>=0.9.2 (from PyWGCNA)
#   Downloading biomart-0.9.2-py3-none-any.whl.metadata (3.3 kB)
# Collecting gseapy>=1.1.3 (from PyWGCNA)
#   Downloading gseapy-1.1.8.tar.gz (112 kB)
#      ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 112.5/112.5 kB 15.9 MB/s eta 0:00:00
#   Installing build dependencies ... done
#   Getting requirements to build wheel ... done
#   Preparing metadata (pyproject.toml) ... done
# Collecting pyvis==0.3.1 (from PyWGCNA)
#   Downloading pyvis-0.3.1.tar.gz (748 kB)
#      ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 748.9/748.9 kB 32.8 MB/s eta 0:00:00
#   Preparing metadata (setup.py) ... done
# Requirement already satisfied: setuptools>=67.4.0 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (69.5.1)
# Collecting reactome2py>=3.0.0 (from PyWGCNA)
#   Downloading reactome2py-3.0.0-py3-none-any.whl.metadata (2.6 kB)
# Collecting anndata>=0.10.8 (from PyWGCNA)
#   Downloading anndata-0.11.4-py3-none-any.whl.metadata (9.3 kB)
# Requirement already satisfied: requests>=2.28.1 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (2.31.0)
# Requirement already satisfied: networkx>=2.8.4 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (3.3)
# Collecting rsrc>=0.1.3 (from PyWGCNA)
#   Downloading rsrc-0.1.3-py3-none-any.whl.metadata (5.6 kB)
# Requirement already satisfied: psutil>=5.9.0 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (5.9.8)
# Requirement already satisfied: jinja2>=2.9.6 in ./mambaforge/lib/python3.10/site-packages (from pyvis==0.3.1->PyWGCNA) (3.1.4)
# Requirement already satisfied: ipython>=5.3.0 in ./mambaforge/lib/python3.10/site-packages (from pyvis==0.3.1->PyWGCNA) (8.25.0)
# Collecting jsonpickle>=1.4.1 (from pyvis==0.3.1->PyWGCNA)
#   Downloading jsonpickle-4.0.5-py3-none-any.whl.metadata (8.2 kB)
# Collecting array-api-compat!=1.5,>1.4 (from anndata>=0.10.8->PyWGCNA)
#   Downloading array_api_compat-1.12.0-py3-none-any.whl.metadata (2.5 kB)
# Requirement already satisfied: exceptiongroup in ./mambaforge/lib/python3.10/site-packages (from anndata>=0.10.8->PyWGCNA) (1.2.0)
# Collecting h5py>=3.7 (from anndata>=0.10.8->PyWGCNA)
#   Downloading h5py-3.13.0-cp310-cp310-macosx_10_9_x86_64.whl.metadata (2.5 kB)
# Collecting natsort (from anndata>=0.10.8->PyWGCNA)
#   Downloading natsort-8.4.0-py3-none-any.whl.metadata (21 kB)
# Collecting packaging>=24.2 (from anndata>=0.10.8->PyWGCNA)
#   Using cached packaging-25.0-py3-none-any.whl.metadata (3.3 kB)
# Requirement already satisfied: contourpy>=1.0.1 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (1.2.1)
# Requirement already satisfied: cycler>=0.10 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (0.12.1)
# Requirement already satisfied: fonttools>=4.22.0 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (4.53.0)
# Requirement already satisfied: kiwisolver>=1.3.1 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (1.4.5)
# Requirement already satisfied: pillow>=8 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (10.3.0)
# Requirement already satisfied: pyparsing>=2.3.1 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (3.1.2)
# Requirement already satisfied: python-dateutil>=2.7 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (2.9.0)
# Requirement already satisfied: pytz>=2020.1 in ./mambaforge/lib/python3.10/site-packages (from pandas>=2.1.0->PyWGCNA) (2024.1)
# Requirement already satisfied: tzdata>=2022.7 in ./mambaforge/lib/python3.10/site-packages (from pandas>=2.1.0->PyWGCNA) (2024.1)
# Requirement already satisfied: json5>=0.8.4 in ./mambaforge/lib/python3.10/site-packages (from reactome2py>=3.0.0->PyWGCNA) (0.9.25)
# Requirement already satisfied: charset-normalizer<4,>=2 in ./mambaforge/lib/python3.10/site-packages (from requests>=2.28.1->PyWGCNA) (3.3.2)
# Requirement already satisfied: idna<4,>=2.5 in ./mambaforge/lib/python3.10/site-packages (from requests>=2.28.1->PyWGCNA) (3.6)
# Requirement already satisfied: urllib3<3,>=1.21.1 in ./mambaforge/lib/python3.10/site-packages (from requests>=2.28.1->PyWGCNA) (2.2.1)
# Requirement already satisfied: certifi>=2017.4.17 in ./mambaforge/lib/python3.10/site-packages (from requests>=2.28.1->PyWGCNA) (2025.4.26)
# Collecting memoir>=0.0.3 (from rsrc>=0.1.3->PyWGCNA)
#   Downloading memoir-0.0.3-py3-none-any.whl.metadata (5.5 kB)
# Collecting reprit>=0.3.0 (from rsrc>=0.1.3->PyWGCNA)
#   Downloading reprit-0.9.0-py3-none-any.whl.metadata (7.6 kB)
# Requirement already satisfied: joblib>=1.2.0 in ./mambaforge/lib/python3.10/site-packages (from scikit-learn>=1.2.2->PyWGCNA) (1.4.2)
# Requirement already satisfied: threadpoolctl>=3.1.0 in ./mambaforge/lib/python3.10/site-packages (from scikit-learn>=1.2.2->PyWGCNA) (3.5.0)
# Requirement already satisfied: patsy>=0.5.6 in ./mambaforge/lib/python3.10/site-packages (from statsmodels>=0.14.0->PyWGCNA) (0.5.6)
# Requirement already satisfied: decorator in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (5.1.1)
# Requirement already satisfied: jedi>=0.16 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.19.1)
# Requirement already satisfied: matplotlib-inline in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.1.7)
# Requirement already satisfied: prompt-toolkit<3.1.0,>=3.0.41 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (3.0.42)
# Requirement already satisfied: pygments>=2.4.0 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (2.18.0)
# Requirement already satisfied: stack-data in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.6.2)
# Requirement already satisfied: traitlets>=5.13.0 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (5.14.3)
# Requirement already satisfied: typing-extensions>=4.6 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (4.12.1)
# Requirement already satisfied: pexpect>4.3 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (4.9.0)
# Requirement already satisfied: MarkupSafe>=2.0 in ./mambaforge/lib/python3.10/site-packages (from jinja2>=2.9.6->pyvis==0.3.1->PyWGCNA) (2.1.5)
# Requirement already satisfied: six in ./mambaforge/lib/python3.10/site-packages (from patsy>=0.5.6->statsmodels>=0.14.0->PyWGCNA) (1.16.0)
# Requirement already satisfied: parso<0.9.0,>=0.8.3 in ./mambaforge/lib/python3.10/site-packages (from jedi>=0.16->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.8.4)
# Requirement already satisfied: ptyprocess>=0.5 in ./mambaforge/lib/python3.10/site-packages (from pexpect>4.3->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.7.0)
# Requirement already satisfied: wcwidth in ./mambaforge/lib/python3.10/site-packages (from prompt-toolkit<3.1.0,>=3.0.41->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.2.13)
# Requirement already satisfied: executing>=1.2.0 in ./mambaforge/lib/python3.10/site-packages (from stack-data->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (2.0.1)
# Requirement already satisfied: asttokens>=2.1.0 in ./mambaforge/lib/python3.10/site-packages (from stack-data->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (2.4.1)
# Requirement already satisfied: pure-eval in ./mambaforge/lib/python3.10/site-packages (from stack-data->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.2.2)
# Downloading anndata-0.11.4-py3-none-any.whl (144 kB)
#    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 144.5/144.5 kB 20.0 MB/s eta 0:00:00
# Downloading biomart-0.9.2-py3-none-any.whl (12 kB)
# Downloading matplotlib-3.10.3-cp310-cp310-macosx_10_12_x86_64.whl (8.2 MB)
#    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 8.2/8.2 MB 60.2 MB/s eta 0:00:00
# Downloading numpy-2.2.6-cp310-cp310-macosx_10_9_x86_64.whl (21.2 MB)
#    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 21.2/21.2 MB 62.4 MB/s eta 0:00:00
# Downloading reactome2py-3.0.0-py3-none-any.whl (26 kB)
# Downloading rsrc-0.1.3-py3-none-any.whl (7.0 kB)
# Downloading array_api_compat-1.12.0-py3-none-any.whl (58 kB)
#    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 58.2/58.2 kB 6.0 MB/s eta 0:00:00
# Downloading h5py-3.13.0-cp310-cp310-macosx_10_9_x86_64.whl (3.4 MB)
#    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 3.4/3.4 MB 54.0 MB/s eta 0:00:00
# Downloading jsonpickle-4.0.5-py3-none-any.whl (46 kB)
#    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 46.4/46.4 kB 17.0 MB/s eta 0:00:00
# Downloading memoir-0.0.3-py3-none-any.whl (4.9 kB)
# Using cached packaging-25.0-py3-none-any.whl (66 kB)
# Downloading reprit-0.9.0-py3-none-any.whl (10 kB)
# Downloading natsort-8.4.0-py3-none-any.whl (38 kB)
# Building wheels for collected packages: PyWGCNA, pyvis, gseapy
#   Building wheel for PyWGCNA (setup.py) ... done
#   Created wheel for PyWGCNA: filename=PyWGCNA-2.2.1-py3-none-any.whl size=54640 sha256=c76bb03750f445b6dc10136b958c32838305291614e1d7dc0e644d94ac7d9abb
#   Stored in directory: /Users/machiinagatoshi/Library/Caches/pip/wheels/ed/3a/5c/1bca9dc6d064986abade541de2cb640d8d1bfa8a6dfd66fea8
#   Building wheel for pyvis (setup.py) ... done
#   Created wheel for pyvis: filename=pyvis-0.3.1-py3-none-any.whl size=755829 sha256=250c8cb8b96855572bb66b20f50b0573efabd56e9eb8f90c55b523f5ff93f2e0
#   Stored in directory: /Users/machiinagatoshi/Library/Caches/pip/wheels/37/f9/93/44dd6cbfb2ead35307b114d27af8c6a14d5762a462af1e04f5
#   Building wheel for gseapy (pyproject.toml) ... error
#   error: subprocess-exited-with-error
#   
#   × Building wheel for gseapy (pyproject.toml) did not run successfully.
#   │ exit code: 1
#   ╰─> [60 lines of output]
#       <string>:10: SetuptoolsDeprecationWarning: The test command is disabled and references to it are deprecated.
#       !!
#       
#               ********************************************************************************
#               Please remove any references to `setuptools.command.test` in all supported versions of the affected package.
#       
#               This deprecation is overdue, please update your project and remove deprecated
#               calls to avoid build errors in the future.
#               ********************************************************************************
#       
#       !!
#       /private/var/folders/cy/fg1wkjzd1fq6mcb6y258n28m0000gn/T/pip-build-env-1lj9d53_/overlay/lib/python3.10/site-packages/setuptools/_distutils/dist.py:289: UserWarning: Unknown distribution option: 'tests_require'
#         warnings.warn(msg)
#       /private/var/folders/cy/fg1wkjzd1fq6mcb6y258n28m0000gn/T/pip-build-env-1lj9d53_/overlay/lib/python3.10/site-packages/setuptools/dist.py:759: SetuptoolsDeprecationWarning: License classifiers are deprecated.
#       !!
#       
#               ********************************************************************************
#               Please consider removing the following classifiers in favor of a SPDX license expression:
#       
#               License :: OSI Approved :: MIT License
#       
#               See https://packaging.python.org/en/latest/guides/writing-pyproject-toml/#license for details.
#               ********************************************************************************
#       
#       !!
#         self._finalize_license_expression()
#       running bdist_wheel
#       running build
#       running build_py
#       creating build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/algorithm.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/ssgsea.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/plot.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/enrichr.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/msigdb.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/__init__.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/parser.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/scipalette.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/utils.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/stats.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/gsva.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/__main__.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/base.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/gsea.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       copying gseapy/biomart.py -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy
#       creating build/lib.macosx-10.9-x86_64-cpython-310/gseapy/data
#       copying gseapy/data/palette.json -> build/lib.macosx-10.9-x86_64-cpython-310/gseapy/data
#       running build_ext
#       running build_rust
#       error: can't find Rust compiler
#       
#       If you are using an outdated pip version, it is possible a prebuilt wheel is available for this package but pip is not able to install from it. Installing from the wheel would avoid the need for a Rust compiler.
#       
#       To update pip, run:
#       
#           pip install --upgrade pip
#       
#       and then retry package installation.
#       
#       If you did intend to build this package from source, try installing a Rust compiler from your system package manager and ensure it is on the PATH during installation. Alternatively, rustup (available at https://rustup.rs) is the recommended way to download and update the Rust compiler toolchain.
#       [end of output]
#   
#   note: This error originates from a subprocess, and is likely not a problem with pip.
#   ERROR: Failed building wheel for gseapy
# Successfully built PyWGCNA pyvis
# Failed to build gseapy
# ERROR: Could not build wheels for gseapy, which is required to install pyproject.toml-based projects
# (base) ~ ❯❯❯ conda install rust                                                                                                                                             ✘ 1 
# Retrieving notices: done
# Channels:
#  - conda-forge
#  - bioconda
# Platform: osx-64
# Collecting package metadata (repodata.json): done
# Solving environment: - failed
#
# CondaError: KeyboardInterrupt
#
# ^C%                                                                                                                                                                              (base) ~ ❯❯❯ conda install rust                                                                                                                                           ✘ 130 
# Channels:
#  - conda-forge
#  - bioconda
# Platform: osx-64
# Collecting package metadata (repodata.json): done
# Solving environment: done
#
#
# ==> WARNING: A newer version of conda exists. <==
#     current version: 24.11.3
#     latest version: 25.3.1
#
# Please update conda by running
#
#     $ conda update -n base -c conda-forge conda
#
#
#
# ## Package Plan ##
#
#   environment location: /Users/machiinagatoshi/mambaforge
#
#   added / updated specs:
#     - rust
#
#
# The following packages will be downloaded:
#
#     package                    |            build
#     ---------------------------|-----------------
#     dask-2025.5.1              |     pyhd8ed1ab_0           8 KB  conda-forge
#     dask-core-2025.5.1         |     pyhd8ed1ab_0         971 KB  conda-forge
#     distributed-2025.5.1       |     pyhd8ed1ab_0         781 KB  conda-forge
#     rust-1.87.0                |       h34a2095_0       230.1 MB  conda-forge
#     rust-std-x86_64-apple-darwin-1.87.0|       h38e4360_0        33.4 MB  conda-forge
#     ------------------------------------------------------------
#                                            Total:       265.3 MB
#
# The following NEW packages will be INSTALLED:
#
#   rust               conda-forge/osx-64::rust-1.87.0-h34a2095_0 
#   rust-std-x86_64-a~ conda-forge/noarch::rust-std-x86_64-apple-darwin-1.87.0-h38e4360_0 
#
# The following packages will be UPDATED:
#
#   dask                                2025.5.0-pyhd8ed1ab_0 --> 2025.5.1-pyhd8ed1ab_0 
#   dask-core                           2025.5.0-pyhd8ed1ab_0 --> 2025.5.1-pyhd8ed1ab_0 
#   distributed                         2025.5.0-pyhd8ed1ab_0 --> 2025.5.1-pyhd8ed1ab_0 
#
#
# Proceed ([y]/n)? y
#
#
# Downloading and Extracting Packages:
#                                                                                                                                                                                  
# Preparing transaction: done                                                                                                                                                      
# Verifying transaction: done                                                                                                                                                      
# Executing transaction: done                                                                                                                                                      
# (base) ~ ❯❯❯ pip install PyWGCNA
# Collecting PyWGCNA
#   Using cached PyWGCNA-2.2.1-py3-none-any.whl
# Requirement already satisfied: pandas>=2.1.0 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (2.2.3)
# Collecting numpy>=2.1.0 (from PyWGCNA)
#   Using cached numpy-2.2.6-cp310-cp310-macosx_10_9_x86_64.whl.metadata (62 kB)
# Requirement already satisfied: scipy>=1.9.1 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (1.13.1)
# Requirement already satisfied: scikit-learn>=1.2.2 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (1.5.1)
# Requirement already satisfied: statsmodels>=0.14.0 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (0.14.2)
# Collecting matplotlib>=3.9.1 (from PyWGCNA)
#   Using cached matplotlib-3.10.3-cp310-cp310-macosx_10_12_x86_64.whl.metadata (11 kB)
# Requirement already satisfied: seaborn>=0.11.2 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (0.13.2)
# Collecting biomart>=0.9.2 (from PyWGCNA)
#   Using cached biomart-0.9.2-py3-none-any.whl.metadata (3.3 kB)
# Collecting gseapy>=1.1.3 (from PyWGCNA)
#   Using cached gseapy-1.1.8.tar.gz (112 kB)
#   Installing build dependencies ... done
#   Getting requirements to build wheel ... done
#   Preparing metadata (pyproject.toml) ... done
# Collecting pyvis==0.3.1 (from PyWGCNA)
#   Using cached pyvis-0.3.1-py3-none-any.whl
# Requirement already satisfied: setuptools>=67.4.0 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (69.5.1)
# Collecting reactome2py>=3.0.0 (from PyWGCNA)
#   Using cached reactome2py-3.0.0-py3-none-any.whl.metadata (2.6 kB)
# Collecting anndata>=0.10.8 (from PyWGCNA)
#   Using cached anndata-0.11.4-py3-none-any.whl.metadata (9.3 kB)
# Requirement already satisfied: requests>=2.28.1 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (2.31.0)
# Requirement already satisfied: networkx>=2.8.4 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (3.3)
# Collecting rsrc>=0.1.3 (from PyWGCNA)
#   Using cached rsrc-0.1.3-py3-none-any.whl.metadata (5.6 kB)
# Requirement already satisfied: psutil>=5.9.0 in ./mambaforge/lib/python3.10/site-packages (from PyWGCNA) (5.9.8)
# Requirement already satisfied: jinja2>=2.9.6 in ./mambaforge/lib/python3.10/site-packages (from pyvis==0.3.1->PyWGCNA) (3.1.4)
# Requirement already satisfied: ipython>=5.3.0 in ./mambaforge/lib/python3.10/site-packages (from pyvis==0.3.1->PyWGCNA) (8.25.0)
# Collecting jsonpickle>=1.4.1 (from pyvis==0.3.1->PyWGCNA)
#   Using cached jsonpickle-4.0.5-py3-none-any.whl.metadata (8.2 kB)
# Collecting array-api-compat!=1.5,>1.4 (from anndata>=0.10.8->PyWGCNA)
#   Using cached array_api_compat-1.12.0-py3-none-any.whl.metadata (2.5 kB)
# Requirement already satisfied: exceptiongroup in ./mambaforge/lib/python3.10/site-packages (from anndata>=0.10.8->PyWGCNA) (1.2.0)
# Collecting h5py>=3.7 (from anndata>=0.10.8->PyWGCNA)
#   Using cached h5py-3.13.0-cp310-cp310-macosx_10_9_x86_64.whl.metadata (2.5 kB)
# Collecting natsort (from anndata>=0.10.8->PyWGCNA)
#   Using cached natsort-8.4.0-py3-none-any.whl.metadata (21 kB)
# Collecting packaging>=24.2 (from anndata>=0.10.8->PyWGCNA)
#   Using cached packaging-25.0-py3-none-any.whl.metadata (3.3 kB)
# Requirement already satisfied: contourpy>=1.0.1 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (1.2.1)
# Requirement already satisfied: cycler>=0.10 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (0.12.1)
# Requirement already satisfied: fonttools>=4.22.0 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (4.53.0)
# Requirement already satisfied: kiwisolver>=1.3.1 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (1.4.5)
# Requirement already satisfied: pillow>=8 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (10.3.0)
# Requirement already satisfied: pyparsing>=2.3.1 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (3.1.2)
# Requirement already satisfied: python-dateutil>=2.7 in ./mambaforge/lib/python3.10/site-packages (from matplotlib>=3.9.1->PyWGCNA) (2.9.0)
# Requirement already satisfied: pytz>=2020.1 in ./mambaforge/lib/python3.10/site-packages (from pandas>=2.1.0->PyWGCNA) (2024.1)
# Requirement already satisfied: tzdata>=2022.7 in ./mambaforge/lib/python3.10/site-packages (from pandas>=2.1.0->PyWGCNA) (2024.1)
# Requirement already satisfied: json5>=0.8.4 in ./mambaforge/lib/python3.10/site-packages (from reactome2py>=3.0.0->PyWGCNA) (0.9.25)
# Requirement already satisfied: charset-normalizer<4,>=2 in ./mambaforge/lib/python3.10/site-packages (from requests>=2.28.1->PyWGCNA) (3.3.2)
# Requirement already satisfied: idna<4,>=2.5 in ./mambaforge/lib/python3.10/site-packages (from requests>=2.28.1->PyWGCNA) (3.6)
# Requirement already satisfied: urllib3<3,>=1.21.1 in ./mambaforge/lib/python3.10/site-packages (from requests>=2.28.1->PyWGCNA) (2.2.1)
# Requirement already satisfied: certifi>=2017.4.17 in ./mambaforge/lib/python3.10/site-packages (from requests>=2.28.1->PyWGCNA) (2025.4.26)
# Collecting memoir>=0.0.3 (from rsrc>=0.1.3->PyWGCNA)
#   Using cached memoir-0.0.3-py3-none-any.whl.metadata (5.5 kB)
# Collecting reprit>=0.3.0 (from rsrc>=0.1.3->PyWGCNA)
#   Using cached reprit-0.9.0-py3-none-any.whl.metadata (7.6 kB)
# Requirement already satisfied: joblib>=1.2.0 in ./mambaforge/lib/python3.10/site-packages (from scikit-learn>=1.2.2->PyWGCNA) (1.4.2)
# Requirement already satisfied: threadpoolctl>=3.1.0 in ./mambaforge/lib/python3.10/site-packages (from scikit-learn>=1.2.2->PyWGCNA) (3.5.0)
# Requirement already satisfied: patsy>=0.5.6 in ./mambaforge/lib/python3.10/site-packages (from statsmodels>=0.14.0->PyWGCNA) (0.5.6)
# Requirement already satisfied: decorator in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (5.1.1)
# Requirement already satisfied: jedi>=0.16 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.19.1)
# Requirement already satisfied: matplotlib-inline in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.1.7)
# Requirement already satisfied: prompt-toolkit<3.1.0,>=3.0.41 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (3.0.42)
# Requirement already satisfied: pygments>=2.4.0 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (2.18.0)
# Requirement already satisfied: stack-data in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.6.2)
# Requirement already satisfied: traitlets>=5.13.0 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (5.14.3)
# Requirement already satisfied: typing-extensions>=4.6 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (4.12.1)
# Requirement already satisfied: pexpect>4.3 in ./mambaforge/lib/python3.10/site-packages (from ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (4.9.0)
# Requirement already satisfied: MarkupSafe>=2.0 in ./mambaforge/lib/python3.10/site-packages (from jinja2>=2.9.6->pyvis==0.3.1->PyWGCNA) (2.1.5)
# Requirement already satisfied: six in ./mambaforge/lib/python3.10/site-packages (from patsy>=0.5.6->statsmodels>=0.14.0->PyWGCNA) (1.16.0)
# Requirement already satisfied: parso<0.9.0,>=0.8.3 in ./mambaforge/lib/python3.10/site-packages (from jedi>=0.16->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.8.4)
# Requirement already satisfied: ptyprocess>=0.5 in ./mambaforge/lib/python3.10/site-packages (from pexpect>4.3->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.7.0)
# Requirement already satisfied: wcwidth in ./mambaforge/lib/python3.10/site-packages (from prompt-toolkit<3.1.0,>=3.0.41->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.2.13)
# Requirement already satisfied: executing>=1.2.0 in ./mambaforge/lib/python3.10/site-packages (from stack-data->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (2.0.1)
# Requirement already satisfied: asttokens>=2.1.0 in ./mambaforge/lib/python3.10/site-packages (from stack-data->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (2.4.1)
# Requirement already satisfied: pure-eval in ./mambaforge/lib/python3.10/site-packages (from stack-data->ipython>=5.3.0->pyvis==0.3.1->PyWGCNA) (0.2.2)
# Using cached anndata-0.11.4-py3-none-any.whl (144 kB)
# Using cached biomart-0.9.2-py3-none-any.whl (12 kB)
# Using cached matplotlib-3.10.3-cp310-cp310-macosx_10_12_x86_64.whl (8.2 MB)
# Using cached numpy-2.2.6-cp310-cp310-macosx_10_9_x86_64.whl (21.2 MB)
# Using cached reactome2py-3.0.0-py3-none-any.whl (26 kB)
# Using cached rsrc-0.1.3-py3-none-any.whl (7.0 kB)
# Using cached array_api_compat-1.12.0-py3-none-any.whl (58 kB)
# Using cached h5py-3.13.0-cp310-cp310-macosx_10_9_x86_64.whl (3.4 MB)
# Using cached jsonpickle-4.0.5-py3-none-any.whl (46 kB)
# Using cached memoir-0.0.3-py3-none-any.whl (4.9 kB)
# Using cached packaging-25.0-py3-none-any.whl (66 kB)
# Using cached reprit-0.9.0-py3-none-any.whl (10 kB)
# Using cached natsort-8.4.0-py3-none-any.whl (38 kB)
# Building wheels for collected packages: gseapy
#   Building wheel for gseapy (pyproject.toml) ... done
#   Created wheel for gseapy: filename=gseapy-1.1.8-cp310-cp310-macosx_10_12_x86_64.whl size=552372 sha256=842471822956645cbc2c078e9c755cddaa18976ae8dcaec8705fd6230454ae97
#   Stored in directory: /Users/machiinagatoshi/Library/Caches/pip/wheels/0d/7c/45/8183c9fc419833a60eed2a863d7a3eba93174dd880c5f6f90c
# Successfully built gseapy
# Installing collected packages: reprit, packaging, numpy, natsort, memoir, jsonpickle, array-api-compat, rsrc, h5py, biomart, reactome2py, matplotlib, anndata, pyvis, gseapy, PyWGCNA
#   Attempting uninstall: packaging
#     Found existing installation: packaging 24.0
#     Uninstalling packaging-24.0:
#       Successfully uninstalled packaging-24.0
#   Attempting uninstall: numpy
#     Found existing installation: numpy 1.26.4
#     Uninstalling numpy-1.26.4:
#       Successfully uninstalled numpy-1.26.4
#   Attempting uninstall: matplotlib
#     Found existing installation: matplotlib 3.8.4
#     Uninstalling matplotlib-3.8.4:
#       Successfully uninstalled matplotlib-3.8.4
# ERROR: pip's dependency resolver does not currently take into account all the packages that are installed. This behaviour is the source of the following dependency conflicts.
# liftoff 1.6.3 requires biopython==1.76, but you have biopython 1.83 which is incompatible.
# liftoff 1.6.3 requires gffutils==0.10.1, but you have gffutils 0.13 which is incompatible.
# liftoff 1.6.3 requires interlap==0.2.6, but you have interlap 0.2.7 which is incompatible.
# liftoff 1.6.3 requires networkx==2.4, but you have networkx 3.3 which is incompatible.
# liftoff 1.6.3 requires numpy==1.21.0, but you have numpy 2.2.6 which is incompatible.
# liftoff 1.6.3 requires parasail==1.2.1, but you have parasail 1.3.4 which is incompatible.
# liftoff 1.6.3 requires pyfaidx==0.5.8, but you have pyfaidx 0.8.1.1 which is incompatible.
# liftoff 1.6.3 requires pysam==0.16.0.1, but you have pysam 0.22.1 which is incompatible.
# liftoff 1.6.3 requires ujson==3.2.0, but you have ujson 5.10.0 which is incompatible.
# Successfully installed PyWGCNA-2.2.1 anndata-0.11.4 array-api-compat-1.12.0 biomart-0.9.2 gseapy-1.1.8 h5py-3.13.0 jsonpickle-4.0.5 matplotlib-3.10.3 memoir-0.0.3 natsort-8.4.0 numpy-2.2.6 packaging-25.0 pyvis-0.3.1 reactome2py-3.0.0 reprit-0.9.0 rsrc-0.1.3
# ```

# %% [markdown]
# ## Execution

# %%
import PyWGCNA
import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import numpy as np

# %%
## read data
path_expression = '/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/LASY-seq/seurat_result/250520/raw/seurat_normalized_data.tsv'
path_sample_metadata = '/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/LASY-seq/seurat_result/250520/raw/seurat_metadata.tsv'
path_gene_metadata = '/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/LASY-seq/seurat_result/250520/raw/Maylandia_zebra.M_zebra_UMD2a.114.gene2transcript_handedit.csv'

df_expression = pd.read_csv(path_expression, sep='\t')
df_sample_metadata = pd.read_csv(path_sample_metadata, sep='\t')
df_gene_metadata = pd.read_csv(path_gene_metadata)

df_expression.iloc[:,1:].applymap(lambda x: math.log10(x+1)).hist(figsize=(20,20))

df_expression = df_expression.set_index('Unnamed: 0')
df_expression = df_expression.T

df_sample_metadata = df_sample_metadata.set_index('column_name')

df_gene_metadata = df_gene_metadata.drop_duplicates(subset='GeneID')
df_gene_metadata = df_gene_metadata.set_index('GeneID')

pyWGCNA_LASYseq = PyWGCNA.WGCNA(name = 'LASY-seq', geneExp = df_expression, sampleInfo = df_sample_metadata, geneInfo = df_gene_metadata,
                                RsquaredCut=0.90, TOMType='signed', networkType='signed', MEDissThres=0.1, TPMcutoff=0.5, minModuleSize=20, species='cichlids', powers=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])

pyWGCNA_LASYseq.geneExpr.to_df().head(5)

# %%
# Preprocessing
print("前処理前:")
print(f"サンプル数: {pyWGCNA_LASYseq.datExpr.n_obs}")
print(f"遺伝子数: {pyWGCNA_LASYseq.datExpr.n_vars}")
print(f"データ形状: {pyWGCNA_LASYseq.datExpr.shape}")

pyWGCNA_LASYseq.preprocess()

print("\n前処理後:")
print(f"サンプル数: {pyWGCNA_LASYseq.datExpr.n_obs}")
print(f"遺伝子数: {pyWGCNA_LASYseq.datExpr.n_vars}")
print(f"データ形状: {pyWGCNA_LASYseq.datExpr.shape}")

# %% [markdown]
# ## Run pyWGCNA

# %%
pyWGCNA_LASYseq.findModules()

# %%
pyWGCNA_LASYseq.analyseWGCNA(order=['species', 'age'])

# %%
import pandas as pd

module_list = pyWGCNA_LASYseq.datExpr.var['moduleColors'].unique()

print(module_list)

list_df_module_tmp = []

for module in module_list:
    
    module_genes = pyWGCNA_LASYseq.datExpr.var[
        pyWGCNA_LASYseq.datExpr.var['moduleColors'] == module
    ].copy()

    df_module_tmp = pyWGCNA_LASYseq.top_n_hub_genes(moduleName=module, n=len(module_genes) + 1)

    list_df_module_tmp.append(df_module_tmp)

df_module_all_genes = pd.concat(list_df_module_tmp, axis=0)

df_module_all_genes.to_csv('/Users/machiinagatoshi/Desktop/Research_Desk/RNAseq/LASY-seq/wgcna/all_genes_with_module.csv')

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


def plot_custom_module_trait_heatmap(wgcna_obj, selected_modules, selected_traits, 
                                   rename_dict, figsize=(12, 8), alternative='two-sided', 
                                   show=True, save_path=None):
    """
    カスタムWGCNAモジュール-表現型関係ヒートマップを描画する関数
    
    Parameters:
    -----------
    wgcna_obj : pyWGCNA_LASYseq object
    selected_modules : list
    selected_traits : list
    figsize : tuple, optional
    alternative : str, optional
    show : bool, optional
    save_path : str, optional
    
    Returns:
    --------
    tuple
        (correlation_matrix, pvalue_matrix, fig)
    """
    
    available_traits = wgcna_obj.datExpr.obs.columns.tolist()
    valid_traits = [trait for trait in selected_traits if trait in available_traits]
    
    if len(valid_traits) != len(selected_traits):
        missing_traits = [trait for trait in selected_traits if trait not in available_traits]
        print(f"警告: 以下の表現型がdatExpr.obsに存在しません: {missing_traits}")
    
    if len(valid_traits) == 0:
        raise ValueError("有効な表現型が見つかりません")
    
    datTraits = wgcna_obj.getDatTraits(valid_traits)
    print(f"ダミー変数に変換された表現型: {list(datTraits.columns)}")
    
    actual_trait_columns = list(datTraits.columns)
    
    me_columns = ['ME' + module for module in selected_modules]
    
    available_me_columns = [col for col in me_columns if col in wgcna_obj.MEs.columns]
    if len(available_me_columns) != len(me_columns):
        missing_modules = [col.replace('ME', '') for col in me_columns if col not in wgcna_obj.MEs.columns]
        print(f"警告: 以下のモジュールが見つかりません: {missing_modules}")
        me_columns = available_me_columns
        selected_modules = [col.replace('ME', '') for col in available_me_columns]
    
    moduleTraitCor = pd.DataFrame(index=me_columns, columns=actual_trait_columns, dtype="float")
    moduleTraitPvalue = pd.DataFrame(index=me_columns, columns=actual_trait_columns, dtype="float")
    
    for me_col in me_columns:
        for trait in actual_trait_columns:
            tmp = stats.pearsonr(wgcna_obj.MEs[me_col], datTraits[trait], alternative=alternative)
            moduleTraitCor.loc[me_col, trait] = tmp[0]
            moduleTraitPvalue.loc[me_col, trait] = tmp[1]
    
    # make fig
    
    xlabels = []
    for module in selected_modules:
        if 'moduleColors' in wgcna_obj.datExpr.var.columns:
            gene_count = sum(wgcna_obj.datExpr.var['moduleColors'] == module)
            rename_label = rename_dict[module]
            xlabels.append(f"{rename_label} ({gene_count})")
        else:
            xlabels.append(module.capitalize())
    
    ylabels = actual_trait_columns
    
    tmp_cor = moduleTraitCor.T.round(decimals=2)
    tmp_pvalue = moduleTraitPvalue.T.round(decimals=3)
    
    labels = (np.asarray(["{0}\n({1})".format(cor, pvalue) 
                         for cor, pvalue in zip(tmp_cor.values.flatten(), 
                                              tmp_pvalue.values.flatten())])).\
             reshape(moduleTraitCor.T.shape)
    
    sns.set(font_scale=1.2)
    res = sns.clustermap(moduleTraitCor.T,
                         row_cluster=False,
                         annot=labels, 
                         fmt="", 
                         cmap='RdBu_r', 
                         vmin=-1, 
                         vmax=1, 
                         annot_kws={'size': 10, "weight": "bold"},
                         xticklabels=xlabels,
                         yticklabels=ylabels,
                         figsize=figsize,
                         dendrogram_ratio=0.2,
                         cbar_kws={'label': 'Correlation'})    
    
    plt.tight_layout()
    
    if save_path:
        plt.rcParams['svg.fonttype'] = 'none'
        res.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"図を保存しました: {save_path}")
    
    if show:
        plt.show()
    
    return moduleTraitCor, moduleTraitPvalue, res


selected_modules = ['black', 'mistyrose', 'brown', 'darkgrey', 'snow', 'dimgrey', 'darksalmon',
                    'seashell', 'orangered', 'red', 'darkred', 'lightgrey', 'salmon', 'white',
                    'lightcoral', 'sandybrown', 'coral', 'rosybrown', 'firebrick', 'gainsboro',
                    'indianred', 'lightsalmon', 'chocolate', 'saddlebrown', 'peachpuff']

rename_dict =  {'black':'M1', 'seashell':'M2', 'darksalmon':'M3', 'lightgrey':'M4', 'dimgrey':'M5', 'gainsboro':'M6',
                'lightcoral':'M7', 'brown':'M8', 'darkgrey':'M9', 'white':'M10', 'mistyrose':'M11', 'snow':'M12', 'darkred':'M13',
                'orangered':'M14', 'salmon':'M15', 'indianred':'M16', 'sandybrown':'M17', 'rosybrown':'M18', 'coral':'M19',
                'red':'M20', 'firebrick':'M21', 'saddlebrown':'M22', 'lightsalmon':'M23', 'chocolate':'M24', 'peachpuff':'M25'}

selected_traits = ['species', 'species_age']

correlation_matrix, pvalue_matrix, fig = plot_custom_module_trait_heatmap(
    wgcna_obj=pyWGCNA_LASYseq,
    selected_modules=selected_modules,
    selected_traits=selected_traits,
    rename_dict=rename_dict,
    figsize=(20, 8),
    save_path='custom_heatmap.svg'
)


