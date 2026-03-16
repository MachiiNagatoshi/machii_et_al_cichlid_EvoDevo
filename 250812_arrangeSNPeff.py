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
# # 250812 analyze SNPeff
#
# reference
# - https://kazumaxneo.hatenablog.com/entry/2021/05/13/120000
# - https://pcingola.github.io/SnpEff/

# %% [markdown]
#
# | Column name | Meaning |
# |-------------|---------|
# | GeneName | Gene name (usually HUGO) |
# | GeneId | Gene's ID |
# | TranscriptId | Transcript's ID |
# | BioType | Transcript's bio-type (if available) |
# | variants_impact_* | Count number of variants for each impact category (HIGH, MODERATE, LOW, MODIFIER) |
# | variants_effect_* | Count number of variants for each effect type (e.g. missense_variant, synonymous_variant, stop_lost, etc.) |
# | bases_affected_* | Number of bases that variants overlap genomic region (repeated for several genomic regions: DOWNSTREAM, EXON, INTRON, UPSTREAM, etc.) |
# | total_score_* | Sum of scores overlapping this genomic region. Note: Scores are only available when input files are type 'BED' (e.g. when annotating ChipSeq experiments) |
# | length_* | Genomic region length |
#
#
# | Impact | Meaning | Example |
# |--------|---------|---------|
# | **HIGH** | The variant is assumed to have high (disruptive) impact in the protein, probably causing protein truncation, loss of function or triggering nonsense mediated decay. | `stop_gained`, `frameshift_variant` |
# | **MODERATE** | A non-disruptive variant that might change protein effectiveness. | `missense_variant`, `inframe_deletion` |
# | **LOW** | Assumed to be mostly harmless or unlikely to change protein behavior. | `synonymous_variant` |
# | **MODIFIER** | Usually non-coding variants or variants affecting non-coding genes, where predictions are difficult or there is no evidence of impact. | `exon_variant`, `downstream_gene_variant` |

# %% [markdown]
# ## Install and run SNPeff

# %% [markdown]
# ```bash
# (base) ~ ❯❯❯ conda create -n snpeff -y                                                          ✘ 1 
# Channels:
#  - conda-forge
# Platform: osx-64
# Collecting package metadata (repodata.json): done
# Solving environment: done
#
#
# ==> WARNING: A newer version of conda exists. <==
#     current version: 24.11.3
#     latest version: 25.7.0
#
# Please update conda by running
#
#     $ conda update -n base -c conda-forge conda
#
#
#
# ## Package Plan ##
#
#   environment location: /Users/machiinagatoshi/mambaforge/envs/snpeff
#
#
#
#
# Downloading and Extracting Packages:
#
# Preparing transaction: done
# Verifying transaction: done
# Executing transaction: done
# #
# # To activate this environment, use
# #
# #     $ conda activate snpeff
# #
# # To deactivate an active environment, use
# #
# #     $ conda deactivate
#
# (base) ~ ❯❯❯ conda activate snpeff
# (snpeff) ~ ❯❯❯ conda install -c bioconda -y snpeff snpsift
# Channels:
#  - bioconda
#  - conda-forge
# Platform: osx-64
# Collecting package metadata (repodata.json): done
# Solving environment: done
#
#
# ==> WARNING: A newer version of conda exists. <==
#     current version: 24.11.3
#     latest version: 25.7.0
#
# Please update conda by running
#
#     $ conda update -n base -c conda-forge conda
#
#
#
# ## Package Plan ##
#
#   environment location: /Users/machiinagatoshi/mambaforge/envs/snpeff
#
#   added / updated specs:
#     - snpeff
#     - snpsift
#
#
# The following packages will be downloaded:
#
#     package                    |            build
#     ---------------------------|-----------------
#     libexpat-2.7.1             |       h21dd04a_0          71 KB  conda-forge
#     libmpdec-4.0.0             |       h6e16a3a_0          76 KB  conda-forge
#     openjdk-23.0.2             |       h18c9476_2       176.7 MB  conda-forge
#     pip-25.2                   |     pyh145f28c_0         1.1 MB  conda-forge
#     python-3.13.5              |hc3a4c56_102_cp313        13.3 MB  conda-forge
#     python_abi-3.13            |          8_cp313           7 KB  conda-forge
#     ------------------------------------------------------------
#                                            Total:       191.3 MB
#
# The following NEW packages will be INSTALLED:
#
#   bzip2              conda-forge/osx-64::bzip2-1.0.8-hfdf4475_7 
#   ca-certificates    conda-forge/noarch::ca-certificates-2025.8.3-hbd8a1cb_0 
#   libexpat           conda-forge/osx-64::libexpat-2.7.1-h21dd04a_0 
#   libffi             conda-forge/osx-64::libffi-3.4.6-h281671d_1 
#   liblzma            conda-forge/osx-64::liblzma-5.8.1-hd471939_2 
#   libmpdec           conda-forge/osx-64::libmpdec-4.0.0-h6e16a3a_0 
#   libsqlite          conda-forge/osx-64::libsqlite-3.50.4-h39a8b3b_0 
#   libzlib            conda-forge/osx-64::libzlib-1.3.1-hd23fc13_2 
#   ncurses            conda-forge/osx-64::ncurses-6.5-h0622a9a_3 
#   openjdk            conda-forge/osx-64::openjdk-23.0.2-h18c9476_2 
#   openssl            conda-forge/osx-64::openssl-3.5.2-h6e31bce_0 
#   perl               conda-forge/osx-64::perl-5.32.1-7_h10d778d_perl5 
#   pip                conda-forge/noarch::pip-25.2-pyh145f28c_0 
#   python             conda-forge/osx-64::python-3.13.5-hc3a4c56_102_cp313 
#   python_abi         conda-forge/noarch::python_abi-3.13-8_cp313 
#   readline           conda-forge/osx-64::readline-8.2-h7cca4af_2 
#   snpeff             bioconda/noarch::snpeff-5.2-hdfd78af_1 
#   snpsift            bioconda/noarch::snpsift-5.2-hdfd78af_0 
#   tk                 conda-forge/osx-64::tk-8.6.13-hf689a15_2 
#   tzdata             conda-forge/noarch::tzdata-2025b-h78e105d_0 
#   zlib               conda-forge/osx-64::zlib-1.3.1-hd23fc13_2 
#
#
#
# Downloading and Extracting Packages:
#                                                                                                      
# Preparing transaction: done                                                                          
# Verifying transaction: done                                                                          
# Executing transaction: done                                                                          
# (snpeff) ~ ❯❯❯ conda list
# # packages in environment at /Users/machiinagatoshi/mambaforge/envs/snpeff:
# #
# # Name                    Version                   Build  Channel
# bzip2                     1.0.8                hfdf4475_7    conda-forge
# ca-certificates           2025.8.3             hbd8a1cb_0    conda-forge
# libexpat                  2.7.1                h21dd04a_0    conda-forge
# libffi                    3.4.6                h281671d_1    conda-forge
# liblzma                   5.8.1                hd471939_2    conda-forge
# libmpdec                  4.0.0                h6e16a3a_0    conda-forge
# libsqlite                 3.50.4               h39a8b3b_0    conda-forge
# libzlib                   1.3.1                hd23fc13_2    conda-forge
# ncurses                   6.5                  h0622a9a_3    conda-forge
# openjdk                   23.0.2               h18c9476_2    conda-forge
# openssl                   3.5.2                h6e31bce_0    conda-forge
# perl                      5.32.1          7_h10d778d_perl5    conda-forge
# pip                       25.2               pyh145f28c_0    conda-forge
# python                    3.13.5          hc3a4c56_102_cp313    conda-forge
# python_abi                3.13                    8_cp313    conda-forge
# readline                  8.2                  h7cca4af_2    conda-forge
# snpeff                    5.2                  hdfd78af_1    bioconda
# snpsift                   5.2                  hdfd78af_0    bioconda
# tk                        8.6.13               hf689a15_2    conda-forge
# tzdata                    2025b                h78e105d_0    conda-forge
# zlib                      1.3.1                hd23fc13_2    conda-forge
# (snpeff) ~ ❯❯❯ conda install pandas matplotlib seaborn numpy -c conda-forge
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
#     latest version: 25.7.0
#
# Please update conda by running
#
#     $ conda update -n base -c conda-forge conda
#
#
#
# ## Package Plan ##
#
#   environment location: /Users/machiinagatoshi/mambaforge/envs/snpeff
#
#   added / updated specs:
#     - matplotlib
#     - numpy
#     - pandas
#     - seaborn
#
#
# The following packages will be downloaded:
#
#     package                    |            build
#     ---------------------------|-----------------
#     brotli-1.1.0               |       h6e16a3a_3          20 KB  conda-forge
#     brotli-bin-1.1.0           |       h6e16a3a_3          17 KB  conda-forge
#     contourpy-1.3.3            |  py313hc551f4f_1         263 KB  conda-forge
#     fonttools-4.59.0           |  py313h4db2fa4_0         2.7 MB  conda-forge
#     kiwisolver-1.4.8           |  py313ha0b1807_1          62 KB  conda-forge
#     libblas-3.9.0              |33_h7f60823_openblas          19 KB  conda-forge
#     libbrotlicommon-1.1.0      |       h6e16a3a_3          66 KB  conda-forge
#     libbrotlidec-1.1.0         |       h6e16a3a_3          30 KB  conda-forge
#     libbrotlienc-1.1.0         |       h6e16a3a_3         288 KB  conda-forge
#     libcblas-3.9.0             |33_hff6cab4_openblas          19 KB  conda-forge
#     libcxx-20.1.8              |       h3d58e20_1         552 KB  conda-forge
#     libdeflate-1.24            |       hcc1b750_0          68 KB  conda-forge
#     liblapack-3.9.0            |33_h236ab99_openblas          19 KB  conda-forge
#     libopenblas-0.3.30         |openmp_hbf64a52_0         5.9 MB  conda-forge
#     libpng-1.6.50              |       h84aeda2_1         291 KB  conda-forge
#     libtiff-4.7.0              |       h1167cee_5         391 KB  conda-forge
#     libwebp-base-1.6.0         |       hb807250_0         357 KB  conda-forge
#     matplotlib-3.10.5          |  py313habf4b1d_0          17 KB  conda-forge
#     matplotlib-base-3.10.5     |  py313h5771d13_0         7.9 MB  conda-forge
#     munkres-1.1.4              |     pyhd8ed1ab_1          15 KB  conda-forge
#     numpy-2.3.2                |  py313hdb1a8e5_0         7.7 MB  conda-forge
#     openjpeg-2.5.3             |       h036ada5_1         328 KB  conda-forge
#     pandas-2.3.1               |  py313h366a99e_0        13.6 MB  conda-forge
#     patsy-1.0.1                |     pyhd8ed1ab_1         182 KB  conda-forge
#     pillow-11.3.0              |  py313h0c4f865_0        40.3 MB  conda-forge
#     pyparsing-3.2.3            |     pyhe01879c_2         100 KB  conda-forge
#     python-dateutil-2.9.0.post0|     pyhe01879c_2         228 KB  conda-forge
#     python-tzdata-2025.2       |     pyhd8ed1ab_0         141 KB  conda-forge
#     pytz-2025.2                |     pyhd8ed1ab_0         185 KB  conda-forge
#     scipy-1.16.0               |  py313h7e69c36_0        14.6 MB  conda-forge
#     seaborn-0.13.2             |       hd8ed1ab_3           7 KB  conda-forge
#     seaborn-base-0.13.2        |     pyhd8ed1ab_3         223 KB  conda-forge
#     six-1.17.0                 |     pyhe01879c_1          18 KB  conda-forge
#     statsmodels-0.14.5         |  py313h1492807_0        11.3 MB  conda-forge
#     tornado-6.5.1              |  py313h63b0ddb_0         854 KB  conda-forge
#     ------------------------------------------------------------
#                                            Total:       108.6 MB
#
# The following NEW packages will be INSTALLED:
#
#   brotli             conda-forge/osx-64::brotli-1.1.0-h6e16a3a_3 
#   brotli-bin         conda-forge/osx-64::brotli-bin-1.1.0-h6e16a3a_3 
#   contourpy          conda-forge/osx-64::contourpy-1.3.3-py313hc551f4f_1 
#   cycler             conda-forge/noarch::cycler-0.12.1-pyhd8ed1ab_1 
#   fonttools          conda-forge/osx-64::fonttools-4.59.0-py313h4db2fa4_0 
#   freetype           conda-forge/osx-64::freetype-2.13.3-h694c41f_1 
#   kiwisolver         conda-forge/osx-64::kiwisolver-1.4.8-py313ha0b1807_1 
#   lcms2              conda-forge/osx-64::lcms2-2.17-h72f5680_0 
#   lerc               conda-forge/osx-64::lerc-4.0.0-hcca01a6_1 
#   libblas            conda-forge/osx-64::libblas-3.9.0-33_h7f60823_openblas 
#   libbrotlicommon    conda-forge/osx-64::libbrotlicommon-1.1.0-h6e16a3a_3 
#   libbrotlidec       conda-forge/osx-64::libbrotlidec-1.1.0-h6e16a3a_3 
#   libbrotlienc       conda-forge/osx-64::libbrotlienc-1.1.0-h6e16a3a_3 
#   libcblas           conda-forge/osx-64::libcblas-3.9.0-33_hff6cab4_openblas 
#   libcxx             conda-forge/osx-64::libcxx-20.1.8-h3d58e20_1 
#   libdeflate         conda-forge/osx-64::libdeflate-1.24-hcc1b750_0 
#   libfreetype        conda-forge/osx-64::libfreetype-2.13.3-h694c41f_1 
#   libfreetype6       conda-forge/osx-64::libfreetype6-2.13.3-h40dfd5c_1 
#   libgfortran        conda-forge/osx-64::libgfortran-5.0.0-14_2_0_h51e75f0_103 
#   libgfortran5       conda-forge/osx-64::libgfortran5-14.2.0-h51e75f0_103 
#   libjpeg-turbo      conda-forge/osx-64::libjpeg-turbo-3.1.0-h6e16a3a_0 
#   liblapack          conda-forge/osx-64::liblapack-3.9.0-33_h236ab99_openblas 
#   libopenblas        conda-forge/osx-64::libopenblas-0.3.30-openmp_hbf64a52_0 
#   libpng             conda-forge/osx-64::libpng-1.6.50-h84aeda2_1 
#   libtiff            conda-forge/osx-64::libtiff-4.7.0-h1167cee_5 
#   libwebp-base       conda-forge/osx-64::libwebp-base-1.6.0-hb807250_0 
#   libxcb             conda-forge/osx-64::libxcb-1.17.0-hf1f96e2_0 
#   llvm-openmp        conda-forge/osx-64::llvm-openmp-20.1.8-hf4e0ed4_1 
#   matplotlib         conda-forge/osx-64::matplotlib-3.10.5-py313habf4b1d_0 
#   matplotlib-base    conda-forge/osx-64::matplotlib-base-3.10.5-py313h5771d13_0 
#   munkres            conda-forge/noarch::munkres-1.1.4-pyhd8ed1ab_1 
#   numpy              conda-forge/osx-64::numpy-2.3.2-py313hdb1a8e5_0 
#   openjpeg           conda-forge/osx-64::openjpeg-2.5.3-h036ada5_1 
#   packaging          conda-forge/noarch::packaging-25.0-pyh29332c3_1 
#   pandas             conda-forge/osx-64::pandas-2.3.1-py313h366a99e_0 
#   patsy              conda-forge/noarch::patsy-1.0.1-pyhd8ed1ab_1 
#   pillow             conda-forge/osx-64::pillow-11.3.0-py313h0c4f865_0 
#   pthread-stubs      conda-forge/osx-64::pthread-stubs-0.4-h00291cd_1002 
#   pyparsing          conda-forge/noarch::pyparsing-3.2.3-pyhe01879c_2 
#   python-dateutil    conda-forge/noarch::python-dateutil-2.9.0.post0-pyhe01879c_2 
#   python-tzdata      conda-forge/noarch::python-tzdata-2025.2-pyhd8ed1ab_0 
#   pytz               conda-forge/noarch::pytz-2025.2-pyhd8ed1ab_0 
#   qhull              conda-forge/osx-64::qhull-2020.2-h3c5361c_5 
#   scipy              conda-forge/osx-64::scipy-1.16.0-py313h7e69c36_0 
#   seaborn            conda-forge/noarch::seaborn-0.13.2-hd8ed1ab_3 
#   seaborn-base       conda-forge/noarch::seaborn-base-0.13.2-pyhd8ed1ab_3 
#   six                conda-forge/noarch::six-1.17.0-pyhe01879c_1 
#   statsmodels        conda-forge/osx-64::statsmodels-0.14.5-py313h1492807_0 
#   tornado            conda-forge/osx-64::tornado-6.5.1-py313h63b0ddb_0 
#   xorg-libxau        conda-forge/osx-64::xorg-libxau-1.0.12-h6e16a3a_0 
#   xorg-libxdmcp      conda-forge/osx-64::xorg-libxdmcp-1.1.5-h00291cd_0 
#   zstd               conda-forge/osx-64::zstd-1.5.7-h8210216_2 
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
# (snpeff) ~ ❯❯❯ snpEff databases
# Genome                                                      	Organism                                                    	Status    	Bundle                        	Database download link
# ------                                                      	--------                                                    	------    	------                        	----------------------
# 129S1_SvImJ_v1.105                                          	Mus_musculus_129s1svimj                                     	          	                              	[https://snpeff.blob.core.windows.net/databases/v5_2/snpEff_v5_2_129S1_SvImJ_v1.105.zip, https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_129S1_SvImJ_v1.105.zip, https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_129S1_SvImJ_v1.105.zip]
#
# (... データベースの一覧が表示される 以下にあるMaylandia_zebraのデータベースをダウンロード...)
#
# M_zebra_UMD2a.105                                           	Maylandia_zebra                                             	          	                              	[https://snpeff.blob.core.windows.net/databases/v5_2/snpEff_v5_2_M_zebra_UMD2a.105.zip, https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_M_zebra_UMD2a.105.zip, https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_M_zebra_UMD2a.105.zip]
# M_zebra_UMD2a.99                                            	Maylandia_zebra                                             	          	                              	[https://snpeff.blob.core.windows.net/databases/v5_2/snpEff_v5_2_M_zebra_UMD2a.99.zip, https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_M_zebra_UMD2a.99.zip, https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_M_zebra_UMD2a.99.zip]
#
# (snpeff) ~ ❯❯❯ snpEff download -v M_zebra_UMD2a.105
# 00:00:00 SnpEff version SnpEff 5.2 (build 2023-09-29 06:17), by Pablo Cingolani
# 00:00:00 Command: 'download'
# 00:00:00 Reading configuration file 'snpEff.config'. Genome: 'M_zebra_UMD2a.105'
# 00:00:00 Reading config file: /Users/machiinagatoshi/snpEff.config
# 00:00:00 Reading config file: /Users/machiinagatoshi/mambaforge/envs/snpeff/share/snpeff-5.2-1/snpEff.config
# 00:00:02 done
# 00:00:02 Downloading database for 'M_zebra_UMD2a.105'
# 00:00:02 Downloading from 'https://snpeff.blob.core.windows.net/databases/v5_2/snpEff_v5_2_M_zebra_UMD2a.105.zip' to local file '/var/folders/cy/fg1wkjzd1fq6mcb6y258n28m0000gn/T//snpEff_v5_2_M_zebra_UMD2a.105.zip'
# 00:00:02 Connecting to https://snpeff.blob.core.windows.net/databases/v5_2/snpEff_v5_2_M_zebra_UMD2a.105.zip
# 00:00:02 Connecting to https://snpeff.blob.core.windows.net/databases/v5_2/snpEff_v5_2_M_zebra_UMD2a.105.zip, using proxy: false
# 00:00:03 ERROR while connecting to https://snpeff.blob.core.windows.net/databases/v5_2/snpEff_v5_2_M_zebra_UMD2a.105.zip
# 00:00:03 Downloading from 'https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_M_zebra_UMD2a.105.zip' to local file '/var/folders/cy/fg1wkjzd1fq6mcb6y258n28m0000gn/T//snpEff_v5_0_M_zebra_UMD2a.105.zip'
# 00:00:03 Connecting to https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_M_zebra_UMD2a.105.zip
# 00:00:03 Connecting to https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_M_zebra_UMD2a.105.zip, using proxy: false
# 00:00:04 ERROR while connecting to https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_M_zebra_UMD2a.105.zip
# 00:00:04 Downloading from 'https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_M_zebra_UMD2a.105.zip' to local file '/var/folders/cy/fg1wkjzd1fq6mcb6y258n28m0000gn/T//snpEff_v5_1_M_zebra_UMD2a.105.zip'
# 00:00:04 Connecting to https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_M_zebra_UMD2a.105.zip
# 00:00:04 Connecting to https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_M_zebra_UMD2a.105.zip, using proxy: false
# 00:00:04 Local file name: '/var/folders/cy/fg1wkjzd1fq6mcb6y258n28m0000gn/T//snpEff_v5_1_M_zebra_UMD2a.105.zip'
# .............................................................................................................................................................00:00:20 
# 00:00:20 Download finished. Total 167685348 bytes.
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.AGTA05000023.1.bin'
# 00:00:20 Creating local directory: '/Users/machiinagatoshi/mambaforge/envs/snpeff/share/snpeff-5.2-1/./data/M_zebra_UMD2a.105'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG1.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG10.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG11.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG12.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG13.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG14.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG15.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG16.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG17.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG18.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG19.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG2.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG20.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG22.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG23.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG3.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG4.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG5.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG6.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG7.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG8.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.LG9.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/sequence.bin'
# 00:00:20 Extracting file 'data/M_zebra_UMD2a.105/snpEffectPredictor.bin'
# 00:00:20 Unzip: OK
# 00:00:20 Deleted local file '/var/folders/cy/fg1wkjzd1fq6mcb6y258n28m0000gn/T//snpEff_v5_1_M_zebra_UMD2a.105.zip'
# 00:00:20 Done
# 00:00:20 Done.
# (snpeff) ~ ❯❯❯ snpEff -v -stats report.html M_zebra_UMD2a.105 /Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset1_missing0.9_maf0.05.vcf.recode.vcf.gz > vcf_dataset1_missing0.9_maf0.05.vcf.recode.annotated.vcf
# 00:00:00 SnpEff version SnpEff 5.2 (build 2023-09-29 06:17), by Pablo Cingolani
# 00:00:00 Command: 'ann'
# 00:00:00 Reading configuration file 'snpEff.config'. Genome: 'M_zebra_UMD2a.105'
# 00:00:00 Reading config file: /Users/machiinagatoshi/snpEff.config
# 00:00:00 Reading config file: /Users/machiinagatoshi/mambaforge/envs/snpeff/share/snpeff-5.2-1/snpEff.config
# 00:00:02 done
# 00:00:02 Reading database for genome version 'M_zebra_UMD2a.105' from file '/Users/machiinagatoshi/mambaforge/envs/snpeff/share/snpeff-5.2-1/./data/M_zebra_UMD2a.105/snpEffectPredictor.bin' (this might take a while)
# 00:00:05 done
# 00:00:05 Loading Motifs and PWMs
# 00:00:05 Building interval forest
# 00:00:07 done.
# 00:00:07 Genome stats :
# #-----------------------------------------------
# # Genome name                : 'Maylandia_zebra'
# # Genome version             : 'M_zebra_UMD2a.105'
# # Genome ID                  : 'M_zebra_UMD2a.105[0]'
# # Has protein coding info    : true
# # Has Tr. Support Level info : true
# # Genes                      : 28622
# # Protein coding genes       : 27187
# #-----------------------------------------------
# # Transcripts                : 39681
# # Avg. transcripts per gene  : 1.39
# # TSL transcripts            : 0
# #-----------------------------------------------
# # Checked transcripts        : 
# #               AA sequences :  38233 ( 99.97% )
# #              DNA sequences :  38295 ( 96.51% )
# #-----------------------------------------------
# # Protein coding transcripts : 38246
# #              Length errors :     58 ( 0.15% )
# #  STOP codons in CDS errors :     12 ( 0.03% )
# #         START codon errors :   5769 ( 15.08% )
# #        STOP codon warnings :    547 ( 1.43% )
# #              UTR sequences :  19303 ( 48.65% )
# #               Total Errors :   5810 ( 15.19% )
# #-----------------------------------------------
# # Cds                        : 375491
# # Exons                      : 384819
# # Exons with sequence        : 384819
# # Exons without sequence     : 0
# # Avg. exons per transcript  : 9.70
# # WARNING!                   : Mitochondrion chromosome 'MT' does not have a mitochondrion codon table (codon table = 'Standard'). You should update the config file.
# #-----------------------------------------------
# # Number of chromosomes      : 1690
# # Chromosomes                : Format 'chromo_name size codon_table'
# #		'LG7'	64916660	Standard
# #		'LG23'	42088218	Standard
# #		'LG6'	39774503	Standard
# #		'LG1'	38676823	Standard
# #		'LG14'	37870038	Standard
# #		'LG3'	37314939	Standard
# #		'LG5'	36170306	Standard
# #		'LG17'	35776103	Standard
# #		'LG16'	34742138	Standard
# #		'LG22'	34725690	Standard
# #		'LG15'	34548429	Standard
# #		'LG12'	34086040	Standard
# #		'LG2'	32660920	Standard
# #		'LG11'	32446080	Standard
# #		'LG10'	32356376	Standard
# #		'LG13'	32072427	Standard
# #		'LG4'	30518969	Standard
# #		'LG20'	29789237	Standard
# #		'LG18'	29505383	Standard
# #		'LG19'	25963618	Standard
# #		'LG8'	23971387	Standard
# #		'LG9'	21023328	Standard
# #		'AGTA05000023.1'	3684755	Standard
# #		'AGTA05000024.1'	1539574	Standard
# #		'AGTA05000025.1'	1394719	Standard
# #		'AGTA05000026.1'	1368916	Standard
# #		'AGTA05000027.1'	1339124	Standard
# #		'AGTA05000028.1'	1268605	Standard
# #		'AGTA05000029.1'	1242648	Standard
# #		'AGTA05000030.1'	1155072	Standard
# #		'AGTA05000031.1'	1146798	Standard
# #		'AGTA05000032.1'	1098226	Standard
# #		'AGTA05000033.1'	989703	Standard
# #		'AGTA05000034.1'	978946	Standard
# #		'AGTA05000035.1'	976411	Standard
# #		'AGTA05000036.1'	975679	Standard
# #		'AGTA05000037.1'	954468	Standard
# #		'AGTA05000038.1'	936231	Standard
# #		'AGTA05000039.1'	933308	Standard
# #		'AGTA05000040.1'	929013	Standard
# #		'AGTA05000041.1'	926824	Standard
# #		'AGTA05000042.1'	913098	Standard
# #		'AGTA05000043.1'	862263	Standard
# #		'AGTA05000044.1'	840938	Standard
# #		'AGTA05000045.1'	806055	Standard
# #		'AGTA05000046.1'	784528	Standard
# #		'AGTA05000047.1'	708617	Standard
# #		'AGTA05000048.1'	695985	Standard
# #		'AGTA05000049.1'	687606	Standard
# #		'AGTA05000050.1'	663057	Standard
# #		'AGTA05000051.1'	653464	Standard
# #		'AGTA05000052.1'	649919	Standard
# #		'AGTA05000053.1'	646767	Standard
# #		'AGTA05000054.1'	632296	Standard
# #		'AGTA05000055.1'	631966	Standard
# #		'AGTA05000056.1'	631000	Standard
# #		'AGTA05000057.1'	625576	Standard
# #		'AGTA05000058.1'	607906	Standard
# #		'AGTA05000059.1'	605074	Standard
# #		'AGTA05000060.1'	604006	Standard
# #		'AGTA05000061.1'	603684	Standard
# #		'AGTA05000062.1'	595956	Standard
# #		'AGTA05000063.1'	592058	Standard
# #		'AGTA05000064.1'	591508	Standard
# #		'AGTA05000065.1'	578447	Standard
# #		'AGTA05000066.1'	575532	Standard
# #		'AGTA05000067.1'	568254	Standard
# #		'AGTA05000068.1'	550753	Standard
# #		'AGTA05000069.1'	545043	Standard
# #		'AGTA05000070.1'	543777	Standard
# #		'AGTA05000071.1'	529070	Standard
# #		'AGTA05000072.1'	525310	Standard
# #		'AGTA05000073.1'	524065	Standard
# #		'AGTA05000074.1'	519457	Standard
# #		'AGTA05000075.1'	509617	Standard
# #		'AGTA05000076.1'	506933	Standard
# #		'AGTA05000077.1'	502303	Standard
# #		'AGTA05000078.1'	499129	Standard
# #		'AGTA05000079.1'	497711	Standard
# #		'AGTA05000080.1'	485592	Standard
# #		'AGTA05000081.1'	483123	Standard
# #		'AGTA05000082.1'	477780	Standard
# #		'AGTA05000083.1'	475462	Standard
# #		'AGTA05000085.1'	474539	Standard
# #		'AGTA05000084.1'	474412	Standard
# #		'AGTA05000086.1'	471871	Standard
# #		'AGTA05000087.1'	469127	Standard
# #		'AGTA05000088.1'	467497	Standard
# #		'AGTA05000089.1'	460424	Standard
# #		'AGTA05000090.1'	459060	Standard
# #		'AGTA05000091.1'	457383	Standard
# #		'AGTA05000092.1'	456270	Standard
# #		'AGTA05000093.1'	452260	Standard
# #		'AGTA05000094.1'	451697	Standard
# #		'AGTA05000095.1'	448427	Standard
# #		'AGTA05000096.1'	447980	Standard
# #		'AGTA05000097.1'	444561	Standard
# #		'AGTA05000098.1'	441931	Standard
# #		'AGTA05000099.1'	440022	Standard
# #		'AGTA05000100.1'	436819	Standard
# #		'AGTA05000101.1'	423697	Standard
# #		'AGTA05000102.1'	423265	Standard
# #		'AGTA05000103.1'	421393	Standard
# #		'AGTA05000104.1'	412762	Standard
# #		'AGTA05000105.1'	405582	Standard
# #		'AGTA05000106.1'	400850	Standard
# #		'AGTA05000107.1'	400406	Standard
# #		'AGTA05000108.1'	398292	Standard
# #		'AGTA05000109.1'	395429	Standard
# #		'AGTA05000111.1'	393896	Standard
# #		'AGTA05000110.1'	393888	Standard
# #		'AGTA05000112.1'	390354	Standard
# #		'AGTA05000113.1'	389928	Standard
# #		'AGTA05000114.1'	386321	Standard
# #		'AGTA05000115.1'	385493	Standard
# #		'AGTA05000116.1'	384620	Standard
# #		'AGTA05000117.1'	383335	Standard
# #		'AGTA05000118.1'	381117	Standard
# #		'AGTA05000119.1'	377351	Standard
# #		'AGTA05000120.1'	376023	Standard
# #		'AGTA05000121.1'	375391	Standard
# #		'AGTA05000122.1'	374273	Standard
# #		'AGTA05000123.1'	372136	Standard
# #		'AGTA05000124.1'	371441	Standard
# #		'AGTA05000125.1'	369756	Standard
# #		'AGTA05000126.1'	368237	Standard
# #		'AGTA05000127.1'	368165	Standard
# #		'AGTA05000128.1'	366926	Standard
# #		'AGTA05000129.1'	364066	Standard
# #		'AGTA05000130.1'	359572	Standard
# #		'AGTA05000131.1'	356976	Standard
# #		'AGTA05000133.1'	353627	Standard
# #		'AGTA05000132.1'	353493	Standard
# #		'AGTA05000134.1'	349145	Standard
# #		'AGTA05000135.1'	348288	Standard
# #		'AGTA05000137.1'	347227	Standard
# #		'AGTA05000136.1'	347042	Standard
# #		'AGTA05000138.1'	345202	Standard
# #		'AGTA05000139.1'	340796	Standard
# #		'AGTA05000140.1'	340776	Standard
# #		'AGTA05000141.1'	339969	Standard
# #		'AGTA05000142.1'	330251	Standard
# #		'AGTA05000143.1'	327265	Standard
# #		'AGTA05000144.1'	323311	Standard
# #		'AGTA05000145.1'	317927	Standard
# #		'AGTA05000146.1'	317505	Standard
# #		'AGTA05000147.1'	315575	Standard
# #		'AGTA05000148.1'	313472	Standard
# #		'AGTA05000149.1'	312791	Standard
# #		'AGTA05000150.1'	311326	Standard
# #		'AGTA05000151.1'	309903	Standard
# #		'AGTA05000152.1'	308849	Standard
# #		'AGTA05000153.1'	308030	Standard
# #		'AGTA05000154.1'	307774	Standard
# #		'AGTA05000155.1'	307036	Standard
# #		'AGTA05000156.1'	305811	Standard
# #		'AGTA05000157.1'	305622	Standard
# #		'AGTA05000158.1'	303986	Standard
# #		'AGTA05000159.1'	302909	Standard
# #		'AGTA05000160.1'	299292	Standard
# #		'AGTA05000161.1'	298353	Standard
# #		'AGTA05000162.1'	296300	Standard
# #		'AGTA05000163.1'	294455	Standard
# #		'AGTA05000164.1'	293967	Standard
# #		'AGTA05000165.1'	293451	Standard
# #		'AGTA05000166.1'	289159	Standard
# #		'AGTA05000167.1'	288789	Standard
# #		'AGTA05000168.1'	288415	Standard
# #		'AGTA05000169.1'	288072	Standard
# #		'AGTA05000170.1'	286313	Standard
# #		'AGTA05000171.1'	284899	Standard
# #		'AGTA05000172.1'	283039	Standard
# #		'AGTA05000173.1'	281721	Standard
# #		'AGTA05000174.1'	280601	Standard
# #		'AGTA05000175.1'	279375	Standard
# #		'AGTA05000176.1'	279119	Standard
# #		'AGTA05000177.1'	279053	Standard
# #		'AGTA05000178.1'	278660	Standard
# #		'AGTA05000179.1'	275015	Standard
# #		'AGTA05000180.1'	274014	Standard
# #		'AGTA05000181.1'	273373	Standard
# #		'AGTA05000182.1'	271199	Standard
# #		'AGTA05000183.1'	266193	Standard
# #		'AGTA05000184.1'	265337	Standard
# #		'AGTA05000185.1'	264660	Standard
# #		'AGTA05000186.1'	264172	Standard
# #		'AGTA05000187.1'	264049	Standard
# #		'AGTA05000188.1'	263244	Standard
# #		'AGTA05000189.1'	262310	Standard
# #		'AGTA05000190.1'	262067	Standard
# #		'AGTA05000191.1'	261836	Standard
# #		'AGTA05000192.1'	260332	Standard
# #		'AGTA05000195.1'	257826	Standard
# #		'AGTA05000193.1'	257697	Standard
# #		'AGTA05000194.1'	257612	Standard
# #		'AGTA05000196.1'	254364	Standard
# #		'AGTA05000197.1'	251481	Standard
# #		'AGTA05000198.1'	251363	Standard
# #		'AGTA05000199.1'	248144	Standard
# #		'AGTA05000200.1'	247221	Standard
# #		'AGTA05000201.1'	244930	Standard
# #		'AGTA05000202.1'	244722	Standard
# #		'AGTA05000203.1'	244113	Standard
# #		'AGTA05000204.1'	243801	Standard
# #		'AGTA05000205.1'	242989	Standard
# #		'AGTA05000206.1'	241555	Standard
# #		'AGTA05000207.1'	238655	Standard
# #		'AGTA05000208.1'	238118	Standard
# #		'AGTA05000209.1'	237953	Standard
# #		'AGTA05000210.1'	236470	Standard
# #		'AGTA05000211.1'	235618	Standard
# #		'AGTA05000212.1'	234584	Standard
# #		'AGTA05000213.1'	234488	Standard
# #		'AGTA05000214.1'	233474	Standard
# #		'AGTA05000215.1'	231729	Standard
# #		'AGTA05000216.1'	231276	Standard
# #		'AGTA05000217.1'	230347	Standard
# #		'AGTA05000218.1'	230081	Standard
# #		'AGTA05000219.1'	229603	Standard
# #		'AGTA05000220.1'	227829	Standard
# #		'AGTA05000221.1'	223282	Standard
# #		'AGTA05000222.1'	222247	Standard
# #		'AGTA05000223.1'	221883	Standard
# #		'AGTA05000224.1'	220559	Standard
# #		'AGTA05000225.1'	220318	Standard
# #		'AGTA05000226.1'	220192	Standard
# #		'AGTA05000227.1'	219083	Standard
# #		'AGTA05000228.1'	218837	Standard
# #		'AGTA05000229.1'	217188	Standard
# #		'AGTA05000230.1'	216570	Standard
# #		'AGTA05000231.1'	215977	Standard
# #		'AGTA05000232.1'	214707	Standard
# #		'AGTA05000233.1'	212563	Standard
# #		'AGTA05000234.1'	211718	Standard
# #		'AGTA05000236.1'	211426	Standard
# #		'AGTA05000237.1'	211296	Standard
# #		'AGTA05000239.1'	211225	Standard
# #		'AGTA05000238.1'	211222	Standard
# #		'AGTA05000235.1'	211213	Standard
# #		'AGTA05000240.1'	209323	Standard
# #		'AGTA05000241.1'	208368	Standard
# #		'AGTA05000242.1'	205461	Standard
# #		'AGTA05000243.1'	204860	Standard
# #		'AGTA05000244.1'	203051	Standard
# #		'AGTA05000245.1'	203035	Standard
# #		'AGTA05000246.1'	202942	Standard
# #		'AGTA05000247.1'	202234	Standard
# #		'AGTA05000248.1'	201766	Standard
# #		'AGTA05000249.1'	201434	Standard
# #		'AGTA05000250.1'	201236	Standard
# #		'AGTA05000251.1'	199383	Standard
# #		'AGTA05000252.1'	199152	Standard
# #		'AGTA05000253.1'	197506	Standard
# #		'AGTA05000255.1'	197025	Standard
# #		'AGTA05000254.1'	196865	Standard
# #		'AGTA05000256.1'	195862	Standard
# #		'AGTA05000257.1'	195220	Standard
# #		'AGTA05000259.1'	194103	Standard
# #		'AGTA05000258.1'	193603	Standard
# #		'AGTA05000260.1'	193243	Standard
# #		'AGTA05000261.1'	192992	Standard
# #		'AGTA05000262.1'	191318	Standard
# #		'AGTA05000264.1'	190609	Standard
# #		'AGTA05000263.1'	190578	Standard
# #		'AGTA05000265.1'	190369	Standard
# #		'AGTA05000267.1'	190316	Standard
# #		'AGTA05000266.1'	190294	Standard
# #		'AGTA05000268.1'	189972	Standard
# #		'AGTA05000269.1'	188875	Standard
# #		'AGTA05000270.1'	188624	Standard
# #		'AGTA05000271.1'	188067	Standard
# #		'AGTA05000272.1'	187448	Standard
# #		'AGTA05000273.1'	187043	Standard
# #		'AGTA05000274.1'	186179	Standard
# #		'AGTA05000276.1'	186088	Standard
# #		'AGTA05000275.1'	186026	Standard
# #		'AGTA05000277.1'	185942	Standard
# #		'AGTA05000278.1'	185448	Standard
# #		'AGTA05000279.1'	185427	Standard
# #		'AGTA05000280.1'	185423	Standard
# #		'AGTA05000282.1'	184615	Standard
# #		'AGTA05000281.1'	184474	Standard
# #		'AGTA05000283.1'	184292	Standard
# #		'AGTA05000284.1'	183184	Standard
# #		'AGTA05000285.1'	183021	Standard
# #		'AGTA05000286.1'	180602	Standard
# #		'AGTA05000288.1'	179770	Standard
# #		'AGTA05000287.1'	179579	Standard
# #		'AGTA05000289.1'	178605	Standard
# #		'AGTA05000290.1'	177212	Standard
# #		'AGTA05000291.1'	176304	Standard
# #		'AGTA05000292.1'	175798	Standard
# #		'AGTA05000293.1'	175495	Standard
# #		'AGTA05000294.1'	174228	Standard
# #		'AGTA05000295.1'	173999	Standard
# #		'AGTA05000296.1'	173859	Standard
# #		'AGTA05000297.1'	173817	Standard
# #		'AGTA05000298.1'	173699	Standard
# #		'AGTA05000299.1'	173643	Standard
# #		'AGTA05000300.1'	173252	Standard
# #		'AGTA05000301.1'	172219	Standard
# #		'AGTA05000302.1'	171255	Standard
# #		'AGTA05000303.1'	170989	Standard
# #		'AGTA05000304.1'	170504	Standard
# #		'AGTA05000305.1'	169596	Standard
# #		'AGTA05000306.1'	169415	Standard
# #		'AGTA05000307.1'	169001	Standard
# #		'AGTA05000308.1'	167742	Standard
# #		'AGTA05000309.1'	166442	Standard
# #		'AGTA05000310.1'	166104	Standard
# #		'AGTA05000311.1'	165799	Standard
# #		'AGTA05000312.1'	165785	Standard
# #		'AGTA05000313.1'	165360	Standard
# #		'AGTA05000315.1'	165205	Standard
# #		'AGTA05000314.1'	165147	Standard
# #		'AGTA05000316.1'	163221	Standard
# #		'AGTA05000317.1'	162701	Standard
# #		'AGTA05000319.1'	162155	Standard
# #		'AGTA05000318.1'	162086	Standard
# #		'AGTA05000320.1'	161722	Standard
# #		'AGTA05000322.1'	160850	Standard
# #		'AGTA05000321.1'	160781	Standard
# #		'AGTA05000323.1'	160561	Standard
# #		'AGTA05000324.1'	160443	Standard
# #		'AGTA05000325.1'	159672	Standard
# #		'AGTA05000326.1'	159534	Standard
# #		'AGTA05000327.1'	157860	Standard
# #		'AGTA05000328.1'	157696	Standard
# #		'AGTA05000329.1'	157614	Standard
# #		'AGTA05000330.1'	156901	Standard
# #		'AGTA05000331.1'	156250	Standard
# #		'AGTA05000332.1'	155961	Standard
# #		'AGTA05000333.1'	155654	Standard
# #		'AGTA05000334.1'	155192	Standard
# #		'AGTA05000335.1'	154980	Standard
# #		'AGTA05000336.1'	154113	Standard
# #		'AGTA05000337.1'	153437	Standard
# #		'AGTA05000339.1'	152812	Standard
# #		'AGTA05000338.1'	152809	Standard
# #		'AGTA05000340.1'	152377	Standard
# #		'AGTA05000341.1'	152310	Standard
# #		'AGTA05000342.1'	152122	Standard
# #		'AGTA05000343.1'	152009	Standard
# #		'AGTA05000344.1'	151733	Standard
# #		'AGTA05000345.1'	151402	Standard
# #		'AGTA05000346.1'	151303	Standard
# #		'AGTA05000347.1'	151162	Standard
# #		'AGTA05000348.1'	149271	Standard
# #		'AGTA05000349.1'	148917	Standard
# #		'AGTA05000350.1'	148681	Standard
# #		'AGTA05000351.1'	148665	Standard
# #		'AGTA05000352.1'	147949	Standard
# #		'AGTA05000353.1'	147234	Standard
# #		'AGTA05000354.1'	146117	Standard
# #		'AGTA05000356.1'	144337	Standard
# #		'AGTA05000355.1'	144298	Standard
# #		'AGTA05000357.1'	143577	Standard
# #		'AGTA05000358.1'	143432	Standard
# #		'AGTA05000359.1'	143013	Standard
# #		'AGTA05000360.1'	142140	Standard
# #		'AGTA05000361.1'	141681	Standard
# #		'AGTA05000362.1'	140481	Standard
# #		'AGTA05000363.1'	140338	Standard
# #		'AGTA05000364.1'	140014	Standard
# #		'AGTA05000365.1'	139714	Standard
# #		'AGTA05000366.1'	139607	Standard
# #		'AGTA05000369.1'	139476	Standard
# #		'AGTA05000368.1'	139453	Standard
# #		'AGTA05000367.1'	139426	Standard
# #		'AGTA05000370.1'	138833	Standard
# #		'AGTA05000371.1'	138661	Standard
# #		'AGTA05000372.1'	137705	Standard
# #		'AGTA05000373.1'	137702	Standard
# #		'AGTA05000375.1'	137602	Standard
# #		'AGTA05000374.1'	137523	Standard
# #		'AGTA05000376.1'	137325	Standard
# #		'AGTA05000377.1'	136890	Standard
# #		'AGTA05000378.1'	136809	Standard
# #		'AGTA05000379.1'	136652	Standard
# #		'AGTA05000380.1'	135678	Standard
# #		'AGTA05000381.1'	135203	Standard
# #		'AGTA05000382.1'	135108	Standard
# #		'AGTA05000384.1'	134990	Standard
# #		'AGTA05000383.1'	134930	Standard
# #		'AGTA05000385.1'	134004	Standard
# #		'AGTA05000386.1'	133635	Standard
# #		'AGTA05000387.1'	133550	Standard
# #		'AGTA05000388.1'	132873	Standard
# #		'AGTA05000389.1'	132827	Standard
# #		'AGTA05000391.1'	132722	Standard
# #		'AGTA05000390.1'	132690	Standard
# #		'AGTA05000392.1'	132682	Standard
# #		'AGTA05000394.1'	131713	Standard
# #		'AGTA05000393.1'	131663	Standard
# #		'AGTA05000395.1'	131574	Standard
# #		'AGTA05000399.1'	131437	Standard
# #		'AGTA05000396.1'	131433	Standard
# #		'AGTA05000398.1'	131405	Standard
# #		'AGTA05000397.1'	131297	Standard
# #		'AGTA05000400.1'	130976	Standard
# #		'AGTA05000401.1'	130877	Standard
# #		'AGTA05000402.1'	130847	Standard
# #		'AGTA05000403.1'	130552	Standard
# #		'AGTA05000404.1'	129643	Standard
# #		'AGTA05000405.1'	128573	Standard
# #		'AGTA05000407.1'	128417	Standard
# #		'AGTA05000406.1'	128256	Standard
# #		'AGTA05000408.1'	127723	Standard
# #		'AGTA05000410.1'	127591	Standard
# #		'AGTA05000411.1'	127441	Standard
# #		'AGTA05000409.1'	127398	Standard
# #		'AGTA05000412.1'	127360	Standard
# #		'AGTA05000413.1'	127144	Standard
# #		'AGTA05000414.1'	127008	Standard
# #		'AGTA05000415.1'	126735	Standard
# #		'AGTA05000417.1'	126480	Standard
# #		'AGTA05000418.1'	126393	Standard
# #		'AGTA05000416.1'	126372	Standard
# #		'AGTA05000419.1'	125998	Standard
# #		'AGTA05000420.1'	125815	Standard
# #		'AGTA05000422.1'	125642	Standard
# #		'AGTA05000421.1'	125631	Standard
# #		'AGTA05000424.1'	125545	Standard
# #		'AGTA05000423.1'	125481	Standard
# #		'AGTA05000426.1'	125134	Standard
# #		'AGTA05000425.1'	125107	Standard
# #		'AGTA05000427.1'	124523	Standard
# #		'AGTA05000428.1'	124272	Standard
# #		'AGTA05000429.1'	123656	Standard
# #		'AGTA05000430.1'	123490	Standard
# #		'AGTA05000431.1'	123391	Standard
# #		'AGTA05000432.1'	123166	Standard
# #		'AGTA05000433.1'	122712	Standard
# #		'AGTA05000434.1'	121737	Standard
# #		'AGTA05000435.1'	121600	Standard
# #		'AGTA05000436.1'	120859	Standard
# #		'AGTA05000437.1'	120701	Standard
# #		'AGTA05000438.1'	120675	Standard
# #		'AGTA05000439.1'	119962	Standard
# #		'AGTA05000440.1'	119512	Standard
# #		'AGTA05000441.1'	119416	Standard
# #		'AGTA05000442.1'	119273	Standard
# #		'AGTA05000443.1'	119110	Standard
# #		'AGTA05000444.1'	119106	Standard
# #		'AGTA05000446.1'	118600	Standard
# #		'AGTA05000445.1'	118501	Standard
# #		'AGTA05000447.1'	118363	Standard
# #		'AGTA05000448.1'	118011	Standard
# #		'AGTA05000449.1'	117550	Standard
# #		'AGTA05000450.1'	117135	Standard
# #		'AGTA05000451.1'	116882	Standard
# #		'AGTA05000452.1'	116594	Standard
# #		'AGTA05000453.1'	116365	Standard
# #		'AGTA05000454.1'	115695	Standard
# #		'AGTA05000455.1'	115493	Standard
# #		'AGTA05000456.1'	115463	Standard
# #		'AGTA05000457.1'	114699	Standard
# #		'AGTA05000458.1'	114421	Standard
# #		'AGTA05000460.1'	114006	Standard
# #		'AGTA05000459.1'	113996	Standard
# #		'AGTA05000461.1'	113925	Standard
# #		'AGTA05000462.1'	113834	Standard
# #		'AGTA05000463.1'	113807	Standard
# #		'AGTA05000464.1'	113415	Standard
# #		'AGTA05000465.1'	113199	Standard
# #		'AGTA05000466.1'	112925	Standard
# #		'AGTA05000467.1'	112546	Standard
# #		'AGTA05000468.1'	112496	Standard
# #		'AGTA05000469.1'	112370	Standard
# #		'AGTA05000470.1'	112343	Standard
# #		'AGTA05000471.1'	111710	Standard
# #		'AGTA05000472.1'	111689	Standard
# #		'AGTA05000473.1'	111448	Standard
# #		'AGTA05000474.1'	111243	Standard
# #		'AGTA05000475.1'	110986	Standard
# #		'AGTA05000476.1'	110924	Standard
# #		'AGTA05000477.1'	110623	Standard
# #		'AGTA05000478.1'	110304	Standard
# #		'AGTA05000480.1'	110098	Standard
# #		'AGTA05000479.1'	109986	Standard
# #		'AGTA05000482.1'	109767	Standard
# #		'AGTA05000481.1'	109738	Standard
# #		'AGTA05000483.1'	109711	Standard
# #		'AGTA05000484.1'	109294	Standard
# #		'AGTA05000485.1'	108793	Standard
# #		'AGTA05000487.1'	108731	Standard
# #		'AGTA05000486.1'	108720	Standard
# #		'AGTA05000488.1'	108421	Standard
# #		'AGTA05000490.1'	108299	Standard
# #		'AGTA05000489.1'	108243	Standard
# #		'AGTA05000491.1'	108151	Standard
# #		'AGTA05000492.1'	107816	Standard
# #		'AGTA05000494.1'	107749	Standard
# #		'AGTA05000493.1'	107683	Standard
# #		'AGTA05000495.1'	107661	Standard
# #		'AGTA05000496.1'	107590	Standard
# #		'AGTA05000497.1'	107310	Standard
# #		'AGTA05000498.1'	107268	Standard
# #		'AGTA05000499.1'	107193	Standard
# #		'AGTA05000500.1'	106974	Standard
# #		'AGTA05000501.1'	106768	Standard
# #		'AGTA05000502.1'	106128	Standard
# #		'AGTA05000503.1'	105744	Standard
# #		'AGTA05000504.1'	105621	Standard
# #		'AGTA05000505.1'	105574	Standard
# #		'AGTA05000507.1'	105251	Standard
# #		'AGTA05000506.1'	105195	Standard
# #		'AGTA05000509.1'	105061	Standard
# #		'AGTA05000508.1'	105049	Standard
# #		'AGTA05000510.1'	104412	Standard
# #		'AGTA05000511.1'	104347	Standard
# #		'AGTA05000512.1'	103612	Standard
# #		'AGTA05000513.1'	103056	Standard
# #		'AGTA05000514.1'	103038	Standard
# #		'AGTA05000515.1'	102844	Standard
# #		'AGTA05000516.1'	102203	Standard
# #		'AGTA05000517.1'	101671	Standard
# #		'AGTA05000518.1'	101648	Standard
# #		'AGTA05000519.1'	101640	Standard
# #		'AGTA05000520.1'	101438	Standard
# #		'AGTA05000521.1'	101434	Standard
# #		'AGTA05000522.1'	100945	Standard
# #		'AGTA05000523.1'	100608	Standard
# #		'AGTA05000525.1'	100522	Standard
# #		'AGTA05000524.1'	100518	Standard
# #		'AGTA05000526.1'	100114	Standard
# #		'AGTA05000527.1'	99976	Standard
# #		'AGTA05000528.1'	99915	Standard
# #		'AGTA05000530.1'	99735	Standard
# #		'AGTA05000529.1'	99710	Standard
# #		'AGTA05000531.1'	99592	Standard
# #		'AGTA05000532.1'	99115	Standard
# #		'AGTA05000533.1'	98952	Standard
# #		'AGTA05000534.1'	98751	Standard
# #		'AGTA05000535.1'	98600	Standard
# #		'AGTA05000536.1'	98147	Standard
# #		'AGTA05000537.1'	97954	Standard
# #		'AGTA05000538.1'	97848	Standard
# #		'AGTA05000539.1'	97192	Standard
# #		'AGTA05000540.1'	97179	Standard
# #		'AGTA05000541.1'	97042	Standard
# #		'AGTA05000542.1'	96931	Standard
# #		'AGTA05000543.1'	96762	Standard
# #		'AGTA05000544.1'	96514	Standard
# #		'AGTA05000545.1'	95890	Standard
# #		'AGTA05000546.1'	95655	Standard
# #		'AGTA05000547.1'	95581	Standard
# #		'AGTA05000549.1'	95509	Standard
# #		'AGTA05000550.1'	95451	Standard
# #		'AGTA05000548.1'	95391	Standard
# #		'AGTA05000551.1'	95169	Standard
# #		'AGTA05000552.1'	95048	Standard
# #		'AGTA05000553.1'	95000	Standard
# #		'AGTA05000554.1'	94682	Standard
# #		'AGTA05000556.1'	94355	Standard
# #		'AGTA05000555.1'	94328	Standard
# #		'AGTA05000557.1'	94232	Standard
# #		'AGTA05000558.1'	94136	Standard
# #		'AGTA05000560.1'	93855	Standard
# #		'AGTA05000559.1'	93816	Standard
# #		'AGTA05000561.1'	93529	Standard
# #		'AGTA05000562.1'	93495	Standard
# #		'AGTA05000563.1'	93400	Standard
# #		'AGTA05000564.1'	93223	Standard
# #		'AGTA05000565.1'	92650	Standard
# #		'AGTA05000566.1'	92204	Standard
# #		'AGTA05000567.1'	91915	Standard
# #		'AGTA05000568.1'	91601	Standard
# #		'AGTA05000569.1'	91565	Standard
# #		'AGTA05000570.1'	91463	Standard
# #		'AGTA05000571.1'	91250	Standard
# #		'AGTA05000572.1'	91131	Standard
# #		'AGTA05000573.1'	91062	Standard
# #		'AGTA05000574.1'	90980	Standard
# #		'AGTA05000576.1'	90679	Standard
# #		'AGTA05000575.1'	90642	Standard
# #		'AGTA05000577.1'	90142	Standard
# #		'AGTA05000578.1'	90077	Standard
# #		'AGTA05000579.1'	89983	Standard
# #		'AGTA05000580.1'	89899	Standard
# #		'AGTA05000581.1'	89875	Standard
# #		'AGTA05000582.1'	89791	Standard
# #		'AGTA05000584.1'	89791	Standard
# #		'AGTA05000583.1'	89747	Standard
# #		'AGTA05000586.1'	89660	Standard
# #		'AGTA05000585.1'	89649	Standard
# #		'AGTA05000587.1'	89389	Standard
# #		'AGTA05000589.1'	89233	Standard
# #		'AGTA05000588.1'	89223	Standard
# #		'AGTA05000590.1'	89078	Standard
# #		'AGTA05000591.1'	89005	Standard
# #		'AGTA05000592.1'	88880	Standard
# #		'AGTA05000594.1'	88757	Standard
# #		'AGTA05000593.1'	88753	Standard
# #		'AGTA05000595.1'	88677	Standard
# #		'AGTA05000596.1'	88532	Standard
# #		'AGTA05000597.1'	88453	Standard
# #		'AGTA05000598.1'	88377	Standard
# #		'AGTA05000599.1'	88282	Standard
# #		'AGTA05000600.1'	88056	Standard
# #		'AGTA05000601.1'	87645	Standard
# #		'AGTA05000603.1'	87223	Standard
# #		'AGTA05000602.1'	87222	Standard
# #		'AGTA05000604.1'	87069	Standard
# #		'AGTA05000605.1'	86577	Standard
# #		'AGTA05000606.1'	86312	Standard
# #		'AGTA05000607.1'	85869	Standard
# #		'AGTA05000608.1'	85868	Standard
# #		'AGTA05000611.1'	85733	Standard
# #		'AGTA05000612.1'	85719	Standard
# #		'AGTA05000610.1'	85711	Standard
# #		'AGTA05000609.1'	85710	Standard
# #		'AGTA05000613.1'	85324	Standard
# #		'AGTA05000614.1'	85221	Standard
# #		'AGTA05000615.1'	84988	Standard
# #		'AGTA05000616.1'	84670	Standard
# #		'AGTA05000617.1'	84568	Standard
# #		'AGTA05000618.1'	84534	Standard
# #		'AGTA05000619.1'	84499	Standard
# #		'AGTA05000620.1'	84316	Standard
# #		'AGTA05000621.1'	84205	Standard
# #		'AGTA05000622.1'	84072	Standard
# #		'AGTA05000623.1'	83893	Standard
# #		'AGTA05000624.1'	83806	Standard
# #		'AGTA05000625.1'	83775	Standard
# #		'AGTA05000626.1'	83489	Standard
# #		'AGTA05000627.1'	83156	Standard
# #		'AGTA05000628.1'	83137	Standard
# #		'AGTA05000629.1'	83130	Standard
# #		'AGTA05000631.1'	82901	Standard
# #		'AGTA05000630.1'	82849	Standard
# #		'AGTA05000632.1'	82848	Standard
# #		'AGTA05000633.1'	82377	Standard
# #		'AGTA05000634.1'	82096	Standard
# #		'AGTA05000635.1'	82074	Standard
# #		'AGTA05000636.1'	81828	Standard
# #		'AGTA05000637.1'	81688	Standard
# #		'AGTA05000638.1'	81649	Standard
# #		'AGTA05000639.1'	81538	Standard
# #		'AGTA05000640.1'	81396	Standard
# #		'AGTA05000642.1'	81351	Standard
# #		'AGTA05000641.1'	81322	Standard
# #		'AGTA05000643.1'	81282	Standard
# #		'AGTA05000644.1'	81270	Standard
# #		'AGTA05000645.1'	81202	Standard
# #		'AGTA05000646.1'	80919	Standard
# #		'AGTA05000647.1'	80786	Standard
# #		'AGTA05000650.1'	80784	Standard
# #		'AGTA05000649.1'	80778	Standard
# #		'AGTA05000648.1'	80764	Standard
# #		'AGTA05000651.1'	80691	Standard
# #		'AGTA05000652.1'	80592	Standard
# #		'AGTA05000653.1'	80574	Standard
# #		'AGTA05000654.1'	80354	Standard
# #		'AGTA05000655.1'	80230	Standard
# #		'AGTA05000657.1'	80077	Standard
# #		'AGTA05000656.1'	80020	Standard
# #		'AGTA05000658.1'	79541	Standard
# #		'AGTA05000660.1'	79336	Standard
# #		'AGTA05000659.1'	79320	Standard
# #		'AGTA05000662.1'	79206	Standard
# #		'AGTA05000661.1'	79141	Standard
# #		'AGTA05000663.1'	78825	Standard
# #		'AGTA05000664.1'	78756	Standard
# #		'AGTA05000667.1'	78391	Standard
# #		'AGTA05000666.1'	78381	Standard
# #		'AGTA05000665.1'	78366	Standard
# #		'AGTA05000668.1'	78322	Standard
# #		'AGTA05000669.1'	78240	Standard
# #		'AGTA05000670.1'	78116	Standard
# #		'AGTA05000673.1'	78046	Standard
# #		'AGTA05000672.1'	78042	Standard
# #		'AGTA05000671.1'	78024	Standard
# #		'AGTA05000674.1'	77829	Standard
# #		'AGTA05000675.1'	77791	Standard
# #		'AGTA05000678.1'	77730	Standard
# #		'AGTA05000677.1'	77729	Standard
# #		'AGTA05000676.1'	77714	Standard
# #		'AGTA05000679.1'	77590	Standard
# #		'AGTA05000680.1'	77247	Standard
# #		'AGTA05000682.1'	77120	Standard
# #		'AGTA05000681.1'	77084	Standard
# #		'AGTA05000683.1'	77048	Standard
# #		'AGTA05000684.1'	76837	Standard
# #		'AGTA05000685.1'	76826	Standard
# #		'AGTA05000686.1'	76732	Standard
# #		'AGTA05000687.1'	76644	Standard
# #		'AGTA05000688.1'	76146	Standard
# #		'AGTA05000689.1'	76040	Standard
# #		'AGTA05000690.1'	75989	Standard
# #		'AGTA05000691.1'	75847	Standard
# #		'AGTA05000692.1'	75746	Standard
# #		'AGTA05000693.1'	75653	Standard
# #		'AGTA05000694.1'	75631	Standard
# #		'AGTA05000695.1'	75553	Standard
# #		'AGTA05000696.1'	75533	Standard
# #		'AGTA05000697.1'	75443	Standard
# #		'AGTA05000698.1'	75391	Standard
# #		'AGTA05000699.1'	75312	Standard
# #		'AGTA05000700.1'	75188	Standard
# #		'AGTA05000701.1'	75144	Standard
# #		'AGTA05000702.1'	75087	Standard
# #		'AGTA05000703.1'	74993	Standard
# #		'AGTA05000704.1'	74789	Standard
# #		'AGTA05000705.1'	74754	Standard
# #		'AGTA05000707.1'	74687	Standard
# #		'AGTA05000706.1'	74669	Standard
# #		'AGTA05000710.1'	74642	Standard
# #		'AGTA05000708.1'	74611	Standard
# #		'AGTA05000709.1'	74603	Standard
# #		'AGTA05000711.1'	74517	Standard
# #		'AGTA05000712.1'	74376	Standard
# #		'AGTA05000713.1'	74219	Standard
# #		'AGTA05000714.1'	74140	Standard
# #		'AGTA05000716.1'	73956	Standard
# #		'AGTA05000715.1'	73930	Standard
# #		'AGTA05000717.1'	73822	Standard
# #		'AGTA05000718.1'	73672	Standard
# #		'AGTA05000719.1'	73608	Standard
# #		'AGTA05000721.1'	73438	Standard
# #		'AGTA05000720.1'	73385	Standard
# #		'AGTA05000722.1'	73307	Standard
# #		'AGTA05000723.1'	73238	Standard
# #		'AGTA05000724.1'	73118	Standard
# #		'AGTA05000725.1'	73031	Standard
# #		'AGTA05000726.1'	72948	Standard
# #		'AGTA05000727.1'	72856	Standard
# #		'AGTA05000728.1'	72718	Standard
# #		'AGTA05000729.1'	72526	Standard
# #		'AGTA05000730.1'	72504	Standard
# #		'AGTA05000731.1'	72413	Standard
# #		'AGTA05000732.1'	72306	Standard
# #		'AGTA05000733.1'	72181	Standard
# #		'AGTA05000734.1'	72166	Standard
# #		'AGTA05000735.1'	72079	Standard
# #		'AGTA05000736.1'	71979	Standard
# #		'AGTA05000737.1'	71959	Standard
# #		'AGTA05000739.1'	71850	Standard
# #		'AGTA05000738.1'	71829	Standard
# #		'AGTA05000740.1'	71686	Standard
# #		'AGTA05000741.1'	71667	Standard
# #		'AGTA05000742.1'	71443	Standard
# #		'AGTA05000743.1'	71291	Standard
# #		'AGTA05000745.1'	71255	Standard
# #		'AGTA05000744.1'	71186	Standard
# #		'AGTA05000746.1'	71137	Standard
# #		'AGTA05000747.1'	71116	Standard
# #		'AGTA05000748.1'	71031	Standard
# #		'AGTA05000749.1'	70831	Standard
# #		'AGTA05000750.1'	70581	Standard
# #		'AGTA05000751.1'	70555	Standard
# #		'AGTA05000752.1'	70470	Standard
# #		'AGTA05000753.1'	70380	Standard
# #		'AGTA05000754.1'	70306	Standard
# #		'AGTA05000756.1'	69831	Standard
# #		'AGTA05000755.1'	69758	Standard
# #		'AGTA05000757.1'	69698	Standard
# #		'AGTA05000758.1'	69672	Standard
# #		'AGTA05000759.1'	69556	Standard
# #		'AGTA05000760.1'	69534	Standard
# #		'AGTA05000761.1'	69494	Standard
# #		'AGTA05000762.1'	69382	Standard
# #		'AGTA05000763.1'	69349	Standard
# #		'AGTA05000764.1'	69170	Standard
# #		'AGTA05000765.1'	69072	Standard
# #		'AGTA05000766.1'	68972	Standard
# #		'AGTA05000767.1'	68798	Standard
# #		'AGTA05000768.1'	68688	Standard
# #		'AGTA05000769.1'	68373	Standard
# #		'AGTA05000770.1'	68360	Standard
# #		'AGTA05000771.1'	68020	Standard
# #		'AGTA05000772.1'	67793	Standard
# #		'AGTA05000773.1'	67760	Standard
# #		'AGTA05000774.1'	67743	Standard
# #		'AGTA05000777.1'	67741	Standard
# #		'AGTA05000775.1'	67737	Standard
# #		'AGTA05000778.1'	67718	Standard
# #		'AGTA05000776.1'	67700	Standard
# #		'AGTA05000779.1'	67694	Standard
# #		'AGTA05000780.1'	67669	Standard
# #		'AGTA05000781.1'	67433	Standard
# #		'AGTA05000783.1'	67408	Standard
# #		'AGTA05000782.1'	67377	Standard
# #		'AGTA05000784.1'	67370	Standard
# #		'AGTA05000785.1'	67178	Standard
# #		'AGTA05000786.1'	66946	Standard
# #		'AGTA05000787.1'	66840	Standard
# #		'AGTA05000789.1'	66745	Standard
# #		'AGTA05000788.1'	66726	Standard
# #		'AGTA05000791.1'	66667	Standard
# #		'AGTA05000792.1'	66633	Standard
# #		'AGTA05000790.1'	66606	Standard
# #		'AGTA05000793.1'	66596	Standard
# #		'AGTA05000794.1'	66574	Standard
# #		'AGTA05000795.1'	66557	Standard
# #		'AGTA05000796.1'	66340	Standard
# #		'AGTA05000798.1'	66014	Standard
# #		'AGTA05000797.1'	66009	Standard
# #		'AGTA05000800.1'	65949	Standard
# #		'AGTA05000799.1'	65935	Standard
# #		'AGTA05000802.1'	65913	Standard
# #		'AGTA05000801.1'	65874	Standard
# #		'AGTA05000803.1'	65631	Standard
# #		'AGTA05000804.1'	65593	Standard
# #		'AGTA05000805.1'	65519	Standard
# #		'AGTA05000806.1'	65501	Standard
# #		'AGTA05000807.1'	65269	Standard
# #		'AGTA05000808.1'	65198	Standard
# #		'AGTA05000809.1'	65178	Standard
# #		'AGTA05000810.1'	65057	Standard
# #		'AGTA05000811.1'	65032	Standard
# #		'AGTA05000812.1'	64998	Standard
# #		'AGTA05000814.1'	64835	Standard
# #		'AGTA05000813.1'	64814	Standard
# #		'AGTA05000815.1'	64779	Standard
# #		'AGTA05000816.1'	64753	Standard
# #		'AGTA05000817.1'	64716	Standard
# #		'AGTA05000819.1'	64648	Standard
# #		'AGTA05000818.1'	64642	Standard
# #		'AGTA05000820.1'	64429	Standard
# #		'AGTA05000821.1'	64249	Standard
# #		'AGTA05000822.1'	64140	Standard
# #		'AGTA05000823.1'	64064	Standard
# #		'AGTA05000825.1'	64019	Standard
# #		'AGTA05000826.1'	64009	Standard
# #		'AGTA05000824.1'	63991	Standard
# #		'AGTA05000828.1'	63898	Standard
# #		'AGTA05000827.1'	63894	Standard
# #		'AGTA05000829.1'	63870	Standard
# #		'AGTA05000830.1'	63837	Standard
# #		'AGTA05000831.1'	63834	Standard
# #		'AGTA05000832.1'	63795	Standard
# #		'AGTA05000833.1'	63739	Standard
# #		'AGTA05000834.1'	63582	Standard
# #		'AGTA05000835.1'	63331	Standard
# #		'AGTA05000837.1'	63228	Standard
# #		'AGTA05000836.1'	63215	Standard
# #		'AGTA05000838.1'	63151	Standard
# #		'AGTA05000839.1'	63033	Standard
# #		'AGTA05000840.1'	62877	Standard
# #		'AGTA05000841.1'	62809	Standard
# #		'AGTA05000842.1'	62690	Standard
# #		'AGTA05000843.1'	62572	Standard
# #		'AGTA05000844.1'	62496	Standard
# #		'AGTA05000846.1'	62399	Standard
# #		'AGTA05000845.1'	62387	Standard
# #		'AGTA05000847.1'	62134	Standard
# #		'AGTA05000848.1'	62073	Standard
# #		'AGTA05000849.1'	61928	Standard
# #		'AGTA05000850.1'	61850	Standard
# #		'AGTA05000851.1'	61781	Standard
# #		'AGTA05000852.1'	61762	Standard
# #		'AGTA05000853.1'	61730	Standard
# #		'AGTA05000854.1'	61605	Standard
# #		'AGTA05000857.1'	61557	Standard
# #		'AGTA05000855.1'	61551	Standard
# #		'AGTA05000858.1'	61533	Standard
# #		'AGTA05000856.1'	61490	Standard
# #		'AGTA05000859.1'	61441	Standard
# #		'AGTA05000860.1'	61427	Standard
# #		'AGTA05000861.1'	61400	Standard
# #		'AGTA05000862.1'	61352	Standard
# #		'AGTA05000863.1'	61246	Standard
# #		'AGTA05000864.1'	61222	Standard
# #		'AGTA05000865.1'	61179	Standard
# #		'AGTA05000866.1'	61179	Standard
# #		'AGTA05000867.1'	61147	Standard
# #		'AGTA05000868.1'	61062	Standard
# #		'AGTA05000869.1'	61029	Standard
# #		'AGTA05000870.1'	60922	Standard
# #		'AGTA05000871.1'	60618	Standard
# #		'AGTA05000872.1'	60596	Standard
# #		'AGTA05000873.1'	60495	Standard
# #		'AGTA05000874.1'	60412	Standard
# #		'AGTA05000875.1'	60412	Standard
# #		'AGTA05000876.1'	60268	Standard
# #		'AGTA05000878.1'	60160	Standard
# #		'AGTA05000877.1'	60152	Standard
# #		'AGTA05000879.1'	59978	Standard
# #		'AGTA05000880.1'	59748	Standard
# #		'AGTA05000881.1'	59685	Standard
# #		'AGTA05000882.1'	59581	Standard
# #		'AGTA05000883.1'	59413	Standard
# #		'AGTA05000884.1'	59338	Standard
# #		'AGTA05000885.1'	59087	Standard
# #		'AGTA05000886.1'	59019	Standard
# #		'AGTA05000887.1'	58844	Standard
# #		'AGTA05000888.1'	58648	Standard
# #		'AGTA05000889.1'	58503	Standard
# #		'AGTA05000890.1'	58501	Standard
# #		'AGTA05000891.1'	58480	Standard
# #		'AGTA05000892.1'	58339	Standard
# #		'AGTA05000893.1'	58276	Standard
# #		'AGTA05000894.1'	58263	Standard
# #		'AGTA05000895.1'	58219	Standard
# #		'AGTA05000896.1'	58188	Standard
# #		'AGTA05000897.1'	58119	Standard
# #		'AGTA05000898.1'	58060	Standard
# #		'AGTA05000900.1'	57984	Standard
# #		'AGTA05000899.1'	57981	Standard
# #		'AGTA05000901.1'	57875	Standard
# #		'AGTA05000903.1'	57759	Standard
# #		'AGTA05000902.1'	57737	Standard
# #		'AGTA05000904.1'	57737	Standard
# #		'AGTA05000905.1'	57722	Standard
# #		'AGTA05000906.1'	57639	Standard
# #		'AGTA05000907.1'	57634	Standard
# #		'AGTA05000908.1'	57553	Standard
# #		'AGTA05000909.1'	57517	Standard
# #		'AGTA05000910.1'	57429	Standard
# #		'AGTA05000911.1'	57365	Standard
# #		'AGTA05000912.1'	57292	Standard
# #		'AGTA05000913.1'	57230	Standard
# #		'AGTA05000914.1'	57166	Standard
# #		'AGTA05000915.1'	57155	Standard
# #		'AGTA05000916.1'	57079	Standard
# #		'AGTA05000917.1'	57043	Standard
# #		'AGTA05000918.1'	56916	Standard
# #		'AGTA05000919.1'	56831	Standard
# #		'AGTA05000920.1'	56780	Standard
# #		'AGTA05000921.1'	56747	Standard
# #		'AGTA05000922.1'	56607	Standard
# #		'AGTA05000923.1'	56500	Standard
# #		'AGTA05000925.1'	56405	Standard
# #		'AGTA05000924.1'	56380	Standard
# #		'AGTA05000926.1'	56324	Standard
# #		'AGTA05000927.1'	56260	Standard
# #		'AGTA05000928.1'	56236	Standard
# #		'AGTA05000929.1'	56166	Standard
# #		'AGTA05000930.1'	56135	Standard
# #		'AGTA05000931.1'	56078	Standard
# #		'AGTA05000932.1'	55930	Standard
# #		'AGTA05000933.1'	55737	Standard
# #		'AGTA05000934.1'	55632	Standard
# #		'AGTA05000935.1'	55529	Standard
# #		'AGTA05000936.1'	55510	Standard
# #		'AGTA05000938.1'	55419	Standard
# #		'AGTA05000937.1'	55418	Standard
# #		'AGTA05000939.1'	55261	Standard
# #		'AGTA05000940.1'	55197	Standard
# #		'AGTA05000943.1'	55003	Standard
# #		'AGTA05000941.1'	54997	Standard
# #		'AGTA05000942.1'	54972	Standard
# #		'AGTA05000944.1'	54968	Standard
# #		'AGTA05000946.1'	54933	Standard
# #		'AGTA05000945.1'	54904	Standard
# #		'AGTA05000947.1'	54848	Standard
# #		'AGTA05000948.1'	54798	Standard
# #		'AGTA05000949.1'	54611	Standard
# #		'AGTA05000950.1'	54544	Standard
# #		'AGTA05000951.1'	54515	Standard
# #		'AGTA05000952.1'	54451	Standard
# #		'AGTA05000953.1'	54415	Standard
# #		'AGTA05000954.1'	54349	Standard
# #		'AGTA05000955.1'	54330	Standard
# #		'AGTA05000956.1'	54237	Standard
# #		'AGTA05000957.1'	54072	Standard
# #		'AGTA05000959.1'	54070	Standard
# #		'AGTA05000960.1'	53993	Standard
# #		'AGTA05000958.1'	53983	Standard
# #		'AGTA05000961.1'	53822	Standard
# #		'AGTA05000962.1'	53802	Standard
# #		'AGTA05000964.1'	53788	Standard
# #		'AGTA05000963.1'	53785	Standard
# #		'AGTA05000966.1'	53722	Standard
# #		'AGTA05000965.1'	53713	Standard
# #		'AGTA05000967.1'	53571	Standard
# #		'AGTA05000969.1'	53477	Standard
# #		'AGTA05000968.1'	53474	Standard
# #		'AGTA05000970.1'	53366	Standard
# #		'AGTA05000971.1'	53056	Standard
# #		'AGTA05000972.1'	53017	Standard
# #		'AGTA05000973.1'	52997	Standard
# #		'AGTA05000974.1'	52993	Standard
# #		'AGTA05000975.1'	52927	Standard
# #		'AGTA05000976.1'	52890	Standard
# #		'AGTA05000977.1'	52888	Standard
# #		'AGTA05000978.1'	52847	Standard
# #		'AGTA05000980.1'	52698	Standard
# #		'AGTA05000979.1'	52682	Standard
# #		'AGTA05000981.1'	52584	Standard
# #		'AGTA05000982.1'	52542	Standard
# #		'AGTA05000983.1'	52534	Standard
# #		'AGTA05000985.1'	52492	Standard
# #		'AGTA05000984.1'	52482	Standard
# #		'AGTA05000986.1'	52433	Standard
# #		'AGTA05000987.1'	52281	Standard
# #		'AGTA05000988.1'	52225	Standard
# #		'AGTA05000990.1'	52221	Standard
# #		'AGTA05000989.1'	52169	Standard
# #		'AGTA05000991.1'	52103	Standard
# #		'AGTA05000993.1'	52037	Standard
# #		'AGTA05000992.1'	52012	Standard
# #		'AGTA05000994.1'	51669	Standard
# #		'AGTA05000995.1'	51666	Standard
# #		'AGTA05000996.1'	51537	Standard
# #		'AGTA05000997.1'	51345	Standard
# #		'AGTA05000998.1'	51293	Standard
# #		'AGTA05001001.1'	51275	Standard
# #		'AGTA05001002.1'	51220	Standard
# #		'AGTA05001000.1'	51213	Standard
# #		'AGTA05000999.1'	51200	Standard
# #		'AGTA05001003.1'	51186	Standard
# #		'AGTA05001004.1'	51111	Standard
# #		'AGTA05001005.1'	50992	Standard
# #		'AGTA05001007.1'	50982	Standard
# #		'AGTA05001006.1'	50950	Standard
# #		'AGTA05001008.1'	50872	Standard
# #		'AGTA05001009.1'	50846	Standard
# #		'AGTA05001011.1'	50789	Standard
# #		'AGTA05001010.1'	50777	Standard
# #		'AGTA05001012.1'	50704	Standard
# #		'AGTA05001013.1'	50691	Standard
# #		'AGTA05001014.1'	50681	Standard
# #		'AGTA05001015.1'	50627	Standard
# #		'AGTA05001016.1'	50599	Standard
# #		'AGTA05001017.1'	50550	Standard
# #		'AGTA05001018.1'	50468	Standard
# #		'AGTA05001019.1'	50456	Standard
# #		'AGTA05001020.1'	50324	Standard
# #		'AGTA05001022.1'	50283	Standard
# #		'AGTA05001021.1'	50272	Standard
# #		'AGTA05001023.1'	50217	Standard
# #		'AGTA05001024.1'	50122	Standard
# #		'AGTA05001025.1'	50071	Standard
# #		'AGTA05001026.1'	50054	Standard
# #		'AGTA05001027.1'	49967	Standard
# #		'AGTA05001028.1'	49930	Standard
# #		'AGTA05001030.1'	49846	Standard
# #		'AGTA05001029.1'	49829	Standard
# #		'AGTA05001031.1'	49696	Standard
# #		'AGTA05001032.1'	49634	Standard
# #		'AGTA05001033.1'	49407	Standard
# #		'AGTA05001034.1'	49311	Standard
# #		'AGTA05001035.1'	49289	Standard
# #		'AGTA05001036.1'	49237	Standard
# #		'AGTA05001037.1'	49231	Standard
# #		'AGTA05001039.1'	49209	Standard
# #		'AGTA05001038.1'	49181	Standard
# #		'AGTA05001040.1'	49167	Standard
# #		'AGTA05001041.1'	49112	Standard
# #		'AGTA05001042.1'	49088	Standard
# #		'AGTA05001043.1'	48981	Standard
# #		'AGTA05001044.1'	48970	Standard
# #		'AGTA05001045.1'	48846	Standard
# #		'AGTA05001046.1'	48555	Standard
# #		'AGTA05001047.1'	48548	Standard
# #		'AGTA05001048.1'	48493	Standard
# #		'AGTA05001050.1'	48369	Standard
# #		'AGTA05001049.1'	48367	Standard
# #		'AGTA05001051.1'	48352	Standard
# #		'AGTA05001052.1'	48299	Standard
# #		'AGTA05001053.1'	48218	Standard
# #		'AGTA05001054.1'	48199	Standard
# #		'AGTA05001055.1'	48030	Standard
# #		'AGTA05001056.1'	47891	Standard
# #		'AGTA05001058.1'	47826	Standard
# #		'AGTA05001057.1'	47823	Standard
# #		'AGTA05001059.1'	47694	Standard
# #		'AGTA05001060.1'	47684	Standard
# #		'AGTA05001061.1'	47606	Standard
# #		'AGTA05001062.1'	47603	Standard
# #		'AGTA05001063.1'	47584	Standard
# #		'AGTA05001064.1'	47583	Standard
# #		'AGTA05001065.1'	47495	Standard
# #		'AGTA05001066.1'	47474	Standard
# #		'AGTA05001067.1'	47452	Standard
# #		'AGTA05001068.1'	47389	Standard
# #		'AGTA05001070.1'	47219	Standard
# #		'AGTA05001069.1'	47156	Standard
# #		'AGTA05001072.1'	47133	Standard
# #		'AGTA05001071.1'	47124	Standard
# #		'AGTA05001073.1'	46994	Standard
# #		'AGTA05001074.1'	46994	Standard
# #		'AGTA05001075.1'	46944	Standard
# #		'AGTA05001077.1'	46910	Standard
# #		'AGTA05001076.1'	46898	Standard
# #		'AGTA05001078.1'	46884	Standard
# #		'AGTA05001079.1'	46866	Standard
# #		'AGTA05001080.1'	46846	Standard
# #		'AGTA05001082.1'	46806	Standard
# #		'AGTA05001081.1'	46791	Standard
# #		'AGTA05001084.1'	46764	Standard
# #		'AGTA05001083.1'	46738	Standard
# #		'AGTA05001086.1'	46671	Standard
# #		'AGTA05001085.1'	46652	Standard
# #		'AGTA05001087.1'	46320	Standard
# #		'AGTA05001089.1'	46318	Standard
# #		'AGTA05001088.1'	46274	Standard
# #		'AGTA05001090.1'	46199	Standard
# #		'AGTA05001091.1'	46148	Standard
# #		'AGTA05001092.1'	46013	Standard
# #		'AGTA05001093.1'	45989	Standard
# #		'AGTA05001094.1'	45937	Standard
# #		'AGTA05001095.1'	45901	Standard
# #		'AGTA05001096.1'	45899	Standard
# #		'AGTA05001098.1'	45840	Standard
# #		'AGTA05001097.1'	45837	Standard
# #		'AGTA05001099.1'	45744	Standard
# #		'AGTA05001100.1'	45704	Standard
# #		'AGTA05001102.1'	45700	Standard
# #		'AGTA05001101.1'	45694	Standard
# #		'AGTA05001104.1'	45692	Standard
# #		'AGTA05001103.1'	45658	Standard
# #		'AGTA05001105.1'	45494	Standard
# #		'AGTA05001106.1'	45463	Standard
# #		'AGTA05001107.1'	45385	Standard
# #		'AGTA05001108.1'	45288	Standard
# #		'AGTA05001109.1'	45285	Standard
# #		'AGTA05001110.1'	45233	Standard
# #		'AGTA05001112.1'	45192	Standard
# #		'AGTA05001111.1'	45172	Standard
# #		'AGTA05001113.1'	45159	Standard
# #		'AGTA05001114.1'	45147	Standard
# #		'AGTA05001115.1'	45092	Standard
# #		'AGTA05001116.1'	45065	Standard
# #		'AGTA05001118.1'	45029	Standard
# #		'AGTA05001117.1'	44987	Standard
# #		'AGTA05001119.1'	44961	Standard
# #		'AGTA05001120.1'	44944	Standard
# #		'AGTA05001121.1'	44861	Standard
# #		'AGTA05001123.1'	44826	Standard
# #		'AGTA05001122.1'	44818	Standard
# #		'AGTA05001124.1'	44685	Standard
# #		'AGTA05001125.1'	44649	Standard
# #		'AGTA05001126.1'	44517	Standard
# #		'AGTA05001127.1'	44487	Standard
# #		'AGTA05001128.1'	44465	Standard
# #		'AGTA05001129.1'	44284	Standard
# #		'AGTA05001130.1'	44101	Standard
# #		'AGTA05001131.1'	44066	Standard
# #		'AGTA05001133.1'	43912	Standard
# #		'AGTA05001132.1'	43911	Standard
# #		'AGTA05001134.1'	43889	Standard
# #		'AGTA05001136.1'	43858	Standard
# #		'AGTA05001135.1'	43848	Standard
# #		'AGTA05001137.1'	43811	Standard
# #		'AGTA05001139.1'	43635	Standard
# #		'AGTA05001138.1'	43615	Standard
# #		'AGTA05001141.1'	43523	Standard
# #		'AGTA05001140.1'	43467	Standard
# #		'AGTA05001142.1'	43461	Standard
# #		'AGTA05001143.1'	43396	Standard
# #		'AGTA05001144.1'	43382	Standard
# #		'AGTA05001145.1'	43268	Standard
# #		'AGTA05001146.1'	43221	Standard
# #		'AGTA05001147.1'	43094	Standard
# #		'AGTA05001148.1'	43038	Standard
# #		'AGTA05001149.1'	43006	Standard
# #		'AGTA05001150.1'	42923	Standard
# #		'AGTA05001151.1'	42876	Standard
# #		'AGTA05001153.1'	42853	Standard
# #		'AGTA05001152.1'	42843	Standard
# #		'AGTA05001155.1'	42792	Standard
# #		'AGTA05001157.1'	42785	Standard
# #		'AGTA05001154.1'	42779	Standard
# #		'AGTA05001156.1'	42779	Standard
# #		'AGTA05001159.1'	42741	Standard
# #		'AGTA05001158.1'	42736	Standard
# #		'AGTA05001160.1'	42654	Standard
# #		'AGTA05001161.1'	42633	Standard
# #		'AGTA05001162.1'	42586	Standard
# #		'AGTA05001163.1'	42527	Standard
# #		'AGTA05001164.1'	42495	Standard
# #		'AGTA05001166.1'	42254	Standard
# #		'AGTA05001165.1'	42252	Standard
# #		'AGTA05001167.1'	42054	Standard
# #		'AGTA05001168.1'	41919	Standard
# #		'AGTA05001169.1'	41892	Standard
# #		'AGTA05001170.1'	41862	Standard
# #		'AGTA05001171.1'	41855	Standard
# #		'AGTA05001172.1'	41796	Standard
# #		'AGTA05001173.1'	41796	Standard
# #		'AGTA05001174.1'	41762	Standard
# #		'AGTA05001175.1'	41760	Standard
# #		'AGTA05001176.1'	41641	Standard
# #		'AGTA05001179.1'	41568	Standard
# #		'AGTA05001178.1'	41559	Standard
# #		'AGTA05001177.1'	41551	Standard
# #		'AGTA05001180.1'	41508	Standard
# #		'AGTA05001181.1'	41476	Standard
# #		'AGTA05001182.1'	41471	Standard
# #		'AGTA05001183.1'	41431	Standard
# #		'AGTA05001184.1'	41403	Standard
# #		'AGTA05001186.1'	41332	Standard
# #		'AGTA05001185.1'	41330	Standard
# #		'AGTA05001188.1'	41306	Standard
# #		'AGTA05001187.1'	41304	Standard
# #		'AGTA05001189.1'	41247	Standard
# #		'AGTA05001190.1'	41199	Standard
# #		'AGTA05001192.1'	41053	Standard
# #		'AGTA05001191.1'	41050	Standard
# #		'AGTA05001193.1'	40991	Standard
# #		'AGTA05001194.1'	40898	Standard
# #		'AGTA05001198.1'	40848	Standard
# #		'AGTA05001195.1'	40824	Standard
# #		'AGTA05001196.1'	40820	Standard
# #		'AGTA05001197.1'	40779	Standard
# #		'AGTA05001200.1'	40727	Standard
# #		'AGTA05001201.1'	40724	Standard
# #		'AGTA05001199.1'	40720	Standard
# #		'AGTA05001203.1'	40549	Standard
# #		'AGTA05001202.1'	40547	Standard
# #		'AGTA05001204.1'	40524	Standard
# #		'AGTA05001205.1'	40507	Standard
# #		'AGTA05001206.1'	40458	Standard
# #		'AGTA05001207.1'	40411	Standard
# #		'AGTA05001208.1'	40389	Standard
# #		'AGTA05001209.1'	40371	Standard
# #		'AGTA05001210.1'	40336	Standard
# #		'AGTA05001212.1'	40325	Standard
# #		'AGTA05001211.1'	40264	Standard
# #		'AGTA05001213.1'	40103	Standard
# #		'AGTA05001214.1'	40042	Standard
# #		'AGTA05001215.1'	39956	Standard
# #		'AGTA05001216.1'	39931	Standard
# #		'AGTA05001218.1'	39825	Standard
# #		'AGTA05001217.1'	39813	Standard
# #		'AGTA05001219.1'	39630	Standard
# #		'AGTA05001220.1'	39580	Standard
# #		'AGTA05001221.1'	39557	Standard
# #		'AGTA05001222.1'	39471	Standard
# #		'AGTA05001223.1'	39467	Standard
# #		'AGTA05001224.1'	39452	Standard
# #		'AGTA05001225.1'	39419	Standard
# #		'AGTA05001226.1'	39225	Standard
# #		'AGTA05001228.1'	39160	Standard
# #		'AGTA05001227.1'	39155	Standard
# #		'AGTA05001229.1'	39094	Standard
# #		'AGTA05001230.1'	39046	Standard
# #		'AGTA05001231.1'	38991	Standard
# #		'AGTA05001232.1'	38961	Standard
# #		'AGTA05001234.1'	38925	Standard
# #		'AGTA05001233.1'	38913	Standard
# #		'AGTA05001235.1'	38897	Standard
# #		'AGTA05001237.1'	38848	Standard
# #		'AGTA05001238.1'	38823	Standard
# #		'AGTA05001236.1'	38817	Standard
# #		'AGTA05001239.1'	38809	Standard
# #		'AGTA05001240.1'	38753	Standard
# #		'AGTA05001241.1'	38654	Standard
# #		'AGTA05001242.1'	38633	Standard
# #		'AGTA05001243.1'	38581	Standard
# #		'AGTA05001244.1'	38555	Standard
# #		'AGTA05001245.1'	38490	Standard
# #		'AGTA05001246.1'	38463	Standard
# #		'AGTA05001247.1'	38380	Standard
# #		'AGTA05001248.1'	38329	Standard
# #		'AGTA05001249.1'	38221	Standard
# #		'AGTA05001250.1'	38185	Standard
# #		'AGTA05001251.1'	38142	Standard
# #		'AGTA05001252.1'	38137	Standard
# #		'AGTA05001253.1'	37997	Standard
# #		'AGTA05001254.1'	37978	Standard
# #		'AGTA05001255.1'	37600	Standard
# #		'AGTA05001256.1'	37545	Standard
# #		'AGTA05001257.1'	37533	Standard
# #		'AGTA05001258.1'	37400	Standard
# #		'AGTA05001260.1'	37345	Standard
# #		'AGTA05001259.1'	37329	Standard
# #		'AGTA05001261.1'	37270	Standard
# #		'AGTA05001262.1'	37220	Standard
# #		'AGTA05001263.1'	37120	Standard
# #		'AGTA05001264.1'	37086	Standard
# #		'AGTA05001265.1'	37060	Standard
# #		'AGTA05001266.1'	37022	Standard
# #		'AGTA05001267.1'	36977	Standard
# #		'AGTA05001268.1'	36922	Standard
# #		'AGTA05001269.1'	36883	Standard
# #		'AGTA05001270.1'	36833	Standard
# #		'AGTA05001271.1'	36776	Standard
# #		'AGTA05001272.1'	36711	Standard
# #		'AGTA05001273.1'	36695	Standard
# #		'AGTA05001274.1'	36491	Standard
# #		'AGTA05001275.1'	36454	Standard
# #		'AGTA05001276.1'	36368	Standard
# #		'AGTA05001277.1'	36298	Standard
# #		'AGTA05001280.1'	36239	Standard
# #		'AGTA05001278.1'	36220	Standard
# #		'AGTA05001279.1'	36213	Standard
# #		'AGTA05001281.1'	36155	Standard
# #		'AGTA05001282.1'	36155	Standard
# #		'AGTA05001285.1'	36123	Standard
# #		'AGTA05001283.1'	36113	Standard
# #		'AGTA05001284.1'	36098	Standard
# #		'AGTA05001286.1'	36087	Standard
# #		'AGTA05001288.1'	36027	Standard
# #		'AGTA05001287.1'	36016	Standard
# #		'AGTA05001289.1'	35877	Standard
# #		'AGTA05001290.1'	35820	Standard
# #		'AGTA05001291.1'	35800	Standard
# #		'AGTA05001292.1'	35647	Standard
# #		'AGTA05001293.1'	35614	Standard
# #		'AGTA05001294.1'	35598	Standard
# #		'AGTA05001295.1'	35588	Standard
# #		'AGTA05001296.1'	35547	Standard
# #		'AGTA05001299.1'	35499	Standard
# #		'AGTA05001297.1'	35460	Standard
# #		'AGTA05001300.1'	35460	Standard
# #		'AGTA05001301.1'	35441	Standard
# #		'AGTA05001298.1'	35438	Standard
# #		'AGTA05001303.1'	35434	Standard
# #		'AGTA05001302.1'	35430	Standard
# #		'AGTA05001304.1'	35312	Standard
# #		'AGTA05001305.1'	35298	Standard
# #		'AGTA05001306.1'	35245	Standard
# #		'AGTA05001308.1'	35219	Standard
# #		'AGTA05001307.1'	35206	Standard
# #		'AGTA05001310.1'	35122	Standard
# #		'AGTA05001309.1'	35112	Standard
# #		'AGTA05001311.1'	35084	Standard
# #		'AGTA05001312.1'	35037	Standard
# #		'AGTA05001313.1'	34997	Standard
# #		'AGTA05001314.1'	34981	Standard
# #		'AGTA05001315.1'	34931	Standard
# #		'AGTA05001316.1'	34885	Standard
# #		'AGTA05001317.1'	34844	Standard
# #		'AGTA05001321.1'	34797	Standard
# #		'AGTA05001319.1'	34762	Standard
# #		'AGTA05001318.1'	34753	Standard
# #		'AGTA05001320.1'	34748	Standard
# #		'AGTA05001322.1'	34689	Standard
# #		'AGTA05001324.1'	34625	Standard
# #		'AGTA05001323.1'	34614	Standard
# #		'AGTA05001325.1'	34594	Standard
# #		'AGTA05001326.1'	34584	Standard
# #		'AGTA05001327.1'	34497	Standard
# #		'AGTA05001328.1'	34394	Standard
# #		'AGTA05001329.1'	34368	Standard
# #		'AGTA05001330.1'	34216	Standard
# #		'AGTA05001332.1'	34142	Standard
# #		'AGTA05001331.1'	34129	Standard
# #		'AGTA05001334.1'	34034	Standard
# #		'AGTA05001333.1'	34011	Standard
# #		'AGTA05001335.1'	34009	Standard
# #		'AGTA05001336.1'	33925	Standard
# #		'AGTA05001337.1'	33863	Standard
# #		'AGTA05001338.1'	33772	Standard
# #		'AGTA05001339.1'	33733	Standard
# #		'AGTA05001340.1'	33676	Standard
# #		'AGTA05001341.1'	33674	Standard
# #		'AGTA05001344.1'	33577	Standard
# #		'AGTA05001345.1'	33554	Standard
# #		'AGTA05001342.1'	33535	Standard
# #		'AGTA05001343.1'	33529	Standard
# #		'AGTA05001346.1'	33440	Standard
# #		'AGTA05001347.1'	33424	Standard
# #		'AGTA05001348.1'	33416	Standard
# #		'AGTA05001349.1'	33334	Standard
# #		'AGTA05001350.1'	33276	Standard
# #		'AGTA05001351.1'	33186	Standard
# #		'AGTA05001352.1'	33146	Standard
# #		'AGTA05001353.1'	33077	Standard
# #		'AGTA05001355.1'	33014	Standard
# #		'AGTA05001354.1'	33001	Standard
# #		'AGTA05001356.1'	32968	Standard
# #		'AGTA05001357.1'	32954	Standard
# #		'AGTA05001358.1'	32866	Standard
# #		'AGTA05001359.1'	32836	Standard
# #		'AGTA05001360.1'	32779	Standard
# #		'AGTA05001362.1'	32626	Standard
# #		'AGTA05001361.1'	32621	Standard
# #		'AGTA05001364.1'	32601	Standard
# #		'AGTA05001363.1'	32592	Standard
# #		'AGTA05001365.1'	32573	Standard
# #		'AGTA05001366.1'	32542	Standard
# #		'AGTA05001367.1'	32530	Standard
# #		'AGTA05001368.1'	32465	Standard
# #		'AGTA05001369.1'	32349	Standard
# #		'AGTA05001370.1'	32278	Standard
# #		'AGTA05001371.1'	32252	Standard
# #		'AGTA05001372.1'	32218	Standard
# #		'AGTA05001373.1'	32174	Standard
# #		'AGTA05001374.1'	32108	Standard
# #		'AGTA05001376.1'	32048	Standard
# #		'AGTA05001375.1'	32038	Standard
# #		'AGTA05001380.1'	31996	Standard
# #		'AGTA05001379.1'	31993	Standard
# #		'AGTA05001377.1'	31985	Standard
# #		'AGTA05001378.1'	31982	Standard
# #		'AGTA05001381.1'	31911	Standard
# #		'AGTA05001382.1'	31815	Standard
# #		'AGTA05001383.1'	31807	Standard
# #		'AGTA05001384.1'	31782	Standard
# #		'AGTA05001385.1'	31721	Standard
# #		'AGTA05001386.1'	31681	Standard
# #		'AGTA05001387.1'	31658	Standard
# #		'AGTA05001389.1'	31646	Standard
# #		'AGTA05001388.1'	31636	Standard
# #		'AGTA05001390.1'	31540	Standard
# #		'AGTA05001391.1'	31496	Standard
# #		'AGTA05001392.1'	31488	Standard
# #		'AGTA05001393.1'	31472	Standard
# #		'AGTA05001394.1'	31406	Standard
# #		'AGTA05001396.1'	31329	Standard
# #		'AGTA05001395.1'	31327	Standard
# #		'AGTA05001398.1'	31272	Standard
# #		'AGTA05001397.1'	31261	Standard
# #		'AGTA05001399.1'	31194	Standard
# #		'AGTA05001401.1'	31193	Standard
# #		'AGTA05001400.1'	31152	Standard
# #		'AGTA05001402.1'	31064	Standard
# #		'AGTA05001403.1'	31047	Standard
# #		'AGTA05001405.1'	31018	Standard
# #		'AGTA05001404.1'	31011	Standard
# #		'AGTA05001406.1'	30930	Standard
# #		'AGTA05001408.1'	30894	Standard
# #		'AGTA05001407.1'	30891	Standard
# #		'AGTA05001409.1'	30799	Standard
# #		'AGTA05001411.1'	30657	Standard
# #		'AGTA05001410.1'	30652	Standard
# #		'AGTA05001412.1'	30633	Standard
# #		'AGTA05001413.1'	30617	Standard
# #		'AGTA05001414.1'	30552	Standard
# #		'AGTA05001415.1'	30506	Standard
# #		'AGTA05001416.1'	30493	Standard
# #		'AGTA05001417.1'	30448	Standard
# #		'AGTA05001418.1'	30331	Standard
# #		'AGTA05001419.1'	30279	Standard
# #		'AGTA05001420.1'	30256	Standard
# #		'AGTA05001421.1'	30172	Standard
# #		'AGTA05001422.1'	30116	Standard
# #		'AGTA05001424.1'	30099	Standard
# #		'AGTA05001423.1'	30090	Standard
# #		'AGTA05001425.1'	30033	Standard
# #		'AGTA05001430.1'	30024	Standard
# #		'AGTA05001426.1'	29986	Standard
# #		'AGTA05001427.1'	29978	Standard
# #		'AGTA05001428.1'	29960	Standard
# #		'AGTA05001429.1'	29957	Standard
# #		'AGTA05001431.1'	29885	Standard
# #		'AGTA05001432.1'	29882	Standard
# #		'AGTA05001433.1'	29809	Standard
# #		'AGTA05001435.1'	29792	Standard
# #		'AGTA05001434.1'	29783	Standard
# #		'AGTA05001437.1'	29711	Standard
# #		'AGTA05001436.1'	29707	Standard
# #		'AGTA05001438.1'	29693	Standard
# #		'AGTA05001439.1'	29652	Standard
# #		'AGTA05001440.1'	29611	Standard
# #		'AGTA05001443.1'	29567	Standard
# #		'AGTA05001441.1'	29544	Standard
# #		'AGTA05001442.1'	29542	Standard
# #		'AGTA05001444.1'	29465	Standard
# #		'AGTA05001445.1'	29381	Standard
# #		'AGTA05001446.1'	29377	Standard
# #		'AGTA05001447.1'	29334	Standard
# #		'AGTA05001448.1'	29318	Standard
# #		'AGTA05001449.1'	29307	Standard
# #		'AGTA05001450.1'	29304	Standard
# #		'AGTA05001451.1'	29229	Standard
# #		'AGTA05001453.1'	29188	Standard
# #		'AGTA05001452.1'	29182	Standard
# #		'AGTA05001454.1'	29125	Standard
# #		'AGTA05001455.1'	29100	Standard
# #		'AGTA05001456.1'	29074	Standard
# #		'AGTA05001457.1'	29072	Standard
# #		'AGTA05001458.1'	29061	Standard
# #		'AGTA05001460.1'	29019	Standard
# #		'AGTA05001459.1'	29017	Standard
# #		'AGTA05001461.1'	28998	Standard
# #		'AGTA05001462.1'	28981	Standard
# #		'AGTA05001463.1'	28976	Standard
# #		'AGTA05001465.1'	28868	Standard
# #		'AGTA05001464.1'	28858	Standard
# #		'AGTA05001466.1'	28758	Standard
# #		'AGTA05001467.1'	28684	Standard
# #		'AGTA05001468.1'	28651	Standard
# #		'AGTA05001469.1'	28608	Standard
# #		'AGTA05001470.1'	28486	Standard
# #		'AGTA05001471.1'	28481	Standard
# #		'AGTA05001473.1'	28473	Standard
# #		'AGTA05001472.1'	28446	Standard
# #		'AGTA05001474.1'	28368	Standard
# #		'AGTA05001477.1'	28355	Standard
# #		'AGTA05001475.1'	28341	Standard
# #		'AGTA05001476.1'	28335	Standard
# #		'AGTA05001478.1'	28277	Standard
# #		'AGTA05001479.1'	28216	Standard
# #		'AGTA05001480.1'	28196	Standard
# #		'AGTA05001481.1'	28096	Standard
# #		'AGTA05001482.1'	28056	Standard
# #		'AGTA05001483.1'	28028	Standard
# #		'AGTA05001484.1'	27954	Standard
# #		'AGTA05001486.1'	27913	Standard
# #		'AGTA05001485.1'	27882	Standard
# #		'AGTA05001487.1'	27857	Standard
# #		'AGTA05001489.1'	27792	Standard
# #		'AGTA05001490.1'	27776	Standard
# #		'AGTA05001488.1'	27773	Standard
# #		'AGTA05001491.1'	27692	Standard
# #		'AGTA05001492.1'	27626	Standard
# #		'AGTA05001493.1'	27512	Standard
# #		'AGTA05001495.1'	27508	Standard
# #		'AGTA05001494.1'	27507	Standard
# #		'AGTA05001496.1'	27406	Standard
# #		'AGTA05001497.1'	27374	Standard
# #		'AGTA05001500.1'	27365	Standard
# #		'AGTA05001498.1'	27312	Standard
# #		'AGTA05001501.1'	27306	Standard
# #		'AGTA05001499.1'	27279	Standard
# #		'AGTA05001504.1'	27265	Standard
# #		'AGTA05001505.1'	27259	Standard
# #		'AGTA05001503.1'	27252	Standard
# #		'AGTA05001502.1'	27241	Standard
# #		'AGTA05001506.1'	27241	Standard
# #		'AGTA05001507.1'	27223	Standard
# #		'AGTA05001508.1'	27214	Standard
# #		'AGTA05001509.1'	27186	Standard
# #		'AGTA05001510.1'	27157	Standard
# #		'AGTA05001511.1'	27130	Standard
# #		'AGTA05001512.1'	27055	Standard
# #		'AGTA05001515.1'	27041	Standard
# #		'AGTA05001513.1'	27029	Standard
# #		'AGTA05001514.1'	27027	Standard
# #		'AGTA05001516.1'	26973	Standard
# #		'AGTA05001517.1'	26965	Standard
# #		'AGTA05001518.1'	26952	Standard
# #		'AGTA05001519.1'	26879	Standard
# #		'AGTA05001520.1'	26745	Standard
# #		'AGTA05001521.1'	26691	Standard
# #		'AGTA05001522.1'	26683	Standard
# #		'AGTA05001523.1'	26664	Standard
# #		'AGTA05001524.1'	26651	Standard
# #		'AGTA05001525.1'	26631	Standard
# #		'AGTA05001526.1'	26587	Standard
# #		'AGTA05001527.1'	26537	Standard
# #		'AGTA05001528.1'	26512	Standard
# #		'AGTA05001529.1'	26487	Standard
# #		'AGTA05001530.1'	26424	Standard
# #		'AGTA05001532.1'	26354	Standard
# #		'AGTA05001531.1'	26353	Standard
# #		'AGTA05001533.1'	26321	Standard
# #		'AGTA05001534.1'	26133	Standard
# #		'AGTA05001535.1'	26104	Standard
# #		'AGTA05001536.1'	26089	Standard
# #		'AGTA05001537.1'	26056	Standard
# #		'AGTA05001538.1'	25973	Standard
# #		'AGTA05001542.1'	25956	Standard
# #		'AGTA05001543.1'	25948	Standard
# #		'AGTA05001540.1'	25940	Standard
# #		'AGTA05001539.1'	25927	Standard
# #		'AGTA05001541.1'	25927	Standard
# #		'AGTA05001544.1'	25886	Standard
# #		'AGTA05001545.1'	25838	Standard
# #		'AGTA05001546.1'	25798	Standard
# #		'AGTA05001547.1'	25793	Standard
# #		'AGTA05001548.1'	25766	Standard
# #		'AGTA05001549.1'	25756	Standard
# #		'AGTA05001550.1'	25732	Standard
# #		'AGTA05001551.1'	25695	Standard
# #		'AGTA05001554.1'	25694	Standard
# #		'AGTA05001552.1'	25682	Standard
# #		'AGTA05001553.1'	25672	Standard
# #		'AGTA05001557.1'	25661	Standard
# #		'AGTA05001555.1'	25657	Standard
# #		'AGTA05001558.1'	25640	Standard
# #		'AGTA05001556.1'	25635	Standard
# #		'AGTA05001560.1'	25558	Standard
# #		'AGTA05001559.1'	25553	Standard
# #		'AGTA05001561.1'	25537	Standard
# #		'AGTA05001562.1'	25513	Standard
# #		'AGTA05001563.1'	25390	Standard
# #		'AGTA05001564.1'	25355	Standard
# #		'AGTA05001565.1'	25331	Standard
# #		'AGTA05001566.1'	25303	Standard
# #		'AGTA05001567.1'	25288	Standard
# #		'AGTA05001570.1'	25244	Standard
# #		'AGTA05001568.1'	25240	Standard
# #		'AGTA05001569.1'	25236	Standard
# #		'AGTA05001571.1'	25229	Standard
# #		'AGTA05001573.1'	25132	Standard
# #		'AGTA05001572.1'	25116	Standard
# #		'AGTA05001574.1'	25040	Standard
# #		'AGTA05001575.1'	25019	Standard
# #		'AGTA05001576.1'	24992	Standard
# #		'AGTA05001577.1'	24899	Standard
# #		'AGTA05001578.1'	24893	Standard
# #		'AGTA05001581.1'	24815	Standard
# #		'AGTA05001579.1'	24806	Standard
# #		'AGTA05001580.1'	24789	Standard
# #		'AGTA05001582.1'	24712	Standard
# #		'AGTA05001583.1'	24665	Standard
# #		'AGTA05001584.1'	24636	Standard
# #		'AGTA05001585.1'	24626	Standard
# #		'AGTA05001586.1'	24572	Standard
# #		'AGTA05001587.1'	24503	Standard
# #		'AGTA05001588.1'	24503	Standard
# #		'AGTA05001590.1'	24455	Standard
# #		'AGTA05001589.1'	24438	Standard
# #		'AGTA05001591.1'	24408	Standard
# #		'AGTA05001592.1'	24314	Standard
# #		'AGTA05001593.1'	24273	Standard
# #		'AGTA05001594.1'	24263	Standard
# #		'AGTA05001595.1'	24131	Standard
# #		'AGTA05001596.1'	24064	Standard
# #		'AGTA05001597.1'	23964	Standard
# #		'AGTA05001598.1'	23912	Standard
# #		'AGTA05001599.1'	23856	Standard
# #		'AGTA05001600.1'	23842	Standard
# #		'AGTA05001601.1'	23826	Standard
# #		'AGTA05001602.1'	23792	Standard
# #		'AGTA05001603.1'	23772	Standard
# #		'AGTA05001604.1'	23756	Standard
# #		'AGTA05001605.1'	23756	Standard
# #		'AGTA05001606.1'	23700	Standard
# #		'AGTA05001607.1'	23551	Standard
# #		'AGTA05001609.1'	23518	Standard
# #		'AGTA05001608.1'	23515	Standard
# #		'AGTA05001610.1'	23465	Standard
# #		'AGTA05001611.1'	23460	Standard
# #		'AGTA05001612.1'	23256	Standard
# #		'AGTA05001613.1'	23240	Standard
# #		'AGTA05001616.1'	23209	Standard
# #		'AGTA05001615.1'	23207	Standard
# #		'AGTA05001614.1'	23198	Standard
# #		'AGTA05001617.1'	23070	Standard
# #		'AGTA05001618.1'	23027	Standard
# #		'AGTA05001619.1'	22982	Standard
# #		'AGTA05001620.1'	22845	Standard
# #		'AGTA05001621.1'	22808	Standard
# #		'AGTA05001623.1'	22704	Standard
# #		'AGTA05001624.1'	22681	Standard
# #		'AGTA05001622.1'	22665	Standard
# #		'AGTA05001626.1'	22622	Standard
# #		'AGTA05001625.1'	22621	Standard
# #		'AGTA05001627.1'	22597	Standard
# #		'AGTA05001628.1'	22547	Standard
# #		'AGTA05001629.1'	22474	Standard
# #		'AGTA05001630.1'	22345	Standard
# #		'AGTA05001631.1'	22320	Standard
# #		'AGTA05001632.1'	22301	Standard
# #		'AGTA05001633.1'	22139	Standard
# #		'AGTA05001634.1'	22021	Standard
# #		'AGTA05001636.1'	21950	Standard
# #		'AGTA05001635.1'	21939	Standard
# #		'AGTA05001637.1'	21893	Standard
# #		'AGTA05001641.1'	21886	Standard
# #		'AGTA05001639.1'	21885	Standard
# #		'AGTA05001640.1'	21872	Standard
# #		'AGTA05001638.1'	21871	Standard
# #		'AGTA05001642.1'	21750	Standard
# #		'AGTA05001644.1'	21731	Standard
# #		'AGTA05001643.1'	21721	Standard
# #		'AGTA05001645.1'	21678	Standard
# #		'AGTA05001646.1'	21676	Standard
# #		'AGTA05001647.1'	21672	Standard
# #		'AGTA05001648.1'	21620	Standard
# #		'AGTA05001649.1'	21587	Standard
# #		'AGTA05001651.1'	21524	Standard
# #		'AGTA05001650.1'	21507	Standard
# #		'AGTA05001652.1'	21384	Standard
# #		'AGTA05001655.1'	21272	Standard
# #		'AGTA05001653.1'	21266	Standard
# #		'AGTA05001654.1'	21254	Standard
# #		'AGTA05001656.1'	21224	Standard
# #		'AGTA05001659.1'	21215	Standard
# #		'AGTA05001657.1'	21177	Standard
# #		'AGTA05001658.1'	21170	Standard
# #		'AGTA05001660.1'	21163	Standard
# #		'AGTA05001661.1'	21145	Standard
# #		'AGTA05001662.1'	21086	Standard
# #		'AGTA05001663.1'	21048	Standard
# #		'AGTA05001664.1'	21042	Standard
# #		'AGTA05001665.1'	20948	Standard
# #		'AGTA05001666.1'	20846	Standard
# #		'AGTA05001667.1'	20717	Standard
# #		'AGTA05001668.1'	20708	Standard
# #		'AGTA05001669.1'	20640	Standard
# #		'AGTA05001670.1'	20565	Standard
# #		'AGTA05001671.1'	20538	Standard
# #		'AGTA05001673.1'	20537	Standard
# #		'AGTA05001672.1'	20527	Standard
# #		'AGTA05001674.1'	20468	Standard
# #		'AGTA05001675.1'	20458	Standard
# #		'AGTA05001677.1'	20443	Standard
# #		'AGTA05001676.1'	20424	Standard
# #		'AGTA05001678.1'	20379	Standard
# #		'AGTA05001679.1'	20303	Standard
# #		'AGTA05001680.1'	20297	Standard
# #		'AGTA05001681.1'	20239	Standard
# #		'AGTA05001682.1'	20207	Standard
# #		'AGTA05001683.1'	20128	Standard
# #		'AGTA05001684.1'	20127	Standard
# #		'AGTA05001685.1'	19993	Standard
# #		'AGTA05001686.1'	19934	Standard
# #		'AGTA05001687.1'	19272	Standard
# #		'AGTA05001688.1'	17223	Standard
# #		'MT'	16582	Standard
# #		'AGTA05001689.1'	7645	Standard
# #-----------------------------------------------
#
# 00:00:09 Predicting variants
# 00:00:13 	10000 variants (2402 variants per second), 9920 VCF entries
# 00:00:16 	20000 variants (2725 variants per second), 19856 VCF entries
# 00:00:19 	30000 variants (2911 variants per second), 29804 VCF entries
# 00:00:22 	40000 variants (3048 variants per second), 39740 VCF entries
# 00:00:25 	50000 variants (3118 variants per second), 49683 VCF entries
# 00:00:27 	60000 variants (3189 variants per second), 59604 VCF entries
# 00:00:30 	70000 variants (3276 variants per second), 69525 VCF entries
# 00:00:33 	80000 variants (3329 variants per second), 79469 VCF entries
# 00:00:35 	90000 variants (3356 variants per second), 89414 VCF entries
# 00:00:38 	100000 variants (3412 variants per second), 99365 VCF entries
# 00:00:41 	110000 variants (3439 variants per second), 109285 VCF entries
# 00:00:43 	120000 variants (3457 variants per second), 119211 VCF entries
# 00:00:46 	130000 variants (3458 variants per second), 129137 VCF entries
# 00:00:49 	140000 variants (3481 variants per second), 139068 VCF entries
# 00:00:51 	150000 variants (3499 variants per second), 149002 VCF entries
# 00:00:54 	160000 variants (3512 variants per second), 158936 VCF entries
# 00:00:57 	170000 variants (3524 variants per second), 168886 VCF entries
# 00:00:59 	180000 variants (3541 variants per second), 178824 VCF entries
# 00:01:02 	190000 variants (3527 variants per second), 188771 VCF entries
# 00:01:05 	200000 variants (3523 variants per second), 198702 VCF entries
# 00:01:08 	210000 variants (3514 variants per second), 208659 VCF entries
# 00:01:11 	220000 variants (3506 variants per second), 218595 VCF entries
# 00:01:14 	230000 variants (3497 variants per second), 228538 VCF entries
# 00:01:17 	240000 variants (3505 variants per second), 238480 VCF entries
# 00:01:20 	250000 variants (3517 variants per second), 248422 VCF entries
# 00:01:22 	260000 variants (3525 variants per second), 258352 VCF entries
# 00:01:25 	270000 variants (3514 variants per second), 268285 VCF entries
# 00:01:28 	280000 variants (3509 variants per second), 278224 VCF entries
# 00:01:31 	290000 variants (3517 variants per second), 288162 VCF entries
# 00:01:34 	300000 variants (3516 variants per second), 298091 VCF entries
# 00:01:37 	310000 variants (3516 variants per second), 308015 VCF entries
# 00:01:40 	320000 variants (3519 variants per second), 317954 VCF entries
# 00:01:42 	330000 variants (3522 variants per second), 327892 VCF entries
# 00:01:46 	340000 variants (3508 variants per second), 337839 VCF entries
# 00:01:48 	350000 variants (3505 variants per second), 347777 VCF entries
# 00:01:51 	360000 variants (3505 variants per second), 357720 VCF entries
# 00:01:54 	370000 variants (3501 variants per second), 367639 VCF entries
# 00:01:57 	380000 variants (3507 variants per second), 377577 VCF entries
# 00:02:00 	390000 variants (3503 variants per second), 387462 VCF entries
# 00:02:03 	400000 variants (3508 variants per second), 397406 VCF entries
# 00:02:05 	410000 variants (3509 variants per second), 407345 VCF entries
# 00:02:08 	420000 variants (3508 variants per second), 417266 VCF entries
# 00:02:11 	430000 variants (3500 variants per second), 427189 VCF entries
# 00:02:14 	440000 variants (3502 variants per second), 437135 VCF entries
# 00:02:17 	450000 variants (3504 variants per second), 447083 VCF entries
# 00:02:20 	460000 variants (3509 variants per second), 457021 VCF entries
# 00:02:22 	470000 variants (3518 variants per second), 466949 VCF entries
# 00:02:25 	480000 variants (3524 variants per second), 476896 VCF entries
# 00:02:27 	490000 variants (3528 variants per second), 486838 VCF entries
# 00:02:30 	500000 variants (3527 variants per second), 496779 VCF entries
# 00:02:33 	510000 variants (3524 variants per second), 506700 VCF entries
# 00:02:37 	520000 variants (3510 variants per second), 516626 VCF entries
# 00:02:40 	530000 variants (3507 variants per second), 526565 VCF entries
# 00:02:43 	540000 variants (3505 variants per second), 536491 VCF entries
# 00:02:45 	550000 variants (3508 variants per second), 546429 VCF entries
# 00:02:48 	560000 variants (3508 variants per second), 556352 VCF entries
# 00:02:52 	570000 variants (3497 variants per second), 566268 VCF entries
# 00:02:55 	580000 variants (3494 variants per second), 576210 VCF entries
# 00:02:58 	590000 variants (3488 variants per second), 586099 VCF entries
# 00:03:01 	600000 variants (3481 variants per second), 596036 VCF entries
# 00:03:04 	610000 variants (3483 variants per second), 605965 VCF entries
# 00:03:07 	620000 variants (3478 variants per second), 615891 VCF entries
# 00:03:10 	630000 variants (3473 variants per second), 625823 VCF entries
# 00:03:13 	640000 variants (3468 variants per second), 635759 VCF entries
# 00:03:16 	650000 variants (3469 variants per second), 645679 VCF entries
# 00:03:19 	660000 variants (3463 variants per second), 655611 VCF entries
# 00:03:22 	670000 variants (3460 variants per second), 665538 VCF entries
# 00:03:25 	680000 variants (3455 variants per second), 675457 VCF entries
# 00:03:28 	690000 variants (3457 variants per second), 685318 VCF entries
# 00:03:31 	700000 variants (3456 variants per second), 695228 VCF entries
# 00:03:34 	710000 variants (3456 variants per second), 705131 VCF entries
# 00:03:37 	720000 variants (3460 variants per second), 714926 VCF entries
# 00:03:40 	730000 variants (3455 variants per second), 724825 VCF entries
# 00:03:43 	740000 variants (3459 variants per second), 734739 VCF entries
# 00:03:45 	750000 variants (3460 variants per second), 744646 VCF entries
# 00:03:48 	760000 variants (3459 variants per second), 754511 VCF entries
# 00:03:51 	770000 variants (3460 variants per second), 764409 VCF entries
# 00:03:54 	780000 variants (3460 variants per second), 774314 VCF entries
# 00:03:57 	790000 variants (3461 variants per second), 784158 VCF entries
# 00:04:00 	800000 variants (3459 variants per second), 793967 VCF entries
# 00:04:03 	810000 variants (3460 variants per second), 803857 VCF entries
# 00:04:06 	820000 variants (3459 variants per second), 813708 VCF entries
# 00:04:09 	830000 variants (3454 variants per second), 823605 VCF entries
# 00:04:12 	840000 variants (3452 variants per second), 833522 VCF entries
# 00:04:14 	850000 variants (3458 variants per second), 843440 VCF entries
# 00:04:17 	860000 variants (3459 variants per second), 853369 VCF entries
# 00:04:20 	870000 variants (3458 variants per second), 863282 VCF entries
# 00:04:23 	880000 variants (3461 variants per second), 873188 VCF entries
# 00:04:26 	890000 variants (3456 variants per second), 883126 VCF entries
# 00:04:29 	900000 variants (3454 variants per second), 893062 VCF entries
# 00:04:32 	910000 variants (3456 variants per second), 902900 VCF entries
# 00:04:35 	920000 variants (3455 variants per second), 912825 VCF entries
# 00:04:38 	930000 variants (3453 variants per second), 922759 VCF entries
# 00:04:41 	940000 variants (3455 variants per second), 932682 VCF entries
# 00:04:44 	950000 variants (3455 variants per second), 942570 VCF entries
# 00:04:47 	960000 variants (3451 variants per second), 952446 VCF entries
# 00:04:49 	970000 variants (3454 variants per second), 962385 VCF entries
# 00:04:52 	980000 variants (3454 variants per second), 972306 VCF entries
# 00:04:55 	990000 variants (3451 variants per second), 982207 VCF entries
# 00:04:58 	1000000 variants (3450 variants per second), 992142 VCF entries
# 00:05:01 	1010000 variants (3451 variants per second), 1002073 VCF entries
# 00:05:05 	1020000 variants (3446 variants per second), 1012022 VCF entries
# 00:05:08 	1030000 variants (3444 variants per second), 1021953 VCF entries
# 00:05:11 	1040000 variants (3442 variants per second), 1031890 VCF entries
# 00:05:14 	1050000 variants (3442 variants per second), 1041828 VCF entries
# 00:05:17 	1060000 variants (3439 variants per second), 1051778 VCF entries
# 00:05:20 	1070000 variants (3438 variants per second), 1061715 VCF entries
# 00:05:23 	1080000 variants (3439 variants per second), 1071680 VCF entries
# 00:05:26 	1090000 variants (3439 variants per second), 1081616 VCF entries
# 00:05:29 	1100000 variants (3434 variants per second), 1091493 VCF entries
# 00:05:32 	1110000 variants (3432 variants per second), 1101396 VCF entries
# 00:05:35 	1120000 variants (3431 variants per second), 1111321 VCF entries
# 00:05:38 	1130000 variants (3430 variants per second), 1121257 VCF entries
# 00:05:41 	1140000 variants (3429 variants per second), 1131197 VCF entries
# 00:05:44 	1150000 variants (3429 variants per second), 1141127 VCF entries
# 00:05:47 	1160000 variants (3427 variants per second), 1151044 VCF entries
# 00:05:50 	1170000 variants (3424 variants per second), 1160976 VCF entries
# 00:05:54 	1180000 variants (3420 variants per second), 1170897 VCF entries
# 00:05:56 	1190000 variants (3422 variants per second), 1180833 VCF entries
# 00:05:59 	1200000 variants (3421 variants per second), 1190790 VCF entries
# 00:06:02 	1210000 variants (3419 variants per second), 1200728 VCF entries
# 00:06:05 	1220000 variants (3420 variants per second), 1210649 VCF entries
# 00:06:08 	1230000 variants (3419 variants per second), 1220581 VCF entries
# 00:06:11 	1240000 variants (3421 variants per second), 1230514 VCF entries
# 00:06:14 	1250000 variants (3423 variants per second), 1240445 VCF entries
# 00:06:17 	1260000 variants (3422 variants per second), 1250392 VCF entries
# 00:06:20 	1270000 variants (3421 variants per second), 1260295 VCF entries
# 00:06:23 	1280000 variants (3418 variants per second), 1270242 VCF entries
# 00:06:26 	1290000 variants (3417 variants per second), 1280172 VCF entries
# 00:06:29 	1300000 variants (3414 variants per second), 1290124 VCF entries
# 00:06:32 	1310000 variants (3415 variants per second), 1300066 VCF entries
# 00:06:35 	1320000 variants (3417 variants per second), 1310016 VCF entries
# 00:06:38 	1330000 variants (3417 variants per second), 1319962 VCF entries
# 00:06:41 	1340000 variants (3416 variants per second), 1329897 VCF entries
# 00:06:44 	1350000 variants (3416 variants per second), 1339835 VCF entries
# 00:06:47 	1360000 variants (3415 variants per second), 1349764 VCF entries
# 00:06:50 	1370000 variants (3412 variants per second), 1359696 VCF entries
# 00:06:53 	1380000 variants (3411 variants per second), 1369618 VCF entries
# 00:06:56 	1390000 variants (3410 variants per second), 1379561 VCF entries
# 00:06:59 	1400000 variants (3409 variants per second), 1389497 VCF entries
# 00:07:02 	1410000 variants (3408 variants per second), 1399414 VCF entries
# 00:07:05 	1420000 variants (3407 variants per second), 1409367 VCF entries
# 00:07:08 	1430000 variants (3406 variants per second), 1419310 VCF entries
# 00:07:11 	1440000 variants (3406 variants per second), 1429221 VCF entries
# 00:07:15 	1450000 variants (3402 variants per second), 1439160 VCF entries
# 00:07:18 	1460000 variants (3399 variants per second), 1449094 VCF entries
# 00:07:21 	1470000 variants (3396 variants per second), 1459021 VCF entries
# 00:07:25 	1480000 variants (3394 variants per second), 1468955 VCF entries
# 00:07:28 	1490000 variants (3392 variants per second), 1478881 VCF entries
# 00:07:31 	1500000 variants (3393 variants per second), 1488818 VCF entries
# 00:07:33 	1510000 variants (3396 variants per second), 1498749 VCF entries
# 00:07:36 	1520000 variants (3396 variants per second), 1508694 VCF entries
# 00:07:39 	1530000 variants (3395 variants per second), 1518626 VCF entries
# 00:07:42 	1540000 variants (3395 variants per second), 1528568 VCF entries
# 00:07:45 	1550000 variants (3396 variants per second), 1538509 VCF entries
# 00:07:48 	1560000 variants (3396 variants per second), 1548456 VCF entries
# 00:07:51 	1570000 variants (3395 variants per second), 1558400 VCF entries
# 00:07:54 	1580000 variants (3397 variants per second), 1568335 VCF entries
# 00:07:57 	1590000 variants (3397 variants per second), 1578266 VCF entries
# 00:08:00 	1600000 variants (3396 variants per second), 1588181 VCF entries
# 00:08:02 	1610000 variants (3399 variants per second), 1598105 VCF entries
# 00:08:05 	1620000 variants (3400 variants per second), 1608046 VCF entries
# 00:08:08 	1630000 variants (3399 variants per second), 1617959 VCF entries
# 00:08:11 	1640000 variants (3398 variants per second), 1627895 VCF entries
# 00:08:14 	1650000 variants (3397 variants per second), 1637827 VCF entries
# 00:08:17 	1660000 variants (3397 variants per second), 1647753 VCF entries
# 00:08:20 	1670000 variants (3400 variants per second), 1657704 VCF entries
# 00:08:23 	1680000 variants (3400 variants per second), 1667640 VCF entries
# 00:08:26 	1690000 variants (3400 variants per second), 1677570 VCF entries
# 00:08:29 	1700000 variants (3400 variants per second), 1687495 VCF entries
# 00:08:31 	1710000 variants (3401 variants per second), 1697427 VCF entries
# 00:08:35 	1720000 variants (3400 variants per second), 1707362 VCF entries
# 00:08:37 	1730000 variants (3399 variants per second), 1717295 VCF entries
# 00:08:40 	1740000 variants (3399 variants per second), 1727236 VCF entries
# 00:08:43 	1750000 variants (3400 variants per second), 1737188 VCF entries
# 00:08:47 	1760000 variants (3398 variants per second), 1747126 VCF entries
# 00:08:50 	1770000 variants (3398 variants per second), 1757061 VCF entries
# 00:08:53 	1780000 variants (3396 variants per second), 1766993 VCF entries
# 00:08:55 	1790000 variants (3398 variants per second), 1776930 VCF entries
# 00:08:58 	1800000 variants (3397 variants per second), 1786867 VCF entries
# 00:09:02 	1810000 variants (3394 variants per second), 1796791 VCF entries
# 00:09:05 	1820000 variants (3393 variants per second), 1806723 VCF entries
# 00:09:08 	1830000 variants (3392 variants per second), 1816623 VCF entries
# 00:09:11 	1840000 variants (3393 variants per second), 1826556 VCF entries
# 00:09:14 	1850000 variants (3393 variants per second), 1836468 VCF entries
# 00:09:17 	1860000 variants (3394 variants per second), 1846380 VCF entries
# 00:09:20 	1870000 variants (3392 variants per second), 1856297 VCF entries
# 00:09:23 	1880000 variants (3390 variants per second), 1866233 VCF entries
# 00:09:26 	1890000 variants (3388 variants per second), 1876120 VCF entries
# 00:09:30 	1900000 variants (3386 variants per second), 1886059 VCF entries
# 00:09:33 	1910000 variants (3385 variants per second), 1896003 VCF entries
# 00:09:36 	1920000 variants (3384 variants per second), 1905935 VCF entries
# 00:09:39 	1930000 variants (3385 variants per second), 1915883 VCF entries
# 00:09:42 	1940000 variants (3384 variants per second), 1925828 VCF entries
# 00:09:45 	1950000 variants (3384 variants per second), 1935757 VCF entries
# 00:09:48 	1960000 variants (3383 variants per second), 1945677 VCF entries
# 00:09:51 	1970000 variants (3384 variants per second), 1955621 VCF entries
# 00:09:54 	1980000 variants (3381 variants per second), 1965509 VCF entries
# 00:09:57 	1990000 variants (3379 variants per second), 1975426 VCF entries
# 00:10:00 	2000000 variants (3379 variants per second), 1985368 VCF entries
# 00:10:03 	2010000 variants (3379 variants per second), 1995318 VCF entries
# 00:10:07 	2020000 variants (3375 variants per second), 2005235 VCF entries
# 00:10:10 	2030000 variants (3375 variants per second), 2015169 VCF entries
# 00:10:13 	2040000 variants (3375 variants per second), 2025082 VCF entries
# 00:10:16 	2050000 variants (3375 variants per second), 2035024 VCF entries
# 00:10:19 	2060000 variants (3376 variants per second), 2044962 VCF entries
# 00:10:22 	2070000 variants (3375 variants per second), 2054901 VCF entries
# 00:10:25 	2080000 variants (3374 variants per second), 2064831 VCF entries
# 00:10:28 	2090000 variants (3373 variants per second), 2074775 VCF entries
# 00:10:31 	2100000 variants (3372 variants per second), 2084722 VCF entries
# 00:10:35 	2110000 variants (3368 variants per second), 2094665 VCF entries
# 00:10:38 	2120000 variants (3367 variants per second), 2104601 VCF entries
# 00:10:42 	2130000 variants (3364 variants per second), 2114535 VCF entries
# 00:10:45 	2140000 variants (3362 variants per second), 2124453 VCF entries
# 00:10:48 	2150000 variants (3362 variants per second), 2134377 VCF entries
# 00:10:51 	2160000 variants (3361 variants per second), 2144325 VCF entries
# 00:10:54 	2170000 variants (3360 variants per second), 2154271 VCF entries
# 00:10:57 	2180000 variants (3360 variants per second), 2164209 VCF entries
# 00:11:00 	2190000 variants (3361 variants per second), 2174142 VCF entries
# 00:11:03 	2200000 variants (3360 variants per second), 2184097 VCF entries
# 00:11:07 	2210000 variants (3358 variants per second), 2194049 VCF entries
# 00:11:09 	2220000 variants (3360 variants per second), 2203996 VCF entries
# 00:11:12 	2230000 variants (3359 variants per second), 2213936 VCF entries
# 00:11:15 	2240000 variants (3359 variants per second), 2223870 VCF entries
# 00:11:19 	2250000 variants (3358 variants per second), 2233806 VCF entries
# 00:11:22 	2260000 variants (3357 variants per second), 2243760 VCF entries
# 00:11:25 	2270000 variants (3356 variants per second), 2253697 VCF entries
# 00:11:28 	2280000 variants (3357 variants per second), 2263640 VCF entries
# 00:11:31 	2290000 variants (3357 variants per second), 2273564 VCF entries
# 00:11:34 	2300000 variants (3357 variants per second), 2283500 VCF entries
# 00:11:37 	2310000 variants (3357 variants per second), 2293433 VCF entries
# 00:11:39 	2320000 variants (3358 variants per second), 2303372 VCF entries
# 00:11:42 	2330000 variants (3359 variants per second), 2313318 VCF entries
# 00:11:45 	2340000 variants (3360 variants per second), 2323264 VCF entries
# 00:11:48 	2350000 variants (3359 variants per second), 2333180 VCF entries
# 00:11:51 	2360000 variants (3359 variants per second), 2343128 VCF entries
# 00:11:54 	2370000 variants (3358 variants per second), 2353066 VCF entries
# 00:11:57 	2380000 variants (3358 variants per second), 2363000 VCF entries
# 00:12:00 	2390000 variants (3358 variants per second), 2372931 VCF entries
# 00:12:03 	2400000 variants (3357 variants per second), 2382874 VCF entries
# 00:12:06 	2410000 variants (3357 variants per second), 2392817 VCF entries
# 00:12:09 	2420000 variants (3357 variants per second), 2402756 VCF entries
# 00:12:12 	2430000 variants (3358 variants per second), 2412682 VCF entries
# 00:12:15 	2440000 variants (3357 variants per second), 2422631 VCF entries
# 00:12:18 	2450000 variants (3356 variants per second), 2432565 VCF entries
# 00:12:22 	2460000 variants (3356 variants per second), 2442501 VCF entries
# 00:12:24 	2470000 variants (3356 variants per second), 2452430 VCF entries
# 00:12:27 	2480000 variants (3356 variants per second), 2462384 VCF entries
# 00:12:31 	2490000 variants (3356 variants per second), 2472315 VCF entries
# 00:12:33 	2500000 variants (3356 variants per second), 2482250 VCF entries
# 00:12:37 	2510000 variants (3356 variants per second), 2492179 VCF entries
# 00:12:40 	2520000 variants (3355 variants per second), 2502113 VCF entries
# 00:12:43 	2530000 variants (3355 variants per second), 2512035 VCF entries
# 00:12:46 	2540000 variants (3355 variants per second), 2521969 VCF entries
# 00:12:49 	2550000 variants (3355 variants per second), 2531903 VCF entries
# 00:12:51 	2560000 variants (3356 variants per second), 2541822 VCF entries
# 00:12:55 	2570000 variants (3354 variants per second), 2551755 VCF entries
# 00:12:58 	2580000 variants (3354 variants per second), 2561694 VCF entries
# 00:13:01 	2590000 variants (3354 variants per second), 2571620 VCF entries
# 00:13:04 	2600000 variants (3354 variants per second), 2581547 VCF entries
# 00:13:07 	2610000 variants (3355 variants per second), 2591480 VCF entries
# 00:13:09 	2620000 variants (3355 variants per second), 2601423 VCF entries
# 00:13:12 	2630000 variants (3356 variants per second), 2611356 VCF entries
# 00:13:15 	2640000 variants (3355 variants per second), 2621300 VCF entries
# 00:13:18 	2650000 variants (3355 variants per second), 2631235 VCF entries
# 00:13:21 	2660000 variants (3355 variants per second), 2641176 VCF entries
# 00:13:25 	2670000 variants (3354 variants per second), 2651118 VCF entries
# 00:13:28 	2680000 variants (3353 variants per second), 2661057 VCF entries
# 00:13:31 	2690000 variants (3352 variants per second), 2670969 VCF entries
# 00:13:34 	2700000 variants (3351 variants per second), 2680875 VCF entries
# 00:13:37 	2710000 variants (3350 variants per second), 2690804 VCF entries
# 00:13:40 	2720000 variants (3350 variants per second), 2700759 VCF entries
# 00:13:44 	2730000 variants (3349 variants per second), 2710691 VCF entries
# 00:13:47 	2740000 variants (3349 variants per second), 2720634 VCF entries
# 00:13:50 	2750000 variants (3348 variants per second), 2730576 VCF entries
# 00:13:53 	2760000 variants (3348 variants per second), 2740500 VCF entries
# 00:13:56 	2770000 variants (3348 variants per second), 2750450 VCF entries
# 00:13:59 	2780000 variants (3348 variants per second), 2760369 VCF entries
# 00:14:02 	2790000 variants (3349 variants per second), 2770303 VCF entries
# 00:14:05 	2800000 variants (3348 variants per second), 2780248 VCF entries
# 00:14:08 	2810000 variants (3348 variants per second), 2790185 VCF entries
# 00:14:11 	2820000 variants (3346 variants per second), 2800107 VCF entries
# 00:14:14 	2830000 variants (3346 variants per second), 2810041 VCF entries
# 00:14:17 	2840000 variants (3348 variants per second), 2819975 VCF entries
# 00:14:20 	2850000 variants (3348 variants per second), 2829925 VCF entries
# 00:14:23 	2860000 variants (3348 variants per second), 2839883 VCF entries
# 00:14:26 	2870000 variants (3348 variants per second), 2849840 VCF entries
# 00:14:29 	2880000 variants (3349 variants per second), 2859782 VCF entries
# 00:14:31 	2890000 variants (3349 variants per second), 2869722 VCF entries
# 00:14:35 	2900000 variants (3349 variants per second), 2879652 VCF entries
# 00:14:37 	2910000 variants (3349 variants per second), 2889582 VCF entries
# 00:14:40 	2920000 variants (3349 variants per second), 2899522 VCF entries
# 00:14:43 	2930000 variants (3349 variants per second), 2909453 VCF entries
# 00:14:46 	2940000 variants (3350 variants per second), 2919384 VCF entries
# 00:14:49 	2950000 variants (3350 variants per second), 2929325 VCF entries
# 00:14:52 	2960000 variants (3351 variants per second), 2939256 VCF entries
# 00:14:55 	2970000 variants (3352 variants per second), 2949190 VCF entries
# 00:14:57 	2980000 variants (3352 variants per second), 2959122 VCF entries
# 00:15:00 	2990000 variants (3352 variants per second), 2969031 VCF entries
# 00:15:03 	3000000 variants (3352 variants per second), 2978967 VCF entries
# 00:15:07 	3010000 variants (3351 variants per second), 2988869 VCF entries
# 00:15:10 	3020000 variants (3349 variants per second), 2998799 VCF entries
# 00:15:13 	3030000 variants (3348 variants per second), 3008743 VCF entries
# 00:15:16 	3040000 variants (3349 variants per second), 3018700 VCF entries
# 00:15:19 	3050000 variants (3349 variants per second), 3028636 VCF entries
# 00:15:22 	3060000 variants (3349 variants per second), 3038606 VCF entries
# 00:15:25 	3070000 variants (3348 variants per second), 3048536 VCF entries
# 00:15:28 	3080000 variants (3348 variants per second), 3058476 VCF entries
# 00:15:31 	3090000 variants (3348 variants per second), 3068436 VCF entries
# 00:15:35 	3100000 variants (3347 variants per second), 3078371 VCF entries
# 00:15:38 	3110000 variants (3347 variants per second), 3088309 VCF entries
# 00:15:41 	3120000 variants (3346 variants per second), 3098253 VCF entries
# 00:15:44 	3130000 variants (3345 variants per second), 3108189 VCF entries
# 00:15:47 	3140000 variants (3345 variants per second), 3118132 VCF entries
# 00:15:50 	3150000 variants (3344 variants per second), 3128057 VCF entries
# 00:15:54 	3160000 variants (3344 variants per second), 3137991 VCF entries
# 00:15:57 	3170000 variants (3343 variants per second), 3147912 VCF entries
# 00:16:00 	3180000 variants (3342 variants per second), 3157833 VCF entries
# 00:16:03 	3190000 variants (3341 variants per second), 3167760 VCF entries
# 00:16:06 	3200000 variants (3341 variants per second), 3177706 VCF entries
# 00:16:09 	3210000 variants (3341 variants per second), 3187623 VCF entries
# 00:16:12 	3220000 variants (3341 variants per second), 3197564 VCF entries
# 00:16:15 	3230000 variants (3340 variants per second), 3207491 VCF entries
# 00:16:19 	3240000 variants (3340 variants per second), 3217432 VCF entries
# 00:16:21 	3250000 variants (3340 variants per second), 3227368 VCF entries
# 00:16:24 	3260000 variants (3341 variants per second), 3237298 VCF entries
# 00:16:27 	3270000 variants (3341 variants per second), 3247223 VCF entries
# 00:16:30 	3280000 variants (3341 variants per second), 3257165 VCF entries
# 00:16:33 	3290000 variants (3343 variants per second), 3267079 VCF entries
# 00:16:35 	3300000 variants (3344 variants per second), 3277016 VCF entries
# 00:16:38 	3310000 variants (3344 variants per second), 3286935 VCF entries
# 00:16:42 	3320000 variants (3343 variants per second), 3296848 VCF entries
# 00:16:45 	3330000 variants (3343 variants per second), 3306796 VCF entries
# 00:16:48 	3340000 variants (3343 variants per second), 3316741 VCF entries
# 00:16:51 	3350000 variants (3341 variants per second), 3326682 VCF entries
# 00:16:54 	3360000 variants (3341 variants per second), 3336612 VCF entries
# 00:16:57 	3370000 variants (3341 variants per second), 3346553 VCF entries
# 00:17:00 	3380000 variants (3342 variants per second), 3356493 VCF entries
# 00:17:03 	3390000 variants (3343 variants per second), 3366428 VCF entries
# 00:17:06 	3400000 variants (3343 variants per second), 3376342 VCF entries
# 00:17:08 	3410000 variants (3343 variants per second), 3386274 VCF entries
# 00:17:11 	3420000 variants (3343 variants per second), 3396212 VCF entries
# 00:17:14 	3430000 variants (3343 variants per second), 3406152 VCF entries
# 00:17:17 	3440000 variants (3343 variants per second), 3416100 VCF entries
# 00:17:21 	3450000 variants (3342 variants per second), 3426014 VCF entries
# 00:17:24 	3460000 variants (3341 variants per second), 3435963 VCF entries
# 00:17:27 	3470000 variants (3341 variants per second), 3445904 VCF entries
# 00:17:30 	3480000 variants (3340 variants per second), 3455855 VCF entries
# 00:17:34 	3490000 variants (3339 variants per second), 3465778 VCF entries
# 00:17:37 	3500000 variants (3337 variants per second), 3475698 VCF entries
# 00:17:40 	3510000 variants (3337 variants per second), 3485561 VCF entries
# 00:17:43 	3520000 variants (3336 variants per second), 3495485 VCF entries
# 00:17:46 	3530000 variants (3337 variants per second), 3505415 VCF entries
# 00:17:49 	3540000 variants (3337 variants per second), 3515342 VCF entries
# 00:17:52 	3550000 variants (3338 variants per second), 3525261 VCF entries
# 00:17:55 	3560000 variants (3338 variants per second), 3535157 VCF entries
# 00:17:58 	3570000 variants (3338 variants per second), 3545072 VCF entries
# 00:18:01 	3580000 variants (3339 variants per second), 3555009 VCF entries
# 00:18:04 	3590000 variants (3339 variants per second), 3564926 VCF entries
# 00:18:06 	3600000 variants (3340 variants per second), 3574852 VCF entries
# 00:18:09 	3610000 variants (3341 variants per second), 3584791 VCF entries
# 00:18:12 	3620000 variants (3340 variants per second), 3594713 VCF entries
# 00:18:15 	3630000 variants (3340 variants per second), 3604640 VCF entries
# 00:18:18 	3640000 variants (3340 variants per second), 3614568 VCF entries
# 00:18:21 	3650000 variants (3340 variants per second), 3624501 VCF entries
# 00:18:24 	3660000 variants (3341 variants per second), 3634437 VCF entries
# 00:18:27 	3670000 variants (3341 variants per second), 3644374 VCF entries
# 00:18:30 	3680000 variants (3340 variants per second), 3654308 VCF entries
# 00:18:33 	3690000 variants (3340 variants per second), 3664248 VCF entries
# 00:18:36 	3700000 variants (3341 variants per second), 3674195 VCF entries
# 00:18:39 	3710000 variants (3340 variants per second), 3684133 VCF entries
# 00:18:42 	3720000 variants (3340 variants per second), 3694068 VCF entries
# 00:18:46 	3730000 variants (3339 variants per second), 3704015 VCF entries
# 00:18:49 	3740000 variants (3338 variants per second), 3713947 VCF entries
# 00:18:52 	3750000 variants (3338 variants per second), 3723873 VCF entries
# 00:18:55 	3760000 variants (3339 variants per second), 3733816 VCF entries
# 00:18:58 	3770000 variants (3339 variants per second), 3743755 VCF entries
# 00:19:01 	3780000 variants (3337 variants per second), 3753687 VCF entries
# 00:19:04 	3790000 variants (3336 variants per second), 3763616 VCF entries
# 00:19:08 	3800000 variants (3335 variants per second), 3773534 VCF entries
# 00:19:11 	3810000 variants (3335 variants per second), 3783460 VCF entries
# 00:19:14 	3820000 variants (3336 variants per second), 3793383 VCF entries
# 00:19:17 	3830000 variants (3335 variants per second), 3803321 VCF entries
# 00:19:20 	3840000 variants (3336 variants per second), 3813258 VCF entries
# 00:19:22 	3850000 variants (3336 variants per second), 3823179 VCF entries
# 00:19:25 	3860000 variants (3336 variants per second), 3833113 VCF entries
# 00:19:29 	3870000 variants (3335 variants per second), 3843041 VCF entries
# 00:19:32 	3880000 variants (3335 variants per second), 3852983 VCF entries
# 00:19:35 	3890000 variants (3336 variants per second), 3862921 VCF entries
# 00:19:38 	3900000 variants (3336 variants per second), 3872833 VCF entries
# 00:19:41 	3910000 variants (3336 variants per second), 3882768 VCF entries
# 00:19:44 	3920000 variants (3336 variants per second), 3892686 VCF entries
# 00:19:46 	3930000 variants (3336 variants per second), 3902621 VCF entries
# 00:19:49 	3940000 variants (3337 variants per second), 3912549 VCF entries
# 00:19:52 	3950000 variants (3337 variants per second), 3922486 VCF entries
# 00:19:55 	3960000 variants (3337 variants per second), 3932416 VCF entries
# 00:19:58 	3970000 variants (3337 variants per second), 3942341 VCF entries
# 00:20:01 	3980000 variants (3337 variants per second), 3952272 VCF entries
# 00:20:04 	3990000 variants (3336 variants per second), 3962211 VCF entries
# 00:20:07 	4000000 variants (3337 variants per second), 3972102 VCF entries
# 00:20:10 	4010000 variants (3338 variants per second), 3982040 VCF entries
# 00:20:13 	4020000 variants (3337 variants per second), 3991968 VCF entries
# 00:20:16 	4030000 variants (3338 variants per second), 4001890 VCF entries
# 00:20:18 	4040000 variants (3339 variants per second), 4011820 VCF entries
# 00:20:21 	4050000 variants (3340 variants per second), 4021756 VCF entries
# 00:20:25 	4060000 variants (3339 variants per second), 4031691 VCF entries
# 00:20:27 	4070000 variants (3339 variants per second), 4041626 VCF entries
# 00:20:30 	4080000 variants (3339 variants per second), 4051550 VCF entries
# 00:20:33 	4090000 variants (3340 variants per second), 4061493 VCF entries
# 00:20:36 	4100000 variants (3340 variants per second), 4071433 VCF entries
# 00:20:39 	4110000 variants (3340 variants per second), 4081372 VCF entries
# 00:20:42 	4120000 variants (3340 variants per second), 4091312 VCF entries
# 00:20:45 	4130000 variants (3340 variants per second), 4101251 VCF entries
# 00:20:48 	4140000 variants (3340 variants per second), 4111188 VCF entries
# 00:20:51 	4150000 variants (3339 variants per second), 4121122 VCF entries
# 00:20:55 	4160000 variants (3338 variants per second), 4131052 VCF entries
# 00:20:57 	4170000 variants (3339 variants per second), 4140993 VCF entries
# 00:21:00 	4180000 variants (3339 variants per second), 4150931 VCF entries
# 00:21:03 	4190000 variants (3339 variants per second), 4160866 VCF entries
# 00:21:06 	4200000 variants (3339 variants per second), 4170803 VCF entries
# 00:21:09 	4210000 variants (3339 variants per second), 4180731 VCF entries
# 00:21:12 	4220000 variants (3339 variants per second), 4190658 VCF entries
# 00:21:15 	4230000 variants (3339 variants per second), 4200605 VCF entries
# 00:21:18 	4240000 variants (3339 variants per second), 4210557 VCF entries
# 00:21:21 	4250000 variants (3339 variants per second), 4220498 VCF entries
# 00:21:24 	4260000 variants (3339 variants per second), 4230445 VCF entries
# 00:21:27 	4270000 variants (3338 variants per second), 4240395 VCF entries
# 00:21:30 	4280000 variants (3339 variants per second), 4250330 VCF entries
# 00:21:33 	4290000 variants (3340 variants per second), 4260272 VCF entries
# 00:21:36 	4300000 variants (3340 variants per second), 4270211 VCF entries
# 00:21:39 	4310000 variants (3340 variants per second), 4280151 VCF entries
# 00:21:42 	4320000 variants (3340 variants per second), 4290097 VCF entries
# 00:21:44 	4330000 variants (3341 variants per second), 4300044 VCF entries
# 00:21:48 	4340000 variants (3341 variants per second), 4309977 VCF entries
# 00:21:50 	4350000 variants (3342 variants per second), 4319911 VCF entries
# 00:21:53 	4360000 variants (3342 variants per second), 4329848 VCF entries
# 00:21:56 	4370000 variants (3342 variants per second), 4339788 VCF entries
# 00:21:59 	4380000 variants (3343 variants per second), 4349727 VCF entries
# 00:22:02 	4390000 variants (3343 variants per second), 4359667 VCF entries
# 00:22:04 	4400000 variants (3343 variants per second), 4369603 VCF entries
# 00:22:07 	4410000 variants (3344 variants per second), 4379530 VCF entries
# 00:22:10 	4420000 variants (3344 variants per second), 4389462 VCF entries
# 00:22:13 	4430000 variants (3344 variants per second), 4399396 VCF entries
# 00:22:16 	4440000 variants (3345 variants per second), 4409320 VCF entries
# 00:22:19 	4450000 variants (3345 variants per second), 4419255 VCF entries
# 00:22:22 	4460000 variants (3345 variants per second), 4429194 VCF entries
# 00:22:25 	4470000 variants (3345 variants per second), 4439125 VCF entries
# 00:22:28 	4480000 variants (3345 variants per second), 4449058 VCF entries
# 00:22:31 	4490000 variants (3345 variants per second), 4459000 VCF entries
# 00:22:34 	4500000 variants (3344 variants per second), 4468946 VCF entries
# 00:22:37 	4510000 variants (3345 variants per second), 4478885 VCF entries
# 00:22:40 	4520000 variants (3345 variants per second), 4488824 VCF entries
# 00:22:43 	4530000 variants (3345 variants per second), 4498760 VCF entries
# 00:22:46 	4540000 variants (3345 variants per second), 4508696 VCF entries
# 00:22:49 	4550000 variants (3345 variants per second), 4518637 VCF entries
# 00:22:51 	4560000 variants (3346 variants per second), 4528572 VCF entries
# 00:22:55 	4570000 variants (3345 variants per second), 4538521 VCF entries
# 00:22:57 	4580000 variants (3345 variants per second), 4548452 VCF entries
# 00:23:00 	4590000 variants (3346 variants per second), 4558365 VCF entries
# 00:23:03 	4600000 variants (3346 variants per second), 4568305 VCF entries
# 00:23:06 	4610000 variants (3346 variants per second), 4578251 VCF entries
# 00:23:09 	4620000 variants (3346 variants per second), 4588205 VCF entries
# 00:23:12 	4630000 variants (3346 variants per second), 4598130 VCF entries
# 00:23:15 	4640000 variants (3346 variants per second), 4608062 VCF entries
# 00:23:18 	4650000 variants (3346 variants per second), 4618000 VCF entries
# 00:23:21 	4660000 variants (3347 variants per second), 4627933 VCF entries
# 00:23:24 	4670000 variants (3346 variants per second), 4637871 VCF entries
# 00:23:27 	4680000 variants (3346 variants per second), 4647807 VCF entries
# 00:23:30 	4690000 variants (3346 variants per second), 4657732 VCF entries
# 00:23:33 	4700000 variants (3346 variants per second), 4667670 VCF entries
# 00:23:36 	4710000 variants (3346 variants per second), 4677610 VCF entries
# 00:23:39 	4720000 variants (3346 variants per second), 4687539 VCF entries
# 00:23:42 	4730000 variants (3346 variants per second), 4697467 VCF entries
# 00:23:45 	4740000 variants (3346 variants per second), 4707400 VCF entries
# 00:23:48 	4750000 variants (3346 variants per second), 4717344 VCF entries
# 00:23:51 	4760000 variants (3346 variants per second), 4727285 VCF entries
# 00:23:54 	4770000 variants (3346 variants per second), 4737219 VCF entries
# 00:23:57 	4780000 variants (3347 variants per second), 4747143 VCF entries
# 00:24:00 	4790000 variants (3347 variants per second), 4757077 VCF entries
# 00:24:03 	4800000 variants (3347 variants per second), 4767005 VCF entries
# 00:24:05 	4810000 variants (3347 variants per second), 4776905 VCF entries
# 00:24:08 	4820000 variants (3348 variants per second), 4786838 VCF entries
# 00:24:11 	4830000 variants (3348 variants per second), 4796759 VCF entries
# 00:24:14 	4840000 variants (3349 variants per second), 4806705 VCF entries
# 00:24:17 	4850000 variants (3348 variants per second), 4816655 VCF entries
# 00:24:20 	4860000 variants (3348 variants per second), 4826591 VCF entries
# 00:24:23 	4870000 variants (3349 variants per second), 4836518 VCF entries
# 00:24:26 	4880000 variants (3349 variants per second), 4846436 VCF entries
# 00:24:29 	4890000 variants (3348 variants per second), 4856338 VCF entries
# 00:24:32 	4900000 variants (3349 variants per second), 4866287 VCF entries
# 00:24:34 	4910000 variants (3349 variants per second), 4876223 VCF entries
# 00:24:37 	4920000 variants (3349 variants per second), 4886170 VCF entries
# 00:24:40 	4930000 variants (3349 variants per second), 4896104 VCF entries
# 00:24:43 	4940000 variants (3350 variants per second), 4906047 VCF entries
# 00:24:46 	4950000 variants (3349 variants per second), 4915976 VCF entries
# 00:24:49 	4960000 variants (3349 variants per second), 4925913 VCF entries
# 00:24:52 	4970000 variants (3349 variants per second), 4935855 VCF entries
# 00:24:56 	4980000 variants (3348 variants per second), 4945789 VCF entries
# 00:24:59 	4990000 variants (3349 variants per second), 4955712 VCF entries
# 00:25:02 	5000000 variants (3348 variants per second), 4965653 VCF entries
# 00:25:05 	5010000 variants (3348 variants per second), 4975599 VCF entries
# 00:25:08 	5020000 variants (3348 variants per second), 4985524 VCF entries
# 00:25:11 	5030000 variants (3348 variants per second), 4995461 VCF entries
# 00:25:14 	5040000 variants (3348 variants per second), 5005399 VCF entries
# 00:25:17 	5050000 variants (3348 variants per second), 5015332 VCF entries
# 00:25:20 	5060000 variants (3348 variants per second), 5025257 VCF entries
# 00:25:23 	5070000 variants (3348 variants per second), 5035175 VCF entries
# 00:25:26 	5080000 variants (3348 variants per second), 5045082 VCF entries
# 00:25:29 	5090000 variants (3348 variants per second), 5054988 VCF entries
# 00:25:32 	5100000 variants (3348 variants per second), 5064929 VCF entries
# 00:25:35 	5110000 variants (3347 variants per second), 5074861 VCF entries
# 00:25:38 	5120000 variants (3347 variants per second), 5084789 VCF entries
# 00:25:41 	5130000 variants (3347 variants per second), 5094730 VCF entries
# 00:25:44 	5140000 variants (3348 variants per second), 5104612 VCF entries
# 00:25:47 	5150000 variants (3348 variants per second), 5114538 VCF entries
# 00:25:49 	5160000 variants (3348 variants per second), 5124477 VCF entries
# 00:25:53 	5170000 variants (3348 variants per second), 5134391 VCF entries
# 00:25:55 	5180000 variants (3348 variants per second), 5144327 VCF entries
# 00:25:58 	5190000 variants (3349 variants per second), 5154265 VCF entries
# 00:26:01 	5200000 variants (3349 variants per second), 5164202 VCF entries
# 00:26:04 	5210000 variants (3349 variants per second), 5174133 VCF entries
# 00:26:07 	5220000 variants (3350 variants per second), 5184090 VCF entries
# 00:26:10 	5230000 variants (3350 variants per second), 5194037 VCF entries
# 00:26:13 	5240000 variants (3350 variants per second), 5203949 VCF entries
# 00:26:16 	5250000 variants (3350 variants per second), 5213894 VCF entries
# 00:26:18 	5260000 variants (3350 variants per second), 5223828 VCF entries
# 00:26:21 	5270000 variants (3350 variants per second), 5233769 VCF entries
# 00:26:24 	5280000 variants (3350 variants per second), 5243709 VCF entries
# 00:26:27 	5290000 variants (3351 variants per second), 5253659 VCF entries
# 00:26:30 	5300000 variants (3351 variants per second), 5263589 VCF entries
# 00:26:33 	5310000 variants (3350 variants per second), 5273527 VCF entries
# 00:26:36 	5320000 variants (3351 variants per second), 5283469 VCF entries
# 00:26:39 	5330000 variants (3351 variants per second), 5293418 VCF entries
# 00:26:42 	5340000 variants (3352 variants per second), 5303361 VCF entries
# 00:26:45 	5350000 variants (3352 variants per second), 5313309 VCF entries
# 00:26:47 	5360000 variants (3352 variants per second), 5323253 VCF entries
# 00:26:50 	5370000 variants (3352 variants per second), 5333196 VCF entries
# 00:26:53 	5380000 variants (3352 variants per second), 5343129 VCF entries
# 00:26:57 	5390000 variants (3352 variants per second), 5353068 VCF entries
# 00:26:59 	5400000 variants (3352 variants per second), 5363015 VCF entries
# 00:27:02 	5410000 variants (3352 variants per second), 5372970 VCF entries
# 00:27:05 	5420000 variants (3352 variants per second), 5382914 VCF entries
# 00:27:08 	5430000 variants (3352 variants per second), 5392846 VCF entries
# 00:27:11 	5440000 variants (3352 variants per second), 5402778 VCF entries
# 00:27:14 	5450000 variants (3352 variants per second), 5412712 VCF entries
# 00:27:17 	5460000 variants (3352 variants per second), 5422644 VCF entries
# 00:27:20 	5470000 variants (3352 variants per second), 5432579 VCF entries
# 00:27:23 	5480000 variants (3352 variants per second), 5442508 VCF entries
# 00:27:26 	5490000 variants (3352 variants per second), 5452446 VCF entries
# 00:27:29 	5500000 variants (3352 variants per second), 5462366 VCF entries
# 00:27:32 	5510000 variants (3352 variants per second), 5472288 VCF entries
# 00:27:35 	5520000 variants (3352 variants per second), 5482226 VCF entries
# 00:27:38 	5530000 variants (3352 variants per second), 5492159 VCF entries
# 00:27:41 	5540000 variants (3352 variants per second), 5502097 VCF entries
# 00:27:44 	5550000 variants (3352 variants per second), 5512038 VCF entries
# 00:27:47 	5560000 variants (3352 variants per second), 5521956 VCF entries
# 00:27:50 	5570000 variants (3352 variants per second), 5531880 VCF entries
# 00:27:53 	5580000 variants (3351 variants per second), 5541803 VCF entries
# 00:27:56 	5590000 variants (3351 variants per second), 5551730 VCF entries
# 00:28:00 	5600000 variants (3350 variants per second), 5561626 VCF entries
# 00:28:03 	5610000 variants (3350 variants per second), 5571551 VCF entries
# 00:28:06 	5620000 variants (3350 variants per second), 5581475 VCF entries
# 00:28:10 	5630000 variants (3349 variants per second), 5591398 VCF entries
# 00:28:13 	5640000 variants (3347 variants per second), 5601340 VCF entries
# 00:28:16 	5650000 variants (3347 variants per second), 5611239 VCF entries
# 00:28:20 	5660000 variants (3347 variants per second), 5621162 VCF entries
# 00:28:23 	5670000 variants (3346 variants per second), 5631095 VCF entries
# 00:28:26 	5680000 variants (3346 variants per second), 5641017 VCF entries
# 00:28:29 	5690000 variants (3346 variants per second), 5650923 VCF entries
# 00:28:32 	5700000 variants (3346 variants per second), 5660855 VCF entries
# 00:28:35 	5710000 variants (3346 variants per second), 5670783 VCF entries
# 00:28:38 	5720000 variants (3346 variants per second), 5680737 VCF entries
# 00:28:41 	5730000 variants (3346 variants per second), 5690680 VCF entries
# 00:28:44 	5740000 variants (3345 variants per second), 5700623 VCF entries
# 00:28:47 	5750000 variants (3345 variants per second), 5710560 VCF entries
# 00:28:50 	5760000 variants (3345 variants per second), 5720505 VCF entries
# 00:28:53 	5770000 variants (3345 variants per second), 5730441 VCF entries
# 00:28:57 	5780000 variants (3345 variants per second), 5740364 VCF entries
# 00:29:00 	5790000 variants (3344 variants per second), 5750288 VCF entries
# 00:29:03 	5800000 variants (3343 variants per second), 5760217 VCF entries
# 00:29:06 	5810000 variants (3343 variants per second), 5770134 VCF entries
# 00:29:10 	5820000 variants (3343 variants per second), 5780061 VCF entries
# 00:29:13 	5830000 variants (3342 variants per second), 5789989 VCF entries
# 00:29:16 	5840000 variants (3341 variants per second), 5799912 VCF entries
# 00:29:19 	5850000 variants (3341 variants per second), 5809864 VCF entries
# 00:29:22 	5860000 variants (3341 variants per second), 5819816 VCF entries
# 00:29:25 	5870000 variants (3342 variants per second), 5829762 VCF entries
# 00:29:28 	5880000 variants (3341 variants per second), 5839692 VCF entries
# 00:29:31 	5890000 variants (3342 variants per second), 5849636 VCF entries
# 00:29:34 	5900000 variants (3342 variants per second), 5859571 VCF entries
# 00:29:37 	5910000 variants (3341 variants per second), 5869501 VCF entries
# 00:29:40 	5920000 variants (3341 variants per second), 5879424 VCF entries
# 00:29:44 	5930000 variants (3341 variants per second), 5889355 VCF entries
# 00:29:47 	5940000 variants (3340 variants per second), 5899285 VCF entries
# 00:29:50 	5950000 variants (3340 variants per second), 5909223 VCF entries
# 00:29:53 	5960000 variants (3340 variants per second), 5919131 VCF entries
# 00:29:56 	5970000 variants (3339 variants per second), 5929067 VCF entries
# 00:29:59 	5980000 variants (3339 variants per second), 5938995 VCF entries
# 00:30:02 	5990000 variants (3339 variants per second), 5948912 VCF entries
# 00:30:05 	6000000 variants (3339 variants per second), 5958850 VCF entries
# 00:30:09 	6010000 variants (3338 variants per second), 5968794 VCF entries
# 00:30:12 	6020000 variants (3337 variants per second), 5978693 VCF entries
# 00:30:16 	6030000 variants (3337 variants per second), 5988619 VCF entries
# 00:30:19 	6040000 variants (3336 variants per second), 5998546 VCF entries
# 00:30:22 	6050000 variants (3336 variants per second), 6008487 VCF entries
# 00:30:25 	6060000 variants (3336 variants per second), 6018412 VCF entries
# 00:30:28 	6070000 variants (3335 variants per second), 6028360 VCF entries
# 00:30:31 	6080000 variants (3335 variants per second), 6038274 VCF entries
# 00:30:34 	6090000 variants (3335 variants per second), 6048225 VCF entries
# 00:30:37 	6100000 variants (3335 variants per second), 6058170 VCF entries
# 00:30:41 	6110000 variants (3335 variants per second), 6068107 VCF entries
# 00:30:44 	6120000 variants (3334 variants per second), 6078050 VCF entries
# 00:30:47 	6130000 variants (3334 variants per second), 6087982 VCF entries
# 00:30:50 	6140000 variants (3334 variants per second), 6097925 VCF entries
# 00:30:53 	6150000 variants (3334 variants per second), 6107854 VCF entries
# 00:30:56 	6160000 variants (3333 variants per second), 6117787 VCF entries
# 00:31:00 	6170000 variants (3333 variants per second), 6127722 VCF entries
# 00:31:03 	6180000 variants (3333 variants per second), 6137656 VCF entries
# 00:31:06 	6190000 variants (3332 variants per second), 6147580 VCF entries
# 00:31:09 	6200000 variants (3332 variants per second), 6157513 VCF entries
# 00:31:13 	6210000 variants (3331 variants per second), 6167445 VCF entries
# 00:31:16 	6220000 variants (3330 variants per second), 6177369 VCF entries
# 00:31:19 	6230000 variants (3330 variants per second), 6187294 VCF entries
# 00:31:22 	6240000 variants (3330 variants per second), 6197232 VCF entries
# 00:31:26 	6250000 variants (3329 variants per second), 6207154 VCF entries
# 00:31:29 	6260000 variants (3329 variants per second), 6217087 VCF entries
# 00:31:32 	6270000 variants (3329 variants per second), 6227038 VCF entries
# 00:31:35 	6280000 variants (3328 variants per second), 6236972 VCF entries
# 00:31:38 	6290000 variants (3328 variants per second), 6246889 VCF entries
# 00:31:42 	6300000 variants (3328 variants per second), 6256821 VCF entries
# 00:31:45 	6310000 variants (3327 variants per second), 6266745 VCF entries
# 00:31:48 	6320000 variants (3326 variants per second), 6276671 VCF entries
# 00:31:52 	6330000 variants (3326 variants per second), 6286617 VCF entries
# 00:31:55 	6340000 variants (3325 variants per second), 6296551 VCF entries
# 00:31:58 	6350000 variants (3325 variants per second), 6306501 VCF entries
# 00:32:01 	6360000 variants (3325 variants per second), 6316434 VCF entries
# 00:32:05 	6370000 variants (3324 variants per second), 6326370 VCF entries
# 00:32:08 	6380000 variants (3324 variants per second), 6336301 VCF entries
# 00:32:11 	6390000 variants (3324 variants per second), 6346227 VCF entries
# 00:32:14 	6400000 variants (3323 variants per second), 6356142 VCF entries
# 00:32:17 	6410000 variants (3323 variants per second), 6366072 VCF entries
# 00:32:20 	6420000 variants (3323 variants per second), 6376012 VCF entries
# 00:32:24 	6430000 variants (3323 variants per second), 6385919 VCF entries
# 00:32:27 	6440000 variants (3323 variants per second), 6395845 VCF entries
# 00:32:30 	6450000 variants (3322 variants per second), 6405741 VCF entries
# 00:32:33 	6460000 variants (3322 variants per second), 6415652 VCF entries
# 00:32:36 	6470000 variants (3322 variants per second), 6425583 VCF entries
# 00:32:40 	6480000 variants (3320 variants per second), 6435456 VCF entries
# 00:32:44 	6490000 variants (3318 variants per second), 6445180 VCF entries
# 00:32:47 	6500000 variants (3318 variants per second), 6455121 VCF entries
# 00:32:50 	6510000 variants (3318 variants per second), 6465020 VCF entries
# 00:32:54 	6520000 variants (3317 variants per second), 6474943 VCF entries
# 00:32:57 	6530000 variants (3317 variants per second), 6484874 VCF entries
# 00:33:00 	6540000 variants (3317 variants per second), 6494805 VCF entries
# 00:33:03 	6550000 variants (3316 variants per second), 6504748 VCF entries
# 00:33:07 	6560000 variants (3316 variants per second), 6514680 VCF entries
# 00:33:10 	6570000 variants (3316 variants per second), 6524597 VCF entries
# 00:33:13 	6580000 variants (3316 variants per second), 6534523 VCF entries
# 00:33:16 	6590000 variants (3316 variants per second), 6544453 VCF entries
# 00:33:19 	6600000 variants (3316 variants per second), 6554393 VCF entries
# 00:33:22 	6610000 variants (3315 variants per second), 6564325 VCF entries
# 00:33:25 	6620000 variants (3315 variants per second), 6574270 VCF entries
# 00:33:29 	6630000 variants (3315 variants per second), 6584216 VCF entries
# 00:33:32 	6640000 variants (3314 variants per second), 6594143 VCF entries
# 00:33:36 	6650000 variants (3313 variants per second), 6604072 VCF entries
# 00:33:39 	6660000 variants (3313 variants per second), 6614030 VCF entries
# 00:33:41 	6670000 variants (3313 variants per second), 6623971 VCF entries
# 00:33:44 	6680000 variants (3313 variants per second), 6633914 VCF entries
# 00:33:48 	6690000 variants (3313 variants per second), 6643847 VCF entries
# 00:33:51 	6700000 variants (3313 variants per second), 6653786 VCF entries
# 00:33:54 	6710000 variants (3313 variants per second), 6663726 VCF entries
# 00:33:57 	6720000 variants (3312 variants per second), 6673657 VCF entries
# 00:34:00 	6730000 variants (3312 variants per second), 6683595 VCF entries
# 00:34:03 	6740000 variants (3312 variants per second), 6693548 VCF entries
# 00:34:06 	6750000 variants (3312 variants per second), 6703485 VCF entries
# 00:34:10 	6760000 variants (3311 variants per second), 6713439 VCF entries
# 00:34:13 	6770000 variants (3311 variants per second), 6723351 VCF entries
# 00:34:17 	6780000 variants (3310 variants per second), 6733279 VCF entries
# 00:34:20 	6790000 variants (3310 variants per second), 6743222 VCF entries
# 00:34:23 	6800000 variants (3310 variants per second), 6753164 VCF entries
# 00:34:27 	6810000 variants (3309 variants per second), 6763088 VCF entries
# 00:34:30 	6820000 variants (3308 variants per second), 6773010 VCF entries
# 00:34:33 	6830000 variants (3308 variants per second), 6782957 VCF entries
# 00:34:36 	6840000 variants (3307 variants per second), 6792902 VCF entries
# 00:34:40 	6850000 variants (3307 variants per second), 6802847 VCF entries
# 00:34:43 	6860000 variants (3306 variants per second), 6812737 VCF entries
# 00:34:46 	6870000 variants (3306 variants per second), 6822680 VCF entries
# 00:34:49 	6880000 variants (3306 variants per second), 6832606 VCF entries
# 00:34:52 	6890000 variants (3306 variants per second), 6842549 VCF entries
# 00:34:55 	6900000 variants (3306 variants per second), 6852490 VCF entries
# 00:34:58 	6910000 variants (3306 variants per second), 6862419 VCF entries
# 00:35:01 	6920000 variants (3307 variants per second), 6872347 VCF entries
# 00:35:04 	6930000 variants (3306 variants per second), 6882268 VCF entries
# 00:35:07 	6940000 variants (3306 variants per second), 6892180 VCF entries
# 00:35:11 	6950000 variants (3306 variants per second), 6902106 VCF entries
# 00:35:14 	6960000 variants (3305 variants per second), 6912043 VCF entries
# 00:35:17 	6970000 variants (3305 variants per second), 6921997 VCF entries
# 00:35:20 	6980000 variants (3306 variants per second), 6931949 VCF entries
# 00:35:23 	6990000 variants (3305 variants per second), 6941891 VCF entries
# 00:35:26 	7000000 variants (3305 variants per second), 6951844 VCF entries
# 00:35:29 	7010000 variants (3305 variants per second), 6961785 VCF entries
# 00:35:33 	7020000 variants (3305 variants per second), 6971734 VCF entries
# 00:35:36 	7030000 variants (3304 variants per second), 6981672 VCF entries
# 00:35:38 	7040000 variants (3305 variants per second), 6991615 VCF entries
# 00:35:41 	7050000 variants (3305 variants per second), 7001547 VCF entries
# 00:35:44 	7060000 variants (3305 variants per second), 7011487 VCF entries
# 00:35:47 	7070000 variants (3305 variants per second), 7021431 VCF entries
# 00:35:51 	7080000 variants (3305 variants per second), 7031341 VCF entries
# 00:35:55 	7090000 variants (3303 variants per second), 7041244 VCF entries
# 00:35:58 	7100000 variants (3303 variants per second), 7051166 VCF entries
# 00:36:01 	7110000 variants (3303 variants per second), 7061096 VCF entries
# 00:36:04 	7120000 variants (3303 variants per second), 7071024 VCF entries
# 00:36:07 	7130000 variants (3303 variants per second), 7080958 VCF entries
# 00:36:11 	7140000 variants (3302 variants per second), 7090896 VCF entries
# 00:36:14 	7150000 variants (3302 variants per second), 7100837 VCF entries
# 00:36:17 	7160000 variants (3302 variants per second), 7110775 VCF entries
# 00:36:20 	7170000 variants (3302 variants per second), 7120679 VCF entries
# 00:36:23 	7180000 variants (3302 variants per second), 7130604 VCF entries
# 00:36:26 	7190000 variants (3301 variants per second), 7140531 VCF entries
# 00:36:30 	7200000 variants (3301 variants per second), 7150485 VCF entries
# 00:36:33 	7210000 variants (3300 variants per second), 7160411 VCF entries
# 00:36:36 	7220000 variants (3300 variants per second), 7170367 VCF entries
# 00:36:39 	7230000 variants (3301 variants per second), 7180300 VCF entries
# 00:36:42 	7240000 variants (3301 variants per second), 7190245 VCF entries
# 00:36:45 	7250000 variants (3300 variants per second), 7200183 VCF entries
# 00:36:48 	7260000 variants (3300 variants per second), 7210082 VCF entries
# 00:36:51 	7270000 variants (3301 variants per second), 7220032 VCF entries
# 00:36:54 	7280000 variants (3300 variants per second), 7229969 VCF entries
# 00:36:57 	7290000 variants (3300 variants per second), 7239876 VCF entries
# 00:37:00 	7300000 variants (3300 variants per second), 7249821 VCF entries
# 00:37:03 	7310000 variants (3300 variants per second), 7259762 VCF entries
# 00:37:07 	7320000 variants (3299 variants per second), 7269694 VCF entries
# 00:37:11 	7330000 variants (3298 variants per second), 7279603 VCF entries
# 00:37:14 	7340000 variants (3298 variants per second), 7289500 VCF entries
# 00:37:17 	7350000 variants (3298 variants per second), 7299428 VCF entries
# 00:37:20 	7360000 variants (3298 variants per second), 7309352 VCF entries
# 00:37:24 	7370000 variants (3297 variants per second), 7319256 VCF entries
# 00:37:27 	7380000 variants (3297 variants per second), 7329144 VCF entries
# 00:37:30 	7390000 variants (3297 variants per second), 7339065 VCF entries
# 00:37:33 	7400000 variants (3297 variants per second), 7348970 VCF entries
# 00:37:36 	7410000 variants (3296 variants per second), 7358904 VCF entries
# 00:37:40 	7420000 variants (3296 variants per second), 7368841 VCF entries
# 00:37:43 	7430000 variants (3296 variants per second), 7378701 VCF entries
# 00:37:46 	7440000 variants (3296 variants per second), 7388626 VCF entries
# 00:37:49 	7450000 variants (3295 variants per second), 7398548 VCF entries
# 00:37:52 	7460000 variants (3295 variants per second), 7408497 VCF entries
# 00:37:55 	7470000 variants (3295 variants per second), 7418432 VCF entries
# 00:37:58 	7480000 variants (3295 variants per second), 7428364 VCF entries
# 00:38:02 	7490000 variants (3295 variants per second), 7438299 VCF entries
# 00:38:05 	7500000 variants (3295 variants per second), 7448229 VCF entries
# 00:38:08 	7510000 variants (3294 variants per second), 7458165 VCF entries
# 00:38:12 	7520000 variants (3293 variants per second), 7468072 VCF entries
# 00:38:15 	7530000 variants (3293 variants per second), 7477992 VCF entries
# 00:38:18 	7540000 variants (3293 variants per second), 7487904 VCF entries
# 00:38:21 	7550000 variants (3292 variants per second), 7497801 VCF entries
# 00:38:25 	7560000 variants (3292 variants per second), 7507700 VCF entries
# 00:38:28 	7570000 variants (3292 variants per second), 7517611 VCF entries
# 00:38:31 	7580000 variants (3291 variants per second), 7527528 VCF entries
# 00:38:35 	7590000 variants (3291 variants per second), 7537417 VCF entries
# 00:38:38 	7600000 variants (3291 variants per second), 7547335 VCF entries
# 00:38:41 	7610000 variants (3290 variants per second), 7557253 VCF entries
# 00:38:44 	7620000 variants (3290 variants per second), 7567183 VCF entries
# 00:38:48 	7630000 variants (3290 variants per second), 7577083 VCF entries
# 00:38:51 	7640000 variants (3289 variants per second), 7586963 VCF entries
# 00:38:54 	7650000 variants (3290 variants per second), 7596884 VCF entries
# 00:38:57 	7660000 variants (3289 variants per second), 7606810 VCF entries
# 00:39:00 	7670000 variants (3289 variants per second), 7616662 VCF entries
# 00:39:04 	7680000 variants (3288 variants per second), 7626586 VCF entries
# 00:39:07 	7690000 variants (3288 variants per second), 7636498 VCF entries
# 00:39:10 	7700000 variants (3288 variants per second), 7646420 VCF entries
# 00:39:14 	7710000 variants (3287 variants per second), 7656345 VCF entries
# 00:39:17 	7720000 variants (3287 variants per second), 7666250 VCF entries
# 00:39:20 	7730000 variants (3286 variants per second), 7676134 VCF entries
# 00:39:24 	7740000 variants (3286 variants per second), 7686061 VCF entries
# 00:39:27 	7750000 variants (3285 variants per second), 7695976 VCF entries
# 00:39:31 	7760000 variants (3285 variants per second), 7705883 VCF entries
# 00:39:34 	7770000 variants (3284 variants per second), 7715768 VCF entries
# 00:39:37 	7780000 variants (3284 variants per second), 7725647 VCF entries
# 00:39:40 	7790000 variants (3284 variants per second), 7735573 VCF entries
# 00:39:44 	7800000 variants (3284 variants per second), 7745467 VCF entries
# 00:39:47 	7810000 variants (3283 variants per second), 7755391 VCF entries
# 00:39:50 	7820000 variants (3283 variants per second), 7765308 VCF entries
# 00:39:53 	7830000 variants (3283 variants per second), 7775219 VCF entries
# 00:39:57 	7840000 variants (3282 variants per second), 7785146 VCF entries
# 00:40:00 	7850000 variants (3282 variants per second), 7795040 VCF entries
# 00:40:04 	7860000 variants (3281 variants per second), 7804948 VCF entries
# 00:40:07 	7870000 variants (3281 variants per second), 7814868 VCF entries
# 00:40:10 	7880000 variants (3281 variants per second), 7824776 VCF entries
# 00:40:13 	7890000 variants (3281 variants per second), 7834682 VCF entries
# 00:40:16 	7900000 variants (3281 variants per second), 7844620 VCF entries
# 00:40:19 	7910000 variants (3281 variants per second), 7854519 VCF entries
# 00:40:23 	7920000 variants (3281 variants per second), 7864438 VCF entries
# 00:40:26 	7930000 variants (3280 variants per second), 7874361 VCF entries
# 00:40:29 	7940000 variants (3280 variants per second), 7884276 VCF entries
# 00:40:32 	7950000 variants (3280 variants per second), 7894218 VCF entries
#
# WARNINGS: Some warning were detected
# Warning type	Number of warnings
# WARNING_TRANSCRIPT_INCOMPLETE	11245
# WARNING_TRANSCRIPT_NO_START_CODON	1207143
# WARNING_TRANSCRIPT_NO_STOP_CODON	89303
#
#
# 00:40:33 Creating summary file: report.html
# 00:40:34 Creating genes file: ./report.genes.txt
# 00:40:35 done.
# 00:40:35 Done.
#
# (snpeff) ~ ❯❯❯ snpEff -v -stats report.html M_zebra_UMD2a.105 /Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset1_missing0.9_maf0.05_Victoria.recode.vcf.gz > /Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset1_missing0.9_maf0.05_Victoria_snpeff.recode.vcf
# 00:00:00 SnpEff version SnpEff 5.2 (build 2023-09-29 06:17), by Pablo Cingolani
# 00:00:00 Command: 'ann'
# 00:00:00 Reading configuration file 'snpEff.config'. Genome: 'M_zebra_UMD2a.105'
# 00:00:00 Reading config file: /Users/machiinagatoshi/snpEff.config
# 00:00:00 Reading config file: /Users/machiinagatoshi/mambaforge/envs/snpeff/share/snpeff-5.2-1/snpEff.config
# 00:00:02 done
# 00:00:02 Reading database for genome version 'M_zebra_UMD2a.105' from file '/Users/machiinagatoshi/mambaforge/envs/snpeff/share/snpeff-5.2-1/./data/M_zebra_UMD2a.105/snpEffectPredictor.bin' (this might take a while)
# 00:00:05 done
# 00:00:05 Loading Motifs and PWMs
# 00:00:05 Building interval forest
# 00:00:07 done.
# 00:00:07 Genome stats :
# #-----------------------------------------------
# # Genome name                : 'Maylandia_zebra'
# # Genome version             : 'M_zebra_UMD2a.105'
# # Genome ID                  : 'M_zebra_UMD2a.105[0]'
# # Has protein coding info    : true
# # Has Tr. Support Level info : true
# # Genes                      : 28622
# # Protein coding genes       : 27187
# #-----------------------------------------------
# # Transcripts                : 39681
# # Avg. transcripts per gene  : 1.39
# # TSL transcripts            : 0
# #-----------------------------------------------
# # Checked transcripts        : 
# #               AA sequences :  38233 ( 99.97% )
# #              DNA sequences :  38295 ( 96.51% )
# #-----------------------------------------------
# # Protein coding transcripts : 38246
# #              Length errors :     58 ( 0.15% )
# #  STOP codons in CDS errors :     12 ( 0.03% )
# #         START codon errors :   5769 ( 15.08% )
# #        STOP codon warnings :    547 ( 1.43% )
# #              UTR sequences :  19303 ( 48.65% )
# #               Total Errors :   5810 ( 15.19% )
# #-----------------------------------------------
# # Cds                        : 375491
# # Exons                      : 384819
# # Exons with sequence        : 384819
# # Exons without sequence     : 0
# # Avg. exons per transcript  : 9.70
# # WARNING!                   : Mitochondrion chromosome 'MT' does not have a mitochondrion codon table (codon table = 'Standard'). You should update the config file.
# #-----------------------------------------------
# # Number of chromosomes      : 1690
# # Chromosomes                : Format 'chromo_name size codon_table'
# #		'LG7'	64916660	Standard
# #		'LG23'	42088218	Standard
# #		'LG6'	39774503	Standard
# #		'LG1'	38676823	Standard
# #		'LG14'	37870038	Standard
# #		'LG3'	37314939	Standard
# #		'LG5'	36170306	Standard
# #		'LG17'	35776103	Standard
# #		'LG16'	34742138	Standard
# #		'LG22'	34725690	Standard
# #		'LG15'	34548429	Standard
# #		'LG12'	34086040	Standard
# #		'LG2'	32660920	Standard
# #		'LG11'	32446080	Standard
# #		'LG10'	32356376	Standard
# #		'LG13'	32072427	Standard
# #		'LG4'	30518969	Standard
# #		'LG20'	29789237	Standard
# #		'LG18'	29505383	Standard
# #		'LG19'	25963618	Standard
# #		'LG8'	23971387	Standard
# #		'LG9'	21023328	Standard
# #		'AGTA05000023.1'	3684755	Standard
# #		'AGTA05000024.1'	1539574	Standard
# #		'AGTA05000025.1'	1394719	Standard
# #		'AGTA05000026.1'	1368916	Standard
# # ...
# #		'AGTA05001688.1'	17223	Standard
# #		'MT'	16582	Standard
# #		'AGTA05001689.1'	7645	Standard
# #-----------------------------------------------
#
# 00:00:08 Predicting variants
# 00:00:11 	10000 variants (3974 variants per second), 9974 VCF entries
# 00:00:12 	20000 variants (4904 variants per second), 19948 VCF entries
# 00:00:17 	30000 variants (3567 variants per second), 29931 VCF entries
# 00:00:18 	40000 variants (3999 variants per second), 39915 VCF entries
# 00:00:20 	50000 variants (4327 variants per second), 49889 VCF entries
# 00:00:21 	60000 variants (4548 variants per second), 59869 VCF entries
# 00:00:23 	70000 variants (4703 variants per second), 69848 VCF entries
# 00:00:24 	80000 variants (4922 variants per second), 79826 VCF entries
# 00:00:26 	90000 variants (5028 variants per second), 89806 VCF entries
# 00:00:28 	100000 variants (5173 variants per second), 99775 VCF entries
# 00:00:29 	110000 variants (5279 variants per second), 109759 VCF entries
# 00:00:31 	120000 variants (5337 variants per second), 119725 VCF entries
# 00:00:32 	130000 variants (5434 variants per second), 129691 VCF entries
# 00:00:34 	140000 variants (5467 variants per second), 139669 VCF entries
# 00:00:35 	150000 variants (5575 variants per second), 149661 VCF entries
# 00:00:37 	160000 variants (5640 variants per second), 159637 VCF entries
# 00:00:38 	170000 variants (5693 variants per second), 169619 VCF entries
# 00:00:40 	180000 variants (5739 variants per second), 179566 VCF entries
# 00:00:41 	190000 variants (5732 variants per second), 189518 VCF entries
# 00:00:43 	200000 variants (5760 variants per second), 199496 VCF entries
# 00:00:45 	210000 variants (5779 variants per second), 209471 VCF entries
# 00:00:46 	220000 variants (5769 variants per second), 219418 VCF entries
# 00:00:48 	230000 variants (5811 variants per second), 229329 VCF entries
# 00:00:49 	240000 variants (5842 variants per second), 239255 VCF entries
# 00:00:51 	250000 variants (5866 variants per second), 249116 VCF entries
# 00:00:52 	260000 variants (5899 variants per second), 259065 VCF entries
# 00:00:54 	270000 variants (5923 variants per second), 268980 VCF entries
# 00:00:55 	280000 variants (5946 variants per second), 278894 VCF entries
# 00:00:57 	290000 variants (5967 variants per second), 288786 VCF entries
# 00:00:58 	300000 variants (5965 variants per second), 298680 VCF entries
# 00:01:00 	310000 variants (5994 variants per second), 308645 VCF entries
# 00:01:01 	320000 variants (6017 variants per second), 318613 VCF entries
# 00:01:03 	330000 variants (6018 variants per second), 328511 VCF entries
# 00:01:05 	340000 variants (6026 variants per second), 338468 VCF entries
# 00:01:06 	350000 variants (6040 variants per second), 348396 VCF entries
# 00:01:08 	360000 variants (6060 variants per second), 358360 VCF entries
# 00:01:09 	370000 variants (6061 variants per second), 368329 VCF entries
# 00:01:11 	380000 variants (6063 variants per second), 378298 VCF entries
# 00:01:12 	390000 variants (6069 variants per second), 388281 VCF entries
# 00:01:14 	400000 variants (6061 variants per second), 398256 VCF entries
# 00:01:16 	410000 variants (6053 variants per second), 408185 VCF entries
# 00:01:17 	420000 variants (6062 variants per second), 418162 VCF entries
# 00:01:19 	430000 variants (6056 variants per second), 428122 VCF entries
# 00:01:21 	440000 variants (6060 variants per second), 438094 VCF entries
# 00:01:22 	450000 variants (6067 variants per second), 448067 VCF entries
# 00:01:24 	460000 variants (6085 variants per second), 458048 VCF entries
# 00:01:25 	470000 variants (6084 variants per second), 468020 VCF entries
# 00:01:27 	480000 variants (6075 variants per second), 478001 VCF entries
# 00:01:29 	490000 variants (6076 variants per second), 487985 VCF entries
# 00:01:31 	500000 variants (6048 variants per second), 497959 VCF entries
# 00:01:33 	510000 variants (6043 variants per second), 507926 VCF entries
# 00:01:34 	520000 variants (6035 variants per second), 517887 VCF entries
# 00:01:36 	530000 variants (6031 variants per second), 527864 VCF entries
# 00:01:38 	540000 variants (6039 variants per second), 537840 VCF entries
# 00:01:39 	550000 variants (6050 variants per second), 547816 VCF entries
# 00:01:41 	560000 variants (6057 variants per second), 557796 VCF entries
# 00:01:42 	570000 variants (6061 variants per second), 567747 VCF entries
# 00:01:44 	580000 variants (6066 variants per second), 577708 VCF entries
# 00:01:45 	590000 variants (6069 variants per second), 587686 VCF entries
# 00:01:47 	600000 variants (6077 variants per second), 597660 VCF entries
# 00:01:48 	610000 variants (6082 variants per second), 607635 VCF entries
# 00:01:50 	620000 variants (6082 variants per second), 617617 VCF entries
# 00:01:52 	630000 variants (6081 variants per second), 627592 VCF entries
# 00:01:53 	640000 variants (6079 variants per second), 637565 VCF entries
# 00:01:55 	650000 variants (6084 variants per second), 647523 VCF entries
# 00:01:57 	660000 variants (6083 variants per second), 657485 VCF entries
# 00:01:59 	670000 variants (6073 variants per second), 667449 VCF entries
# 00:02:00 	680000 variants (6074 variants per second), 677430 VCF entries
# 00:02:02 	690000 variants (6079 variants per second), 687409 VCF entries
# 00:02:03 	700000 variants (6080 variants per second), 697366 VCF entries
# 00:02:05 	710000 variants (6081 variants per second), 707332 VCF entries
# 00:02:06 	720000 variants (6096 variants per second), 717303 VCF entries
# 00:02:08 	730000 variants (6104 variants per second), 727286 VCF entries
# 00:02:09 	740000 variants (6109 variants per second), 737261 VCF entries
# 00:02:11 	750000 variants (6107 variants per second), 747241 VCF entries
# 00:02:13 	760000 variants (6098 variants per second), 757208 VCF entries
# 00:02:14 	770000 variants (6097 variants per second), 767192 VCF entries
# 00:02:16 	780000 variants (6102 variants per second), 777175 VCF entries
# 00:02:18 	790000 variants (6099 variants per second), 787159 VCF entries
# 00:02:19 	800000 variants (6105 variants per second), 797130 VCF entries
# 00:02:21 	810000 variants (6112 variants per second), 807119 VCF entries
# 00:02:22 	820000 variants (6117 variants per second), 817095 VCF entries
# 00:02:24 	830000 variants (6115 variants per second), 827076 VCF entries
# 00:02:25 	840000 variants (6121 variants per second), 837054 VCF entries
# 00:02:27 	850000 variants (6118 variants per second), 847034 VCF entries
# 00:02:29 	860000 variants (6122 variants per second), 857014 VCF entries
# 00:02:30 	870000 variants (6121 variants per second), 866983 VCF entries
# 00:02:32 	880000 variants (6126 variants per second), 876954 VCF entries
# 00:02:34 	890000 variants (6123 variants per second), 886926 VCF entries
# 00:02:35 	900000 variants (6129 variants per second), 896906 VCF entries
# 00:02:37 	910000 variants (6131 variants per second), 906891 VCF entries
# 00:02:38 	920000 variants (6132 variants per second), 916872 VCF entries
# 00:02:40 	930000 variants (6126 variants per second), 926838 VCF entries
# 00:02:42 	940000 variants (6124 variants per second), 936799 VCF entries
# 00:02:43 	950000 variants (6117 variants per second), 946773 VCF entries
# 00:02:45 	960000 variants (6113 variants per second), 956734 VCF entries
# 00:02:47 	970000 variants (6117 variants per second), 966718 VCF entries
# 00:02:48 	980000 variants (6117 variants per second), 976698 VCF entries
# 00:02:50 	990000 variants (6125 variants per second), 986674 VCF entries
# 00:02:51 	1000000 variants (6130 variants per second), 996657 VCF entries
# 00:02:53 	1010000 variants (6134 variants per second), 1006638 VCF entries
# 00:02:54 	1020000 variants (6136 variants per second), 1016609 VCF entries
# 00:02:56 	1030000 variants (6140 variants per second), 1026577 VCF entries
# 00:02:57 	1040000 variants (6144 variants per second), 1036553 VCF entries
# 00:02:59 	1050000 variants (6145 variants per second), 1046497 VCF entries
# 00:03:01 	1060000 variants (6145 variants per second), 1056483 VCF entries
# 00:03:02 	1070000 variants (6143 variants per second), 1066470 VCF entries
# 00:03:04 	1080000 variants (6140 variants per second), 1076440 VCF entries
# 00:03:06 	1090000 variants (6136 variants per second), 1086419 VCF entries
# 00:03:08 	1100000 variants (6126 variants per second), 1096396 VCF entries
# 00:03:09 	1110000 variants (6124 variants per second), 1106364 VCF entries
# 00:03:11 	1120000 variants (6123 variants per second), 1116343 VCF entries
# 00:03:13 	1130000 variants (6123 variants per second), 1126317 VCF entries
# 00:03:14 	1140000 variants (6126 variants per second), 1136289 VCF entries
# 00:03:16 	1150000 variants (6129 variants per second), 1146265 VCF entries
# 00:03:17 	1160000 variants (6128 variants per second), 1156218 VCF entries
# 00:03:19 	1170000 variants (6127 variants per second), 1166199 VCF entries
# 00:03:21 	1180000 variants (6133 variants per second), 1176176 VCF entries
# 00:03:22 	1190000 variants (6134 variants per second), 1186157 VCF entries
# 00:03:24 	1200000 variants (6133 variants per second), 1196135 VCF entries
# 00:03:26 	1210000 variants (6129 variants per second), 1206118 VCF entries
# 00:03:27 	1220000 variants (6122 variants per second), 1216049 VCF entries
# 00:03:29 	1230000 variants (6123 variants per second), 1226006 VCF entries
# 00:03:31 	1240000 variants (6124 variants per second), 1235961 VCF entries
# 00:03:32 	1250000 variants (6126 variants per second), 1245939 VCF entries
# 00:03:34 	1260000 variants (6128 variants per second), 1255915 VCF entries
# 00:03:35 	1270000 variants (6132 variants per second), 1265890 VCF entries
# 00:03:37 	1280000 variants (6130 variants per second), 1275869 VCF entries
# 00:03:39 	1290000 variants (6132 variants per second), 1285853 VCF entries
# 00:03:40 	1300000 variants (6129 variants per second), 1295823 VCF entries
# 00:03:42 	1310000 variants (6125 variants per second), 1305802 VCF entries
# 00:03:44 	1320000 variants (6122 variants per second), 1315776 VCF entries
# 00:03:45 	1330000 variants (6124 variants per second), 1325749 VCF entries
# 00:03:47 	1340000 variants (6121 variants per second), 1335717 VCF entries
# 00:03:49 	1350000 variants (6119 variants per second), 1345685 VCF entries
# 00:03:50 	1360000 variants (6118 variants per second), 1355657 VCF entries
# 00:03:52 	1370000 variants (6121 variants per second), 1365630 VCF entries
# 00:03:54 	1380000 variants (6122 variants per second), 1375585 VCF entries
# 00:03:55 	1390000 variants (6125 variants per second), 1385552 VCF entries
# 00:03:57 	1400000 variants (6124 variants per second), 1395524 VCF entries
# 00:03:58 	1410000 variants (6123 variants per second), 1405495 VCF entries
# 00:04:00 	1420000 variants (6126 variants per second), 1415466 VCF entries
# 00:04:02 	1430000 variants (6123 variants per second), 1425437 VCF entries
# 00:04:03 	1440000 variants (6122 variants per second), 1435414 VCF entries
# 00:04:05 	1450000 variants (6123 variants per second), 1445393 VCF entries
# 00:04:07 	1460000 variants (6118 variants per second), 1455362 VCF entries
# 00:04:08 	1470000 variants (6117 variants per second), 1465345 VCF entries
# 00:04:10 	1480000 variants (6116 variants per second), 1475318 VCF entries
# 00:04:12 	1490000 variants (6118 variants per second), 1485298 VCF entries
# 00:04:13 	1500000 variants (6118 variants per second), 1495269 VCF entries
# 00:04:15 	1510000 variants (6121 variants per second), 1505246 VCF entries
# 00:04:17 	1520000 variants (6121 variants per second), 1515228 VCF entries
# 00:04:18 	1530000 variants (6122 variants per second), 1525211 VCF entries
# 00:04:20 	1540000 variants (6123 variants per second), 1535194 VCF entries
# 00:04:21 	1550000 variants (6121 variants per second), 1545178 VCF entries
# 00:04:23 	1560000 variants (6121 variants per second), 1555148 VCF entries
# 00:04:25 	1570000 variants (6120 variants per second), 1565136 VCF entries
# 00:04:26 	1580000 variants (6123 variants per second), 1575112 VCF entries
# 00:04:28 	1590000 variants (6124 variants per second), 1585089 VCF entries
# 00:04:29 	1600000 variants (6123 variants per second), 1595060 VCF entries
# 00:04:31 	1610000 variants (6120 variants per second), 1605033 VCF entries
# 00:04:33 	1620000 variants (6117 variants per second), 1615008 VCF entries
# 00:04:35 	1630000 variants (6115 variants per second), 1624982 VCF entries
# 00:04:36 	1640000 variants (6115 variants per second), 1634962 VCF entries
# 00:04:38 	1650000 variants (6118 variants per second), 1644934 VCF entries
# 00:04:39 	1660000 variants (6122 variants per second), 1654893 VCF entries
# 00:04:41 	1670000 variants (6124 variants per second), 1664873 VCF entries
# 00:04:43 	1680000 variants (6121 variants per second), 1674803 VCF entries
# 00:04:44 	1690000 variants (6119 variants per second), 1684783 VCF entries
# 00:04:46 	1700000 variants (6117 variants per second), 1694764 VCF entries
# 00:04:48 	1710000 variants (6115 variants per second), 1704735 VCF entries
# 00:04:49 	1720000 variants (6114 variants per second), 1714712 VCF entries
# 00:04:51 	1730000 variants (6115 variants per second), 1724683 VCF entries
# 00:04:53 	1740000 variants (6111 variants per second), 1734656 VCF entries
# 00:04:55 	1750000 variants (6110 variants per second), 1744596 VCF entries
# 00:04:56 	1760000 variants (6107 variants per second), 1754575 VCF entries
# 00:04:58 	1770000 variants (6107 variants per second), 1764514 VCF entries
# 00:05:00 	1780000 variants (6107 variants per second), 1774474 VCF entries
# 00:05:01 	1790000 variants (6109 variants per second), 1784462 VCF entries
# 00:05:03 	1800000 variants (6108 variants per second), 1794443 VCF entries
# 00:05:05 	1810000 variants (6106 variants per second), 1804414 VCF entries
# 00:05:06 	1820000 variants (6106 variants per second), 1814387 VCF entries
# 00:05:08 	1830000 variants (6105 variants per second), 1824368 VCF entries
# 00:05:09 	1840000 variants (6108 variants per second), 1834351 VCF entries
# 00:05:11 	1850000 variants (6107 variants per second), 1844328 VCF entries
# 00:05:13 	1860000 variants (6108 variants per second), 1854313 VCF entries
# 00:05:15 	1870000 variants (6104 variants per second), 1864289 VCF entries
# 00:05:16 	1880000 variants (6100 variants per second), 1874261 VCF entries
# 00:05:18 	1890000 variants (6100 variants per second), 1884238 VCF entries
# 00:05:20 	1900000 variants (6100 variants per second), 1894218 VCF entries
# 00:05:21 	1910000 variants (6100 variants per second), 1904196 VCF entries
# 00:05:23 	1920000 variants (6101 variants per second), 1914171 VCF entries
# 00:05:25 	1930000 variants (6097 variants per second), 1924130 VCF entries
# 00:05:27 	1940000 variants (6090 variants per second), 1934101 VCF entries
# 00:05:29 	1950000 variants (6085 variants per second), 1944046 VCF entries
# 00:05:31 	1960000 variants (6080 variants per second), 1954008 VCF entries
# 00:05:32 	1970000 variants (6080 variants per second), 1963978 VCF entries
# 00:05:34 	1980000 variants (6079 variants per second), 1973958 VCF entries
# 00:05:36 	1990000 variants (6079 variants per second), 1983937 VCF entries
# 00:05:37 	2000000 variants (6074 variants per second), 1993908 VCF entries
# 00:05:39 	2010000 variants (6071 variants per second), 2003884 VCF entries
# 00:05:41 	2020000 variants (6071 variants per second), 2013864 VCF entries
# 00:05:42 	2030000 variants (6072 variants per second), 2023842 VCF entries
# 00:05:44 	2040000 variants (6067 variants per second), 2033799 VCF entries
# 00:05:46 	2050000 variants (6063 variants per second), 2043781 VCF entries
# 00:05:48 	2060000 variants (6062 variants per second), 2053752 VCF entries
# 00:05:50 	2070000 variants (6057 variants per second), 2063728 VCF entries
# 00:05:52 	2080000 variants (6053 variants per second), 2073700 VCF entries
# 00:05:53 	2090000 variants (6056 variants per second), 2083674 VCF entries
# 00:05:55 	2100000 variants (6054 variants per second), 2093646 VCF entries
# 00:05:57 	2110000 variants (6051 variants per second), 2103629 VCF entries
# 00:05:59 	2120000 variants (6049 variants per second), 2113606 VCF entries
# 00:06:01 	2130000 variants (6045 variants per second), 2123571 VCF entries
# 00:06:02 	2140000 variants (6040 variants per second), 2133541 VCF entries
# 00:06:04 	2150000 variants (6040 variants per second), 2143512 VCF entries
# 00:06:06 	2160000 variants (6041 variants per second), 2153489 VCF entries
# 00:06:08 	2170000 variants (6035 variants per second), 2163464 VCF entries
# 00:06:09 	2180000 variants (6034 variants per second), 2173446 VCF entries
# 00:06:11 	2190000 variants (6031 variants per second), 2183423 VCF entries
# 00:06:13 	2200000 variants (6028 variants per second), 2193385 VCF entries
# 00:06:15 	2210000 variants (6023 variants per second), 2203359 VCF entries
# 00:06:17 	2220000 variants (6020 variants per second), 2213322 VCF entries
# 00:06:19 	2230000 variants (6009 variants per second), 2223259 VCF entries
# 00:06:22 	2240000 variants (5995 variants per second), 2233046 VCF entries
# 00:06:23 	2250000 variants (5995 variants per second), 2243011 VCF entries
# 00:06:25 	2260000 variants (5995 variants per second), 2252981 VCF entries
# 00:06:27 	2270000 variants (5994 variants per second), 2262951 VCF entries
# 00:06:29 	2280000 variants (5993 variants per second), 2272927 VCF entries
# 00:06:30 	2290000 variants (5994 variants per second), 2282907 VCF entries
# 00:06:32 	2300000 variants (5993 variants per second), 2292883 VCF entries
# 00:06:34 	2310000 variants (5992 variants per second), 2302855 VCF entries
# 00:06:35 	2320000 variants (5994 variants per second), 2312832 VCF entries
# 00:06:37 	2330000 variants (5992 variants per second), 2322805 VCF entries
# 00:06:39 	2340000 variants (5990 variants per second), 2332783 VCF entries
# 00:06:41 	2350000 variants (5987 variants per second), 2342755 VCF entries
# 00:06:43 	2360000 variants (5984 variants per second), 2352733 VCF entries
# 00:06:45 	2370000 variants (5979 variants per second), 2362694 VCF entries
# 00:06:46 	2380000 variants (5982 variants per second), 2372671 VCF entries
# 00:06:48 	2390000 variants (5982 variants per second), 2382640 VCF entries
# 00:06:49 	2400000 variants (5983 variants per second), 2392592 VCF entries
# 00:06:51 	2410000 variants (5983 variants per second), 2402571 VCF entries
# 00:06:53 	2420000 variants (5984 variants per second), 2412558 VCF entries
# 00:06:54 	2430000 variants (5986 variants per second), 2422540 VCF entries
# 00:06:56 	2440000 variants (5977 variants per second), 2432481 VCF entries
# 00:06:58 	2450000 variants (5978 variants per second), 2442457 VCF entries
# 00:07:00 	2460000 variants (5978 variants per second), 2452428 VCF entries
# 00:07:01 	2470000 variants (5979 variants per second), 2462373 VCF entries
# 00:07:03 	2480000 variants (5977 variants per second), 2472351 VCF entries
# 00:07:05 	2490000 variants (5978 variants per second), 2482325 VCF entries
# 00:07:06 	2500000 variants (5980 variants per second), 2492282 VCF entries
# 00:07:08 	2510000 variants (5979 variants per second), 2502248 VCF entries
# 00:07:10 	2520000 variants (5975 variants per second), 2512205 VCF entries
# 00:07:12 	2530000 variants (5975 variants per second), 2522156 VCF entries
# 00:07:14 	2540000 variants (5970 variants per second), 2532100 VCF entries
# 00:07:15 	2550000 variants (5969 variants per second), 2542006 VCF entries
# 00:07:17 	2560000 variants (5968 variants per second), 2551973 VCF entries
# 00:07:19 	2570000 variants (5968 variants per second), 2561879 VCF entries
# 00:07:20 	2580000 variants (5968 variants per second), 2571858 VCF entries
# 00:07:22 	2590000 variants (5968 variants per second), 2581815 VCF entries
# 00:07:24 	2600000 variants (5965 variants per second), 2591776 VCF entries
# 00:07:26 	2610000 variants (5962 variants per second), 2601710 VCF entries
# 00:07:28 	2620000 variants (5963 variants per second), 2611630 VCF entries
# 00:07:29 	2630000 variants (5963 variants per second), 2621560 VCF entries
# 00:07:31 	2640000 variants (5962 variants per second), 2631512 VCF entries
# 00:07:33 	2650000 variants (5961 variants per second), 2641451 VCF entries
# 00:07:34 	2660000 variants (5961 variants per second), 2651373 VCF entries
# 00:07:36 	2670000 variants (5958 variants per second), 2661266 VCF entries
# 00:07:38 	2680000 variants (5956 variants per second), 2671223 VCF entries
# 00:07:40 	2690000 variants (5957 variants per second), 2681153 VCF entries
# 00:07:42 	2700000 variants (5955 variants per second), 2691089 VCF entries
# 00:07:43 	2710000 variants (5952 variants per second), 2701017 VCF entries
# 00:07:45 	2720000 variants (5950 variants per second), 2710937 VCF entries
# 00:07:47 	2730000 variants (5950 variants per second), 2720878 VCF entries
# 00:07:49 	2740000 variants (5949 variants per second), 2730815 VCF entries
# 00:07:51 	2750000 variants (5947 variants per second), 2740748 VCF entries
# 00:07:52 	2760000 variants (5946 variants per second), 2750684 VCF entries
# 00:07:54 	2770000 variants (5944 variants per second), 2760624 VCF entries
# 00:07:56 	2780000 variants (5944 variants per second), 2770547 VCF entries
# 00:07:57 	2790000 variants (5944 variants per second), 2780485 VCF entries
# 00:07:59 	2800000 variants (5944 variants per second), 2790435 VCF entries
#
# WARNINGS: Some warning were detected
# Warning type	Number of warnings
# WARNING_TRANSCRIPT_INCOMPLETE	4234
# WARNING_TRANSCRIPT_NO_START_CODON	486344
# WARNING_TRANSCRIPT_NO_STOP_CODON	33824
#
#
# 00:08:00 Creating summary file: report.html
# 00:08:01 Creating genes file: ./report.genes.txt
# 00:08:02 done.
# 00:08:02 Done.
# ```
#

# %% [markdown]
# ## 解釈
#
# 今回はアウトプットの/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/snpeff/report.genes.txtを用いて解析をする

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %%
# read file

path_report_genes = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/snpeff/Victoria/report.genes_handedit.txt'
df_report_genes = pd.read_csv(path_report_genes, sep='\t')

path_Gene_IDs_in_HDR = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/result/Gene_IDs_in_HDR.txt'
df_Gene_IDs_in_HDR = pd.read_csv(path_Gene_IDs_in_HDR, sep='\t')

df_report_genes['Missense/Silent ratio'] = df_report_genes['variants_effect_missense_variant']/df_report_genes['variants_effect_synonymous_variant']

df_report_genes = pd.merge(df_report_genes, df_Gene_IDs_in_HDR, how='left', left_on='GeneId', right_on='Gene_IDs_in_HDR')
df_report_genes['Gene_IDs_in_HDR'] = df_report_genes['Gene_IDs_in_HDR'].fillna('nonHDR')
df_report_genes['Gene_IDs_in_HDR_tag'] = df_report_genes['Gene_IDs_in_HDR'].apply(lambda HDR: 'HDR' if HDR!='nonHDR' else 'nonHDR')
df_report_genes = df_report_genes.sort_values(['Gene_IDs_in_HDR_tag', 'TranscriptId'], ascending=[True, True])
df_report_genes = df_report_genes.drop_duplicates('GeneId')

df_report_genes

# %%
fig, axes = plt.subplots(figsize=(13, 5), nrows=1, ncols=4, tight_layout=True)

snptypes = ['variants_effect_missense_variant', 'variants_effect_synonymous_variant', 'Missense/Silent ratio']
colors = ['#D4748D', '#CCCCCC']  # HDR, nonHDR

for i, snptype in enumerate(snptypes):
    sns.violinplot(data=df_report_genes, x='Gene_IDs_in_HDR_tag', y=snptype, 
                   palette=colors, ax=axes[i], cut=0)
    axes[i].set_title(f'{snptype}')


plt.rc('svg', fonttype='none')
plt.savefig('/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/snpeff/Victoria/MissenseSilentRatio.svg', bbox_inches='tight')
plt.show()

# %%

# %%

# %%
