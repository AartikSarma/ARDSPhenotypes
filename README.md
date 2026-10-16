# ARDS Phenotypes

Analysis code for the transcriptomic study of hyperinflammatory (HyperARDS) and hypoinflammatory (HypoARDS) phenotypes of acute respiratory distress syndrome in tracheal aspirates, published in the *American Journal of Respiratory and Critical Care Medicine* (2025;211(12):2352).

`scripts/ajrccm_2025/` reproduces the paper's transcriptomic analyses from the data deposited in GEO. `scripts/preprint_2022/` holds the scripts behind the 2022 medRxiv preprint, kept for reference; they read data that are not public.

## Data

All inputs are public and are downloaded by script 0 from GEO SuperSeries [GSE200849](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200849):

| Series | Data |
|---|---|
| GSE200847 | Non-COVID ARDS, bulk RNA-seq, 41 subjects |
| GSE350140 | COVID ARDS, bulk RNA-seq, 72 samples from 39 subjects on repeated ventilator days |
| GSE200848 | ARDS, single-cell RNA-seq, 9 subjects (one Seurat object) |

The repository also holds one published table, `published_results/Supplementary Data E1 - Bulk RNASeq.xlsx`. The published bulk models used exact ages and voom weights from samples that are not in GEO, so the gene set analyses start from this table rather than from refits.

GEO holds no raw reads (UCSF IRB restriction). Age, a covariate of the bulk models, comes from GEO, where patients older than 89 are recorded as ">89"; the analyses set those ages to 90. Additional data are available on reasonable request.

## Setup

The code runs in R 4.5.2 (pinned in `.r-version`). R packages are managed with [uvr](https://github.com/nbafrank/uvr) and locked in `uvr.lock`. Install them with:

```sh
uvr sync
```

### CellChat on Apple silicon

Script 6 uses CellChat 1.5.0, which `uvr sync` installs from its GitHub tag. The tag ships compiled x86_64 objects, so on Apple silicon the installed package fails to load ("incompatible architecture"). Rebuild it from source after `uvr sync`, from the repository root:

1. Install the Xcode command line tools, if they are not already installed:

   ```sh
   xcode-select --install
   ```

2. Download the tagged source outside the repository:

   ```sh
   git clone --depth 1 --branch v1.5.0 https://github.com/sqjin/CellChat /tmp/CellChat-1.5.0
   ```

3. Delete the prebuilt objects:

   ```sh
   rm /tmp/CellChat-1.5.0/src/*.o /tmp/CellChat-1.5.0/src/*.so
   ```

4. Build and install the package into the project library with R 4.5.2:

   ```sh
   R CMD INSTALL --preclean -l .uvr/library /tmp/CellChat-1.5.0
   ```

5. Check that it loads:

   ```sh
   Rscript -e 'library(CellChat); packageVersion("CellChat")'
   ```

   The project's `.Rprofile` puts `.uvr/library` first, so this should print `1.5.0`. A later `uvr sync` may reinstall the prebuilt package; if so, repeat steps 3 to 5.

## Running the analyses

Run all the scripts in order from the repository root:

```sh
uvr run run_all.R
```

`run_all.R` runs each script in a fresh R session and stops at the first one that fails. To resume from a given script, pass its number, for example `uvr run run_all.R -- 4`. To run one script on its own:

```sh
uvr run "scripts/ajrccm_2025/ARDS Phenotypes - 0 - Import and clean GEO data.R"
```

| Script | Analysis | Writes |
|---|---|---|
| 0 | Download and clean the GEO data; parse age | `data/geo/clean/` |
| 1 | Non-COVID bulk RNA-seq: limma-voom | bulk differential expression |
| 2 | COVID bulk RNA-seq: dream mixed model | bulk differential expression |
| 3 | Single-cell differential expression by cell type: MAST | Supplementary Data E5 |
| 4 | Reproducible genes; Reactome gene set enrichment | Supplementary Data E2, E4 |
| 5 | HyperARDS signature scores: GSVA | Supplementary Data E3 |
| 6 | Ligand-receptor analysis: CellChat | Supplementary Data E6 |
| 7 | Figures | `data/results/ajrccm_2025/figures/` |

Results go to `data/results/ajrccm_2025/`. Script 3 takes over an hour. Script 7 needs scripts 0, 3, 4, 5 and 6. Each script's header states its model, its inputs and its outputs.

Fold changes and enrichment scores are HyperARDS relative to HypoARDS throughout.

## License

See `LICENSE`.
