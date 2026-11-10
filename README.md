# simple‑ribo
Automated Ribo‑seq & Bulk RNA‑seq end‑to‑end analysis pipeline.

Workflow
- Ribo‑seq: QC, adapter‑trim, rRNA‑tRNA depletion, STAR alignment, riboWaltz P‑site calibration, periodicity, codon occupancy
- Bulk RNA‑seq: QC → adapter trimming → alignment → gene quantification → PCA / DEG / GO‑KEGG enrichment
- Combined mode: compute Translation Efficiency(TE) and differential TE analysis

> Tested environment: Linux, conda.
> Reference genome / index files are NOT included in this repository due to large size.

## ⚡ Quick install
```bash
git clone https://github.com/GY-Zhang2023/simple-ribo.git
cd simple-ribo
bash install.sh
mkdir -p ~/conda_pkgs
export CONDA_PKGS_DIRS=~/conda_pkgs
conda env create -f environment.yml
```
## Activate environmentbash
```bash
source env/ribo_LAN_env.sh
bash tools/check_env.sh
ribo_LAN -h
rna_LAN -h
```
## Custom resource limit (environment variable override)
```bash
export N_CORES=8
export MEMORY_LIMIT_GB=50
# run pipeline commands ...
```
## Supported species
- human (GRCh38.p14, GENCODE v49)
- mouse (GRCm39, GENCODE vM35)
- other species in future development

## Reference preparation
Read document: [`docs/reference_prepare.md`](docs/reference_prepare.md)

> 
> You must download genome fasta, GTF annotation and build STAR / bowtie2 index manually, fill absolute paths in pipeline script `SPECIES_CONFIG` dictionary.

## Output directory structure
Read document: docs/output_structure.md
## Input template
Sample csv templates are under folder `example/`
- `ribo_sample.csv` for ribo‑seq
- `rna_sample.csv` for bulk RNA‑seq
- `group.csv` for differential analysis
## Usage example
```bash
mkdir -p your_ribo_project
cd your_ribo_project
mkdir -p 1.raw_fastq
cp YOUR_DATA_SAMPLE*.fastq.gz 1.raw_fastq/
cp ~/simple-ribo/example/ribo_sample.csv 1.raw_fastq/    # sample_name,adapter_3,adapter_5
cp ~/simple-ribo/example/group.csv 1.raw_fastq/          # sample_name,group_name
```
Ribo‑seq
```bash
ribo_LAN --step1 ribo_sample.csv
ribo_LAN --step2
ribo_LAN --step3
ribo_LAN --step4
ribo_LAN --step5
ribo_LAN --step6 group.csv gene_list.txt
```
Bulk RNA‑seq
```bash
rna_LAN --step1 sample.csv
rna_LAN --step2
rna_LAN --step3 group.csv
rna_LAN --step4
rna_LAN --step5 group.csv
```
## Software version list
| Software  | Version   |
|-----------|-----------|
| R         | 4.4.x     |
| cutadapt  | >=4.0     |
| STAR      | >=2.7.11  |
| HISAT2    | >=2.2.1   |
| bowtie2   | >=2.5     |
| samtools  | >=1.18    |
| riboWaltz | >=1.2.0   |

## Notes

- Do NOT run analysis inside source‑code git repository folder, use independent project workspace.
- All log / plot outputs are in project sub‑folders.
