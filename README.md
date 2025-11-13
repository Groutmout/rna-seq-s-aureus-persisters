# RNA-seq Reproducibility Pipeline (Nextflow + Docker)

Ce dépôt contient un pipeline **Nextflow DSL2** entièrement reproductible permettant de reproduire l'analyse RNA-seq décrite dans l'article :

> *Lechner et al., Nature Communications (2020)* — Analyse transcriptomique de **Staphylococcus aureus USA300**.

Le pipeline suit les étapes :
FASTQ → trimming → download du génome → indexation → alignement → comptage → analyse différentielle.

---

##  Structure du dépôt

```
rna-clean/
├── main.nf
├── nextflow.config
├── run.sh
├── README.md
├── samples.tsv
├── containers/
│   ├── bowtie/
│   │     └── Dockerfile
│   ├── cutadapt/
│   │     └── Dockerfile
│   ├── deseq2/
│   │     ├── Dockerfile
│   │     └── scripts/
│   │          └── run_deseq2.R
│   ├── sratoolkit/
│   │     └── Dockerfile
│   └── subread/
│         └── Dockerfile
└── data/ (créé automatiquement)
```

---

##  Conteneurs utilisés

Toutes les étapes utilisent des images Docker **construites localement** à partir des Dockerfiles du dossier `containers/`.

Versions exactes :

* **SRA Toolkit** : latest
* **Cutadapt** : 1.11
* **Bowtie (Bowtie1)** : 0.12.7
* **Samtools** : 0.1.19
* **Subread / featureCounts** : 1.4.6-p3
* **DESeq2** : 1.16 (R 3.4.1 via micromamba)

Build des images :

```
cd containers/bowtie && docker build -t alantrbt/bowtie:latest .
cd containers/cutadapt && docker build -t alantrbt/cutadapt:1.11 .
cd containers/sratoolkit && docker build -t alantrbt/sratoolkit:latest .
cd containers/subread && docker build -t alantrbt/subread:latest .
cd containers/deseq2 && docker build -t alantrbt/deseq2:latest .
```

---

##  Données analysées

Les identifiants SRA proviennent de l'étude originale. Ils sont listés dans :

```
SRA_ref/sra_ids.txt
```

Ce sont des **lectures single-end**.

---

##  Exécution du pipeline

Assurez-vous que Nextflow et Docker sont installés.

### 1. Lancer le pipeline

```
bash run.sh
```

Ou directement :

```
nextflow run main.nf -with-docker
```

### 2. Résultats produits

Tous les résultats seront écrits dans `results/` :

```
results/
├── raw_fastq/
├── trimmed/
├── index/
├── aligned/
├── counts/
└── deseq2/
```

Le fichier final principal est :

```
results/deseq2/deseq2_results.csv
```

---

##  Description du workflow

Le pipeline contient 7 étapes :

1. **Téléchargement FASTQ** (SRA Toolkit)
2. **Trimming** (Cutadapt 1.11)
3. **Téléchargement du génome & GFF3** (NCBI — CP000253.1)
4. **Indexation** (Bowtie1)
5. **Alignement single-end** (Bowtie1 + Samtools)
6. **Comptage** (featureCounts)
7. **Analyse différentielle** (DESeq2 1.16)

Le code du pipeline complet est dans `main.nf`.

---

##  Analyse différentielle (DESeq2)

Le script R utilisé est embarqué dans le conteneur :

```
containers/deseq2/scripts/run_deseq2.R
```

Il :

* fusionne les fichiers `counts_*`
* génère la matrice d'expression
* applique DESeq2
* écrit le fichier `deseq2_results.csv`

---

##  Reproductibilité

* Toutes les versions logicielles sont figées dans les Dockerfiles.
* Le pipeline est entièrement décrit par `main.nf`.
* Le workflow peut être relancé avec :

```
nextflow run main.nf -resume
```

---

##  Références

* Article original : Nature Communications (2020)

---

## Equipe

- Turbot Alan
- Beaumatin Yohan
- Faramus François
- Bouteville Tristan

