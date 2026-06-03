**Вариант:** 12 - Kedro + Illumina + BWA MEM

### 1. Данные

Для работы использован набор прочтений:

```text
SRR001666
```

Ссылка на данные:

```text
https://www.ebi.ac.uk/ena/browser/view/SRR001666
```

Использованы paired-end FASTQ-файлы:

```text
data/fastq/SRR001666_1.fastq.gz
data/fastq/SRR001666_2.fastq.gz
```

Скачивание:

```bash
wget -O data/fastq/SRR001666_1.fastq.gz \
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001666/SRR001666_1.fastq.gz

wget -O data/fastq/SRR001666_2.fastq.gz \
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001666/SRR001666_2.fastq.gz
```

FASTQ-файлы не добавлены в репозиторий из-за размера.

### 2. Референсный геном

Использован референсный геном *Escherichia coli*:

```text
GCF_000005845.2
```

Скачивание:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
```

Файл референса:

```text
data/reference/ecoli_reference.fna
```

Референс и индексные файлы не добавлены в репозиторий.

### 3. Установка инструментов

Использованные инструменты:

```text
FastQC
BWA
Samtools
FreeBayes
Kedro
Kedro-Viz
```

Установка:

```bash
conda create -n bio_hw3 -c conda-forge -c bioconda \
  python=3.11 fastqc bwa samtools freebayes pip -y

conda activate bio_hw3

pip install kedro kedro-viz
```

### 4. Bash-реализация алгоритма

Bash-скрипт:

```text
scripts/run_mapping_quality_pipeline.sh
```

Скрипт выполняет:

```text
FastQC
bwa index
bwa mem
samtools view
samtools flagstat
parse %mapped
OK / not OK
samtools sort
samtools index
freebayes
```

Результат проверки качества картирования:

```text
Mapped reads: 98.02%
OK
```

Результат `samtools flagstat`:

```text
results/flagstat/SRR001666_flagstat.txt
```

Результат разбора `% mapped`:

```text
results/flagstat/SRR001666_mapping_quality.txt
```

VCF-файл:

```text
results/variants/SRR001666.vcf
```

### 5. Скрипт разбора `samtools flagstat`

```text
scripts/parse_flagstat.sh
```


### 6. Kedro: Hello World pipeline

Код тестового пайплайна:

```text
pipeline/src/pipeline/pipelines/hello/
```

Запуск:

```bash
cd pipeline
kedro run --pipelines hello
```

Лог запуска:

```text
logs/kedro_hello.log
```

### 7. Kedro: пайплайн оценки качества картирования

Код основного Kedro-пайплайна:

```text
pipeline/src/pipeline/pipelines/mapping/
```

Запуск:

```bash
cd pipeline
kedro run --pipelines mapping
```

Kedro-пайплайн выполняет:

```text
run_fastqc
index_reference
run_bwa_mem
convert_sam_to_bam
run_flagstat
parse_mapping_quality
sort_bam
index_sorted_bam
run_freebayes
```

Результаты Kedro-пайплайна:

```text
results/flagstat/SRR001666_kedro_flagstat.txt
results/flagstat/SRR001666_kedro_mapping_quality.txt
results/variants/SRR001666.kedro.vcf
```

### 8. Визуализация пайплайна

Визуализация выполнена с помощью Kedro-Viz:

```bash
cd pipeline
kedro viz
```

Файл с визуализацией:

```text
results/visualization/kedro-pipeline.png
```

### 9. Отличие DAG от блок-схемы алгоритма

Блок-схема алгоритма показывает порядок действий и условие %mapped > 90%.


Если условие выполняется, выводится `OK`, выполняются `samtools sort` и `freebayes`.
Если условие не выполняется, выводится `not OK`.

DAG в Kedro показывает зависимости между узлами пайплайна. Каждый узел соответствует отдельному шагу обработки данных. Условная логика `OK / not OK` реализована внутри функций пайплайна, поэтому на DAG она отображается не так явно, как на блок-схеме алгоритма.

### 10. Структура проекта

```text
.
├── README.md
├── .gitignore
├── scripts/
│   ├── run_mapping_quality_pipeline.sh
│   └── parse_flagstat.sh
├── pipeline/
│   ├── README.md
│   ├── pyproject.toml
│   ├── requirements.txt
│   ├── conf/
│   ├── src/
│   │   └── pipeline/
│   │       ├── pipeline_registry.py
│   │       └── pipelines/
│   │           ├── hello/
│   │           └── mapping/
│   └── tests/
├── results/
│   ├── flagstat/
│   │   ├── SRR001666_flagstat.txt
│   │   ├── SRR001666_mapping_quality.txt
│   │   ├── SRR001666_kedro_flagstat.txt
│   │   └── SRR001666_kedro_mapping_quality.txt
│   ├── variants/
│   │   ├── SRR001666.vcf
│   │   └── SRR001666.kedro.vcf
│   └── visualization/
│       └── kedro-pipeline.png
└── logs/
    ├── fastqc.log
    ├── bwa_mem.log
    ├── kedro_fastqc.log
    ├── kedro_bwa_mem.log
    └── kedro_hello.log
```

### 11. Файлы, не добавленные в репозиторий

```text
data/fastq/
data/reference/
results/alignment/
```
