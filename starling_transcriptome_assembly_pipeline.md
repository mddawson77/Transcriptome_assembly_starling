## Starling Transcriptome Assembly Pipeline

### Initial Quality Check with Solexa QA++ (v3.1.4)

```
SolexaQA++ analysis ../rawdata/MM_GB_1_128.1.fastq.gz ../rawdata/MM_GB_1_128.2.fastq.gz ../rawdata/MM_GB_1_75.1.fastq.gz ../rawdata/MM_GB_1_75.2.fastq.gz ../rawdata/MM_GB_1_X46.1.fastq.gz ../rawdata/MM_GB_1_X46.2.fastq.gz
```

### Error Correct with Rcorrector (v1.01)
```
perl /home/molly/bin/Rcorrector/run_rcorrector.pl -k 31 -t 8 \
-1 ../rawdata/MM_GB_1_128.1.fastq.gz,../rawdata/MM_GB_1_75.1.fastq.gz,../rawdata/MM_GB_1_X46.1.fastq.gz \
-2 ../rawdata/MM_GB_1_128.2.fastq.gz,../rawdata/MM_GB_1_75.2.fastq.gz,../rawdata/MM_GB_1_X46.2.fastq.gz
```

### Transcriptome Assembly and Trimming with Trinity Trimmomatic (v2.2.0)
```
Trinity --seqType fq --max_memory 40G --trimmomatic --CPU 10 --full_cleanup --output Rcorr_trinity_starling \
--left /home/molly/starling/rcorrector_v1.01/MM_GB_1_128.1.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_75.1.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_X46.1.cor.fq.gz \
--right /home/molly/starling/rcorrector_v1.01/MM_GB_1_128.2.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_75.2.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_X46.2.cor.fq.gz \
--quality_trimming_params "ILLUMINACLIP:/opt/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
```

### Generate Optimized Assembly Including Quality Check with Transrate
```
s
```

### Evaluate Original Assembly Completeness with BUSCO 
```
s
```

##Continue Pipeline with Optimized Assembly
```
s
```

### Estimate Gene Expression with Kallisto (v0.42.4) 
```
q
```

### Estimate Gene Expression with Salmon (v0.3.0)
```
w
```

### Filter out Contigs Based on Gene Expression < TPM=1 
```
q
```

### Generate 2nd Optimized Assembly Including Quality Check with Transrate (v1.0.1)
```
`
```

### Evaluate TPM Filtration Assembly Completeness with BUSCO (v1.1b1)
```
`
```

### Evaluate 2nd Optimized Assembly Completeness with BUSCO (v1.1b1)
```
`
```

## Continue Pipeline with Original Assembly
### Estimate Gene Expression with Kallisto (v0.42.4)
```
`
```

### Estimate Gene Expression with Salmon (v0.3.0)
```
`
```

### Filter out Contigs Based on Gene Expression < TPM=1
```
`
```

### Generate Optimized Assembly Including Quality Check with Transrate (v1.0.1)
```
`
```

### Evaluate TPM Filtration Assembly Completeness with BUSCO (v1.1b1)
```
`
```

### Evaluate Optimized Assembly Comleteness with BUSCO (v1.1b1)
```
`
```

### Compare all Metrics and Scores to Determine Best Assembly

### Annotate with dammit (v0.3)
```
``
```
