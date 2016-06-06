## Starling Transcriptome Assembly Pipeline

### Initial Quality Check with SolexaQA++ (v3.1.4)

```
SolexaQA++ analysis ../rawdata/MM_GB_1_128.1.fastq.gz ../rawdata/MM_GB_1_128.2.fastq.gz ../rawdata/MM_GB_1_75.1.fastq.gz ../rawdata/MM_GB_1_75.2.fastq.gz ../rawdata/MM_GB_1_X46.1.fastq.gz ../rawdata/MM_GB_1_X46.2.fastq.gz
```
### Plot Results using R (v)


### Error Correct with Rcorrector (v1.01)
```
perl /home/molly/bin/Rcorrector/run_rcorrector.pl -k 31 -t 8 \
-1 /home/molly/starling/rawdata/MM_GB_1_128.1.fastq.gz,/home/molly/starling/rawdata/MM_GB_1_75.1.fastq.gz,/home/molly/starling/rawdata/MM_GB_1_X46.1.fastq.gz \
-2 /home/molly/starling/rawdata/MM_GB_1_128.2.fastq.gz,/home/molly/starling/rawdata/MM_GB_1_75.2.fastq.gz,/home/molly/starling/rawdata/MM_GB_1_X46.2.fastq.gz
```

### Transcriptome Assembly and Trimming with Trinity Trimmomatic (v2.2.0)
```
Trinity --seqType fq --max_memory 40G --trimmomatic --CPU 30 --full_cleanup --output Rcorr_trinity_starling \
--left /home/ubuntu/starling/rcorrector_v1.01/MM_GB_1_128.1.cor.fq.gz,/home/ubuntu/starling/rcorrector_v1.01/MM_GB_1_75.1.cor.fq.gz,/home/ubuntu/starling/rcorrector_v1.01/MM_GB_1_X46.1.cor.fq.gz \
--right /home/ubuntu/starling/rcorrector_v1.01/MM_GB_1_128.2.cor.fq.gz,/home/ubuntu/starling/rcorrector_v1.01/MM_GB_1_75.2.cor.fq.gz,/home/ubuntu/starling/rcorrector_v1.01/MM_GB_1_X46.2.cor.fq.gz \
--quality_trimming_params "ILLUMINACLIP:/home/ubuntu/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
```

### Generate Optimized Assembly Including Quality Check with Transrate
```
transrate -o starling_1 -t 16 \
-a Rcorr_trinity.Trinity.fasta \
--left /home/ubuntu/starling/rawdata/MM_GB_1_128.1.fastq.gz,/home/ubuntu/starling/rawdata/MM_GB_1_75.1.fastq.gz,/home/ubuntu/starling/rawdata/MM_GB_1_X46.1.fastq.gz \
--right /home/ubuntu/starling/rawdata/MM_GB_1_128.2.fastq.gz,/home/ubuntu/starling/rawdata/MM_GB_1_75.2.fastq.gz,/home/ubuntu/starling/rawdata/MM_GB_1_X46.2.fastq.gz
```

### Evaluate Original Assembly Completeness with BUSCO (v)
```
python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 16 -l ~/BUSCO_v1.1b1/vertebrata \
-o starling_1.1 -in /home/ubuntu/starling/transrate_v/Rcorr_trinity.Trinity.fasta
```

### Evaluate Optimized Assembly Completeness with BUSCO (v)
```
python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 16 -l ~/BUSCO_v1.1b1/vertebrata \
-o starling_1.2 -in /home/ubuntu/starling/transrate_v/Rcorr_trinity.Trinity.fasta
```

##Continue Pipeline with Optimized Assembly
### Estimate Gene Expression with Kallisto (v0.42.4) 
```
q
```

### Estimate Gene Expression with Salmon (v0.3.0)
```
w
```

### Filter out Contigs Based on Gene Expression < TPM =1 
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

## Continue Pipeline with Optimized Assembly
### Estimate Gene Expression with Kallisto (v0.42.4)
```
`
```

### Estimate Gene Expression with Salmon (v0.3.0)
```
`
```

### Filter out Contigs Based on Gene Expression < TPM =1
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
