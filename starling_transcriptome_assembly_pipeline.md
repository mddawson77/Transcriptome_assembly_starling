## Starling Transcriptome Assembly Pipeline

### Initial Quality Check with SolexaQA++ (v3.1.4)

```
SolexaQA++ analysis ../rawdata/MM_GB_1_128.1.fastq.gz ../rawdata/MM_GB_1_128.2.fastq.gz ../rawdata/MM_GB_1_75.1.fastq.gz ../rawdata/MM_GB_1_75.2.fastq.gz ../rawdata/MM_GB_1_X46.1.fastq.gz ../rawdata/MM_GB_1_X46.2.fastq.gz
```
### Plot Results using R (v0.99.491)


### Error Correct with Rcorrector (v1.01)
```
perl /home/molly/bin/Rcorrector/run_rcorrector.pl -k 31 -t 8 \
-1 /home/molly/starling/rawdata/MM_GB_1_128.1.fastq.gz,/home/molly/starling/rawdata/MM_GB_1_75.1.fastq.gz,/home/molly/starling/rawdata/MM_GB_1_X46.1.fastq.gz \
-2 /home/molly/starling/rawdata/MM_GB_1_128.2.fastq.gz,/home/molly/starling/rawdata/MM_GB_1_75.2.fastq.gz,/home/molly/starling/rawdata/MM_GB_1_X46.2.fastq.gz
```

### Transcriptome Assembly and Trimming with Trinity and Trimmomatic (v2.2.0)
```
Trinity --seqType fq --max_memory 40G --trimmomatic --CPU 10 --full_cleanup --output Rcorr_trinity_tucoKidney \
--left /home/molly/starling/rcorrector_v1.01/MM_GB_1_128.1.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_75.1.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_X46.1.cor.fq.gz \
--right /home/molly/starling/rcorrector_v1.01/MM_GB_1_128.2.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_75.2.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_X46.2.cor.fq.gz \
--quality_trimming_params "ILLUMINACLIP:/opt/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
```

### Seperately Compile Left and Right Corrected Reads
```
cat MM_GB_1_128.1.cor.fq.gz,MM_GB_1_75.1.cor.fq.gz,MM_GB_1_X46.1.cor.fq.gz > compiled.1.cor.fq
cat MM_GB_1_128.2.cor.fq.gz,MM_GB_1_75.2.cor.fq.gz,MM_GB_1_X46.2.cor.fq.gz > compiled.2.cor.fq
```

### Pull out 40M Reads for BinPacker (v3) with Seqtk (v1.0-r82-dirty) 
```
seqtk sample -s1025340 /home/molly/starling/rcorrector_v1.01/compiled.1.cor.fq 40000000 > starling_compiled_40M-1.1.fq |
seqtk sample -s1025340 /home/molly/starling/rcorrector_v1.01/compiled.2.cor.fq 40000000 > starling_compiled_40M-1.2.fq
```

### Transcriptome Assembly with BinPacker (v1.0)
```
BinPacker -d -q -s fq -p pair -m RF -k 25 -g 200 -o Rcorr_binpacker_starling \
-l /home/molly/starling/seqtk_v1.0-r82-dirty/starling_compiled_40M-1.1.fq \
-r /home/molly/starling/seqtk_v1.0-r82-dirty/starling_compiled_40M-1.2.fq 
```

### Merge Assemblies with TransFuse (v0.4.6)
```
transfuse -t 10 -i 0.98 -o transfuse_starling \
-l /home/molly/starling/rcorrector_v1.01/MM_GB_1_128.1.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_75.1.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_X46.1.cor.fq.gz \
-r /home/molly/starling/rcorrector_v1.01/MM_GB_1_128.2.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_75.2.cor.fq.gz,/home/molly/starling/rcorrector_v1.01/MM_GB_1_X46.2.cor.fq.gz \
-a /home/molly/starling/binpacker_v1.0/Rcorr_binpacker/BinPacker.fa,/home/molly/starling/trinity_v2.2.0/Rcorr_trinity_starling.Trinity.fasta
```

### Evaluate Transfuse Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \

```

### Evaluate Transrate Transfuse Optimized Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \

```

##Continue Pipeline with Merged Transfuse Assembly
### Estimate Gene Expression with Kallisto (v0.42.4) 
```
q
```

### Estimate Gene Expression with Salmon (v0.5.1)
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

## Continue Pipeline with Optimized Merged Transfuse Assembly
### Estimate Gene Expression with Kallisto (v0.42.4)
```
`
```

### Estimate Gene Expression with Salmon (v0.5.1)
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
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \

```

### Evaluate Optimized Assembly Comleteness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \

```

### Compare all Metrics and Scores to Determine Best Assembly

### Annotate with dammit (v0.3)
```
``
```
