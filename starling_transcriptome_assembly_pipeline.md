## Starling Transcriptome Assembly Pipeline

### Initial Quality Check with SolexaQA++ (v3.1.4)

```
SolexaQA++ analysis ../rawdata/MM_GB_1_128.1.fastq.gz ../rawdata/MM_GB_1_128.2.fastq.gz ../rawdata/MM_GB_1_75.1.fastq.gz ../rawdata/MM_GB_1_75.2.fastq.gz ../rawdata/MM_GB_1_X46.1.fastq.gz ../rawdata/MM_GB_1_X46.2.fastq.gz
```
### Plot Results using R 

### Error Correct with Rcorrector (v1.0.1)
```
perl /home/molly/bin/Rcorrector/run_rcorrector.pl -k 31 -t 8 \
-1 ../rawdata/MM_GB_1_128.1.fastq.gz,../rawdata/MM_GB_1_75.1.fastq.gz,../rawdata/MM_GB_1_X46.1.fastq.gz \
-2 ../rawdata/MM_GB_1_128.2.fastq.gz,../rawdata/MM_GB_1_75.2.fastq.gz,../rawdata/MM_GB_1_X46.2.fastq.gz
```

### Concatenate Left and Right Corrected Reads
```
cat MM_GB_1_128.1.cor.fq.gz MM_GB_1_75.1.cor.fq.gz MM_GB_1_X46.1.cor.fq.gz > compiled.1.cor.fq.gz
cat MM_GB_1_128.2.cor.fq.gz MM_GB_1_75.2.cor.fq.gz MM_GB_1_X46.2.cor.fq.gz > compiled.2.cor.fq.gz
```

### Unzip Concatenated Files
```
gunzip compiled.1.cor.fq.gz
gunzip compiled.2.cor.fq.gz
```

### Pullout 40M Reads at Random with Seqtk (v1.0-r82-dirty)
```
seqtk sample -s1025340 /home/molly/starling/rcorrector_v1.01/compiled.1.cor.fq 40000000 > starling_new_compiled_40M-1.1.fq |
seqtk sample -s1025340 /home/molly/starling/rcorrector_v1.01/compiled.2.cor.fq 40000000 > starling_new_compiled_40M-1.2.fq
```

### Trim and Assemble with Trimmomatic and Trinity (v2.2.0)
```
Trinity --seqType fq --max_memory 40G --trimmomatic --CPU 10 --full_cleanup --output Rcorr_trinity_starling_40M \
--left /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq \
--right /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.2.fq \
--quality_trimming_params "ILLUMINACLIP:/opt/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25
```

### Assemble with BinPacker (v1.0)
```
/opt/BinPacker/BinPacker -d -q -s fq -p pair -m RF -k 25 -g 200 -o Rcorr_binpacker_starling \
-l /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq \
-r /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.2.fq 
```

### Generate Optimized Trinity Assembly Including Quality Check with Transrate (v1.0.4beta)
```
/opt/transrate-1.0.4beta/transrate -o starling_1 -t 10 \
-a /home/molly/starling/trinity_v2.2.0/Rcorr_trinity_starling_40M.Trinity.fasta \
--left /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq \
--right /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.2.fq 
```

### Generate Optimized BinPacker Assembly Including Quality Check with Transrate (v1.0.4beta)
```
/opt/transrate-1.0.4beta/transrate -o starling_1 -t 10 \
-a /home/molly/starling/binpacker_v1.0/Rcorr_binpacker_starling/BinPacker.fa \
--left /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq \
--right /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.2.fq 
```

### Evaluate Original Trinity Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o starling_1.1 -in /home/molly/starling/trinity_v2.2.0/Rcorr_trinity_starling_40M.Trinity.fasta
```

### Evaluate Optimized Trinity Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o starling_1.2 -in /home/molly/starling/transrate_v1.0.1/40M_trinity_starling/starling_1/Rcorr_trinity_starling_40M.Trinity/good.Rcorr_trinity_starling_40M.Trinity.fasta
```
### Evaluate Original BinPacker Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o starling_1.1 -in /home/molly/starling/binpacker_v1.0/Rcorr_binpacker_starling/BinPacker.fa 
```

### Evaluate Optimized BinPacker Assembly Completeness with BUSCO (v1.1b1) 
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o starling_1.2 -in /home/molly/starling/transrate_v1.0.1/40M_binpacker_starling/starling_1/BinPacker/good.BinPacker.fa
```

### Compare all Metrics and Scores to Determine Best Trinity Assembly and Best BinPacker Assembly

### Merge Best 40M Trinity Assembly and Best 40M BinPacker Assembly with Transfuse (v0.5.0)
```
/opt/transfuse-0.5.0-linux-x86_64/transfuse -t 10 -i 0.98 -o transfuse_starling \
-l /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq \
-r /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.2.fq \
-a /home/molly/starling/binpacker_v1.0/Rcorr_binpacker_starling/BinPacker.fa,/home/molly/starling/trinity_v2.2.0/Rcorr_trinity_starling_40M.Trinity.fasta
```

### Evaluate Original Transfuse Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o starling_1.1 -in /home/molly/starling/transfuse_v0.5.0/transfuse_starling_cons.fa
```

### Evaluate Optimized Transfuse Assembly Completeness with BUSCO (v1.1b1)(Remember: Transrate is a subprogram of Transfuse)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o starling_1.2 -in /home/molly/starling/transfuse_v0.5.0/transrate_transfuse_starling_cons/good.transfuse_starling_cons.fa
```

##Continue Pipeline with Original Transfuse Assembly

### Estimate Gene Expression with Kallisto (v0.42.4) 
```
kallisto index -i kallisto.idx /home/molly/starling/transfuse_v0.5.0/transfuse_starling_cons.fa
kallisto quant -t 10 -i kallisto.idx -o kallisto_orig_transfuse_assembly /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.2.fq  
```

### Estimate Gene Expression with Salmon (v0.5.1)
```
/opt/salmon-0.5.1/bin/salmon index -t /home/molly/starling/transfuse_v0.5.0/transfuse_starling_cons.fa -i salmon.idx --type quasi -k 31
/opt/salmon-0.5.1/bin/salmon quant -p 32 -i salmon.idx -l IU -1 /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq -2 /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq -o salmon_orig_transfuse_assembly
```

### Filter out Contigs Based on Gene Expression < TPM =1 from Kallisto and Salmon
```
awk '1>$5{next}1' /home/molly/starling/kallisto_v0.42.4/kallisto_orig_transfuse_assembly/abundance.tsv | awk '{print $1}' > kallist
awk '1>$4{next}1' /home/molly/starling/salmon_v0.5.1/salmon_orig_transfuse_assembly/quant.sf | sed  '1,10d' | awk '{print $1}' > salist
cat kallist salist | sort -u > uniq_list

python filter.py /home/molly/starling/transfuse_v0.5.0/transfuse_starling_cons.fa uniq_list > Highexp.fasta
```

### Evaluate TPM Filtration Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o starling_2.1 -in /home/molly/starling/TPM/TPM_orig_transfuse_assembly/Highexp.fasta
```

### Generate Optimized TPM Filtration Assembly Including Quality Check with Transrate (v1.0.4beta)
```
/opt/transrate-1.0.4beta/transrate -o starling_2 -t 10 \
-a /home/molly/starling/TPM/TPM_orig_transfuse_assembly/Highexp.fasta \
--left /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq \
--right /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.2.fq 
```

### Evaluate Optimized TPM Filtration Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o starling_2.2 -in /home/molly/starling/transrate_v1.0.1/transfuse_starling/starling_2/Highexp/good.Highexp.fasta
```

## Continue Pipeline with Optimized Transfuse Assembly

### Estimate Gene Expression with Kallisto (v0.42.4)
```
kallisto index -i kallisto.idx /home/molly/starling/transfuse_v0.5.0/transrate_transfuse_starling_cons/good.transfuse_starling_cons.fa
kallisto quant -t 10 -i kallisto.idx -o kallisto_opt_transfuse_assembly /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.2.fq  
```

### Estimate Gene Expression with Salmon (v0.5.1)
```
/opt/salmon-0.5.1/bin/salmon index -t /home/molly/starling/transfuse_v0.5.0/transrate_transfuse_starling_cons/good.transfuse_starling_cons.fa -i salmon.idx --type quasi -k 31
/opt/salmon-0.5.1/bin/salmon quant -p 32 -i salmon.idx -l IU -1 /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq -2 /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq -o salmon_opt_transfuse_assembly 
```

### Filter out Contigs Based on Gene Expression < TPM =1 from Kallisto and Salmon
```
awk '1>$5{next}1' /home/molly/starling/kallisto_v0.42.4/kallisto_opt_transfuse_assembly/abundance.tsv | awk '{print $1}' > kallist
awk '1>$4{next}1' /home/molly/starling/salmon_v0.5.1/salmon_opt_transfuse_assembly/quant.sf | sed  '1,10d' | awk '{print $1}' > salist
cat kallist salist | sort -u > uniq_list

python filter.py /home/molly/starling/transfuse_v0.5.0/transrate_transfuse_starling_cons/good.transfuse_starling_cons.fa uniq_list > Highexp.fasta
```

### Evaluate TPM Filtration Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o starling_3.1 -in /home/molly/starling/TPM/TPM_opt_transfuse_assembly/Highexp.fasta
```

### Generate Optimized TPM Filtration Assembly Including Quality Check with Transrate (v1.0.4beta)
```
/opt/transrate-1.0.4beta/transrate -o starling_3 -t 10 \
-a /home/molly/starling/TPM/TPM_opt_transfuse_assembly/Highexp.fasta \
--left /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.1.fq \
--right /home/molly/starling/seqtk_v1.0-r82-dirty/starling_new_compiled_40M-1.2.fq 
```

### Evaluate Optimized TPM Filtration Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o starling_3.2 -in /home/molly/starling/transrate_v1.0.1/transfuse_starling/starling_3/Highexp/good.Highexp.fasta
```

## Compare all Metrics and Scores to Determine Best Assembly

### Annotate with dammit (v0.3)
```
/usr/local/bin/dammit annotate /home/molly/starling/transfuse_v0.5.0/transrate_transfuse_starling_cons/good.transfuse_starling_cons.fa --busco-group vertebrata --n_threads 10 --database-dir /data/dammit_databases --full
```
