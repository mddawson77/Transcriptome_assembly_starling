## Ranitomeya imitator Transciptome Assembly Pipeline

### Initial Quality Check with SolexaQA++ (v3.1.4)
```
SolexaQA++ analysis /home/molly/imitator/rawdata/ /home/molly/imitator/rawdata/ 
```
### Error Correct with Rcorrector (v1.0.1)
```
perl /home/ubuntu/Rcorrector/run_rcorrector.pl -k 31 -t 30 \
-1 /home/ubuntu/imitator/rawdata/KS001F_2_S37_L008_R1_001.fastq.gz \
-2 /home/ubuntu/imitator/rawdata/KS001F_2_S37_L008_R2_001.fastq.gz
```
### Trim and Assemble with Trimmomatic and Trinity (v2.2.0)
```
Trinity --seqType fq --max_memory 40G --trimmomatic --CPU 30 --full_cleanup --output Rcorr_trinity_imitator \
--left /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R1_001.cor.fq.gz \
--right /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R2_001.cor.fq.gz \
--quality_trimming_params "ILLUMINACLIP:/opt/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
```
### Assemble with BinPacker (v1.0)
```
/home/ubuntu/BinPacker/BinPacker -d -q -s fq -p pair -m RF -k 25 -g 200 -o Rcorr_binpacker_imitator \
-l /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R1_001.cor.fq.gz \
-r /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R2_001.cor.fq.gz
```
### Generate Optimized Trinity Assembly Including Quality Check with Transrate (v1.0.4beta)
```
/home/ubuntu/transrate-1.0.4beta.tar.gz -o imitator_1 -t 16 \
-a /home/ubuntu/imitator/trinity_v2.2.0/Rcorr_trinity_imitator.Trinity.fasta \
--left /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R1_001.cor.fq.gz \
--right /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R2_001.cor.fq.gz
```
### Generate Optimized BinPacker Assembly Including Quality Check with Transrate (v1.0.4beta)
```
/home/ubuntu/transrate-1.0.4beta.tar.gz -o imitator_1 -t 16 \
-a /home/ubuntu/imitator/binpacker_v1.0/Rcorr_binpacker_imitator/BinPacker.fa \
--left /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R1_001.cor.fq.gz \
--right /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R2_001.cor.fq.gz
```
### Evaluate Original Trinity Assembly Completeness with BUSCO (v1.22)
```
python3 /home/ubuntu/BUSCO_v1.22/BUSCO_v1.22.py -m trans --cpu 16 -l /home/ubuntu/vertebrata \
-o imitator_1.1 -in /home/ubuntu/imitator/trinity_v2.2.0/Rcorr_trinity_imitator.Trinity.fasta
```
### Evaluate Optimized Trinity Assembly Completenes with BUSCO (v1.22)
```
python3 /home/ubuntu/BUSCO_v1.22/BUSCO_v1.22.py -m trans --cpu 16 -l /home/ubuntu/vertebrata \
-o imitator_1.2 -in /home/ubuntu/imitator/transrate_v1.0.4beta/imitator_1/Rcorr_trinity_imitator.Trinity/good.Rcorr_trinity_imitator.Trinity.fasta
```
### Evaluate Original BinPacker Assembly Completeness with BUSCO (v1.22)
```
python3 /home/ubuntu/BUSCO_v1.22/BUSCO_v1.22.py -m trans --cpu 16 -l /home/ubuntu/vertebrata \
-o imitator_1.1 -in /home/ubuntu/imitator/binpacker_v1.0/Rcorr_binpacker_imitator/BinPacker.fa
```
### Evaluate Optimized BinPacker Assembly Completeness with BUSCO (v1.22)
```
python3 /home/ubuntu/BUSCO_v1.22/BUSCO_v1.22.py -m trans --cpu 18 -l /home/ubuntu/vertebrata \
-o imitator_1.2 -in /home/ubuntu/imitator/transrate_v1.0.4beta/binpacker_imitator/imitator_1/BinPacker/good.BinPacker.fa
```
### Compare all Metrics and Scores to Determine Best Trinity Assembly and Best BinPacker Assembly

### Merge Best Trinity Assembly and Best BinPacker Assembly with Transfuse (v0.5.0)
```
transfuse -t 40 -i 0.98 -o transfuse_imitator \
-l /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R1_001.cor.fq.gz \
-r /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R2_001.cor.fq.gz \
-a /home/ubuntu/imitator/transrate_v1.0.4beta/binpacker_imitator/imitator_1/BinPacker/good.BinPacker.fa,/home/ubuntu/imitator/transrate_v1.0.4beta/trinity_imitator/imitator_1/Rcorr_trinity_imitator.Trinity/good.Rcorr_trinity_imitator.Trinity.fasta
```
### Evaluate Original Transfuse Assembly Completenes with BUSCO (v1.22)
```
python3 /home/ubuntu/BUSCO_v1.22/BUSCO_v1.22.py -m trans --cpu 16 -l /home/ubuntu/vertebrata \
-o imitator_1.1 -in /home/ubuntu/imitator/transfuse_v0.5.0/transfuse_imitator_cons.fa 
```
### Evaluate Optimized Transfuse Assembly Completeness with BUSCO (v1.22) (Remember: Transrate is subprogram of Transfuse)
```
python3 /home/ubuntu/BUSCO_v1.22/BUSCO_v1.22.py -m trans --cpu 16 -l /home/ubuntu/vertebrata \
-o imitator_1.2 -in /home/ubuntu/imitator/transfuse_v0.5.0/transrate_transfuse_imitator_cons/good.transfuse_imitator_cons.fa
```
### Continue Pipeline with Original Transfuse Assembly

### Estimate Gene Expression with Kallisto (v0.43.0)
```
kallisto index -i kallisto.idx /home/ubuntu/imitator/transfuse_v0.5.0/transfuse_imitator_cons.fa
kallisto quant -t 32 -i kallisto.idx -o kallisto_orig_transfuse_assembly /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R1_001.cor.fq.gz /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R2_001.cor.fq.gz
```
### Estimate Gene Expression with Salmon (v0.6.0)
```
/home/ubuntu/salmon-0.6.0/bin/salmon index -t /home/ubuntu/imitator/transfuse_v0.5.0/transfuse_imitator_cons.fa -i salmon.idx --type quasi -k 31
/home/ubuntu/salmon-0.6.0/bin/salmon quant -p 32 -i salmon.idx -l IU -1 /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R1_001.cor.fq -2 /home/ubuntu/imitator/rcorrector_v1.0.1/KS001F_2_S37_L008_R2_001.cor.fq -o salmon_orig_transfuse_assembly 
```
### Filter Out Contigs Based on Gene Expression <TPM = 1
```
awk '1>$5{next}1' /home/ubuntu/imitator/kallisto_v0.43.0/kallisto_orig_transfuse_assembly/abundance.tsv | awk '{print $1}' > kallist
awk '1>$4{next}1' /home/ubuntu/imitator/salmon_v0.6.0/salmon_orig_transfuse_assembly/quant.sf | sed  '1,10d' | awk '{print $1}' > salist
cat kallist salist | sort -u > uniq_list

python filter.py /home/ubuntu/imitator/transfuse_v0.5.0/transfuse_imitator_cons.fa uniq_list > Highexp.fasta
```
### 

