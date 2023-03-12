# spipe.v1.0.0 is already installed as conda env
conda activate spipe.v1.0.0

PBS='/data2/ivanir/Feline2023/ParseBS'

PBS2='/data2/hanna/synaptogenesis/newvolume'

PATH="/data2/ivanir/Feline2023/ParseBS/ParseBiosciences-Pipeline.1.0.3p:$PATH"

cd $PBS/newvolume/genomes/
cat Homo_sapiens.GRCh38.108.gtf EmGFP.gtf > Homo_sapiens.GRCh38.108.EmGFP.gtf
cat Homo_sapiens.GRCh38.dna.primary_assembly.fa EmGFP.fa > Homo_sapiens.GRCh38.dna.primary_assembly.EmGFP.fa

cat Mus_musculus.GRCm39.108.gtf EmGFP.gtf > Mus_musculus.GRCm38.108.EmGFP.gtf
cat Mus_musculus.GRCm39.dna.primary_assembly.fa EmGFP.fa > Mus_musculus.GRCm38.dna.primary_assembly.EmGFP.fa

split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.EmGFP.fa \
--genes $PBS/newvolume/genomes/Homo_sapiens.GRCh38.108.EmGFP.gtf \
--output_dir $PBS/newvolume/genomes/hg38

split-pipe \
--mode mkref \
--genome_name mm10 \
--fasta $PBS/newvolume/genomes/Mus_musculus.GRCm38.dna.primary_assembly.EmGFP.fa \
--genes $PBS/newvolume/genomes/Mus_musculus.GRCm38.108.EmGFP.gtf \
--output_dir $PBS/newvolume/genomes/mm10

cat Homo_sapiens.GRCh38.108.gtf hEmGFP.gtf > Homo_sapiens.GRCh38.108.hEmGFP.gtf
cat Homo_sapiens.GRCh38.dna.primary_assembly.fa hEmGFP.fa > Homo_sapiens.GRCh38.dna.primary_assembly.hEmGFP.fa

cat Mus_musculus.GRCm39.108.gtf mEmGFP.gtf > Mus_musculus.GRCm38.108.mEmGFP.gtf
cat Mus_musculus.GRCm39.dna.primary_assembly.fa mEmGFP.fa > Mus_musculus.GRCm38.dna.primary_assembly.mEmGFP.fa

# Genome indexing
split-pipe \
--mode mkref \
--genome_name hg38 mm10 \
--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.EmGFP.fa $PBS/newvolume/genomes/Mus_musculus.GRCm38.dna.primary_assembly.EmGFP.fa \
--genes $PBS/newvolume/genomes/Homo_sapiens.GRCh38.108.EmGFP.gtf $PBS/newvolume/genomes/Mus_musculus.GRCm38.108.EmGFP.gtf \
--output_dir $PBS/newvolume/genomes/hg38_mm10

cd $PBS/newvolume/expdata/

cat SLX-22602.DNAA007.HGMLNDMXY.s_1.r_1.fq.gz SLX-22602.DNAA007.HGMLNDMXY.s_2.r_1.fq.gz > SLX-22602.tmp.r_1.fq.gz
cat SLX-22602.DNAA007.HGMLNDMXY.s_1.r_2.fq.gz SLX-22602.DNAA007.HGMLNDMXY.s_2.r_2.fq.gz > SLX-22602.tmp.r_2.fq.gz

cat SLX-22602.HGMLNDMXY.s_2.r_1.lostreads.fq.gz SLX-22602.HGMLNDMXY.s_2.r_1.lostreads.fq.gz > SLX-22602.lostreads.tmp.r_1.fq.gz
cat SLX-22602.HGMLNDMXY.s_2.r_1.lostreads.fq.gz SLX-22602.HGMLNDMXY.s_2.r_2.lostreads.fq.gz > SLX-22602.lostreads.tmp.r_2.fq.gz

cat SLX-22602.tmp.r_1.fq.gz SLX-22602.lostreads.tmp.r_1.fq.gz > SLX-22602.r_1.fq.gz
cat SLX-22602.tmp.r_2.fq.gz SLX-22602.lostreads.tmp.r_2.fq.gz > SLX-22602.r_2.fq.gz

nohup ./demultiplexer.rhel/demuxFQ \
    -c -d -e -i -t 1 -r 0.01 \
    -o correctFastq \
    -b SLX-22602.lostreads.r_1.fq.gz \
    -s SLX-22602.demultiplexsummary.r1.txt \
    SLX-22602.r_1.index.txt \
    SLX-22602.r_1.fq.gz &
    
nohup ./demultiplexer.rhel/demuxFQ \
    -c -d -i -e -t 1 -r 0.01 \
    -o correctFastq \
    -b SLX-22602.lostreads.r_2.fq.gz \
    -s SLX-22602.demultiplexsummary.r2.txt \
    SLX-22602.r_2.index.txt \
    SLX-22602.r_2.fq.gz &
    
rm SLX-22602.r_1.fq.gz SLX-22602.r_2.fq.gz
  
# Pipeline running
#single cell
split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/ACTTGA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/AGTCAA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/AGTTCC

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/ATGTCA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/CAGATC

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/CTTGTA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/GATCAG

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/TAGCTT

split-pipe \
    --mode comb \
    --sublibraries $PBS/newvolume/analysis/sCell/ACTTGA $PBS/newvolume/analysis/sCell/AGTCAA $PBS/newvolume/analysis/sCell/AGTTCC $PBS/newvolume/analysis/sCell/ATGTCA $PBS/newvolume/analysis/sCell/CAGATC $PBS/newvolume/analysis/sCell/CTTGTA $PBS/newvolume/analysis/sCell/GATCAG  $PBS/newvolume/analysis/sCell/TAGCTT \
    --output_dir $PBS/newvolume/analysis/sCell/combined
  
  

#human reference genome (written to '/data2/hanna/synaptogenesis/newvolume/analysis')
split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/ACTTGA_h

nohup split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/AGTCAA_h &
[1] 51084

nohup split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/AGTTCC_h &
[1] 59438

nohup split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/ATGTCA_h &
[2] 59748

nohup split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/CAGATC_h &
[3] 60324

nohup split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/CTTGTA_h &
[4] 60558

nohup split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/GATCAG_h &
[5] 60937

nohup split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/TAGCTT_h &
[6] 61208

nohup split-pipe \
    --mode comb \
    --sublibraries $PBS2/analysis/ACTTGA_h $PBS2/analysis/AGTCAA_h $PBS2/analysis/AGTTCC_h $PBS2/analysis/ATGTCA_h $PBS2/analysis/CAGATC_h $PBS2/analysis/CTTGTA_h $PBS2/analysis/GATCAG_h  $PBS2/analysis/TAGCTT_h \
    --output_dir $PBS2/analysis/combined_h &
[7] 61563
    
    

#mouse reference genome (written to '/data2/hanna/synaptogenesis/newvolume/analysis')
split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/ACTTGA_m

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/AGTCAA_m

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/AGTTCC_m

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/ATGTCA_m

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/CAGATC_m

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/CTTGTA_m

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/GATCAG_m

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_2.fq.gz \
--output_dir $PBS2/analysis/TAGCTT_m

split-pipe \
    --mode comb \
    --sublibraries $PBS2/analysis/ACTTGA_m $PBS2/analysis/AGTCAA_m $PBS2/analysis/AGTTCC_m $PBS2/analysis/ATGTCA_m $PBS2/analysis/CAGATC_m $PBS2/analysis/CTTGTA_m $PBS2/analysis/GATCAG_m  $PBS2/analysis/TAGCTT_m \
    --output_dir $PBS2/analysis/combined_m
