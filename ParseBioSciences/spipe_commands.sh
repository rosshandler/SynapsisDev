# spipe.v1.0.0 is already installed as conda env
conda activate spipe.v1.0.0

PBS='/data2/ivanir/Feline2023/ParseBS'

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
cat SLX-22602.DNAA007.HGMLNDMXY.s_1.r_1.fq.gz SLX-22602.DNAA007.HGMLNDMXY.s_2.r_1.fq.gz > SLX-22602.r_1.fq.gz
cat SLX-22602.DNAA007.HGMLNDMXY.s_1.r_2.fq.gz SLX-22602.DNAA007.HGMLNDMXY.s_2.r_2.fq.gz > SLX-22602.r_2.fq.gz

# Pipeline running
#single cell
split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10/ \
--fq1 $PBS/newvolume/expdata/SLX-22602.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/SLX-22602.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell
