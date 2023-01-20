# spipe.v1.0.0 is already installed as conda env
conda activate spipe.v1.0.0

PBS='/data2/ivanir/Feline2023/ParseBS'

PATH="/data2/ivanir/Feline2023/ParseBS/ParseBiosciences-Pipeline.1.0.3p:$PATH"

# Genome indexing
split-pipe \
--mode mkref \
--genome_name hg38 mm10 \
--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz $PBS/newvolume/genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
--genes $PBS/newvolume/genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz $PBS/newvolume/genomes/Mus_musculus.GRCm38.93.gtf.gz \
--output_dir $PBS/newvolume/genomes/hg38_mm10

split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--genes $PBS/newvolume/genomes/Homo_sapiens.GRCh38.108.gtf.gz \
--output_dir $PBS/newvolume/genomes/hg38

split-pipe \
--mode mkref \
--genome_name mm10 \
--fasta $PBS/newvolume/genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
--genes $PBS/newvolume/genomes/Mus_musculus.GRCm38.93.gtf.gz \
--output_dir $PBS/newvolume/genomes/mm10
