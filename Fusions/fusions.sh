#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err


COHORT=$1

index=/ceph/project/otmc/asriniva/index
genome=${index}/hg38_clean.fa
gtf=${index}/gencode46/gencode.v46.annotation.gtf
gtf2=$index/Homo_sapiens.GRCh38.104.gtf
#short=$index/gencode.v37.annotation-short.gtf
#short=$index/gencode.v38.annotation-short.gtf
short="$index/$(basename "$gtf" .gtf)-short.gtf"

root=/ceph/project/otmc/asriniva/long_read/rna/samples/$COHORT
root=/ceph/project/cribbslab/asriniva/$COHORT
cd $root

files=($root/data.dir/dRNA/*.fastq.gz)
filename=${files[$2]}
SAMPLE=$(basename "$filename")
SAMPLE="${SAMPLE%.fastq.gz}"

fq=$root/processed_fastq.dir/${SAMPLE}.filtered.fastq.gz
sam=$root/mapped_files.dir/${SAMPLE}_gene.sam
bam=${sam%.sam}.sorted.bam
bai=${bam}.bai
jaf_res="$(basename "$fq" .fastq.gz)"
fa=$root/fusions/jaffal/${SAMPLE}/${jaf_res}.fastq/${jaf_res}.fastq.fasta

##JAFFAL
export APPTAINER_CACHEDIR=/project/otmc/asriniva/mycontainers

outpath_jaf=$root/fusions/jaffal/$SAMPLE
mkdir -p ${outpath_jaf}
cd ${outpath_jaf}

singularity exec -B /ceph/package/ARCHIVE/u22 -B $root -B /ceph/project/otmc/asriniva/long_read/rna docker://trinityctat/jaffal:latest bpipe run /opt/JAFFA/JAFFAL.groovy $fq

jaf_res="$(basename "$fq" .fastq.gz)"
cp jaffa_results.csv ${jaf_res}.jaffal_results.csv
cp jaffa_results.fasta ${jaf_res}.jaffal_results.fasta

cd ${root}

##CTAT_LR
export APPTAINER_CACHEDIR=/project/otmc/asriniva/mycontainers

ctat_path=/ceph/project/otmc/asriniva/long_read/rna/fusions/ctat_lr/
outpath_ctat=$root/fusions/ctat_lr/v1.1.0_results/$SAMPLE
mkdir -p ${outpath_ctat}

singularity exec -B /ceph/package/ARCHIVE/u22 -B $root -B /usr/lib/locale/ -B /ceph/project/otmc/asriniva/ docker://trinityctat/ctat_lr_fusion:1.1.0 ctat-LR-fusion -T $fq --genome_lib_dir ${ctat_path}/ctat_genome_lib_build_dir --output ${outpath_ctat} --vis --CPU 8 --prep_reference --verbose_level 2

##FUSIONSEEKER
outpath_fs=fusions/fusionseeker/${SAMPLE}
mkdir -p ${outpath_fs}

fusionseeker_path=/ceph/project/otmc/asriniva/long_read/rna/fusions/fusionseeker_test/
export PATH=${fusionseeker_path}/FusionSeeker/:$PATH
export PATH=${fusionseeker_path}/bsalign:$PATH

fusionseeker --bam $bam -o ${outpath_fs} --datatype nanopore --ref $genome --gtf $gtf -s 3

##FLAIR_FUSION

flair_fusion=/ceph/project/otmc/asriniva/long_read/rna/fusions/flair_fusion_test/FLAIR-fusion

outpath_fl=fusions/flair_fusion/${SAMPLE}
mkdir -p ${outpath_fl}


ln -s $sam $root/${outpath_fl}/${SAMPLE}.fastq.aligned.sam -f
ln -s $fa $root/${outpath_fl}/${SAMPLE}.fastq.fasta -f
ln -s $bam $root/${outpath_fl}/${SAMPLE}.fastq.aligned.bam -f
ln -s $bai $root/${outpath_fl}/${SAMPLE}.fastq.aligned.bam.bai -f

bedtools bamtobed -bed12 -i ${outpath_fl}/${SAMPLE}.fastq.aligned.sam > ${outpath_fl}/${SAMPLE}.fastq.aligned.bed

if [ ! -f $short ]; then cd $index; python ${flair_fusion}/makeShortAnno.py $gtf; cd $root; fi
if [ ! -f ${flair_fusion}/intropolis.liftover.hg38.junctions.sorted.txt ]; then gdown https://drive.google.com/file/d/10Kz7lzVQlNF2ANoEKLcYIXPgfRubxQCQ/view?usp=sharing -O /ceph/project/otmc/asriniva/long_read/rna/fusions/flair_fusion_test/FLAIR-fusion/ --fuzzy; fi

cd ${outpath_fl}

python3 ${flair_fusion}/19-03-2021-fasta-to-fusions-pipe.py -r ${SAMPLE}.fastq.fasta -f $CONDA_PREFIX/bin/flair -g $genome -t $gtf -a $short -u -m ${SAMPLE}.fastq.aligned.bed -l 3

cd $root

##LONGGF

outpath_lg=fusions/longgf/${SAMPLE}
mkdir -p ${outpath_lg}

samtools view -S -b $sam | samtools sort -n -o ${outpath_lg}/${SAMPLE}_name_sorted.bam

LongGF ${outpath_lg}/${SAMPLE}_name_sorted.bam $gtf 100 50 100 0 1 2 64 > ${outpath_lg}/LongGF.${SAMPLE}.log && grep "SumGF" ${outpath_lg}/LongGF.${SAMPLE}.log > ${outpath_lg}/${SAMPLE}_longgf_results.tsv


##GENION
genion_path=/ceph/project/otmc/asriniva/long_read/rna/fusions/genion_test

outpath_gn=fusions/genion/${SAMPLE}
mkdir -p ${outpath_gn}

paf=${outpath_gn}/${SAMPLE}.fastq.aligned.paf
paftools.js sam2paf $sam > $paf

genion -i $fa --gtf $gtf2 --gpaf $paf -s ${genion_path}/small_example/cdna.self.tsv -d ${genion_path}/small_example/genomicSuperDups.txt -o ${outpath_gn}/${SAMPLE}_output.tsv --min-support 3 --non-coding
