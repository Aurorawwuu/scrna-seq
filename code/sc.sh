#! bin/sh 

### download SRA toolkit
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -vxzf sratoolkit.tar.gz

mkdir module  # move SRA toolkit to module file

module/sratoolkit/bin/prefetch --option-file srr.txt --max-size 60g
nohup module/sratoolkit/bin/fastq-dump --split-files --gzip --outdir scrna-seq/data/ SRR* &

### download cellranger and 10x refernce genome
curl -o cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1699944730&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=gyN3dZvw1tx3l1JS2-FcLrHRbIqOgGC-WOpqdaTyhH1K7db1AOUEggiMPOLCjuz1lBppGcJvWjouCPwt9oKxCEhX6Acbtuh8V1YUQ6nLr8mRJPHYleUIDy4AcCRF0LT0-~spWF8f2Fb8L5v~H7q-CoKD2e0im-~Jwf7SbcJMuo2OFiiaWY9Y4yLz~kDVmIkWX44N2UbEIV-9JvRNPBuEl0OV4YaHv6i0r~GdodSxlopXGEcR0e5woWGCrDyM~BNWjOpp7kqWCXfYsV8U9xIEiWof1bBdmquk58fsK0IDXgR~v71tZ~jw15RKdABZTL8qOd9TNpq1PdAcfyF4VLeYRg__"
tar -vxzf cellranger-7.2.0.tar.gz
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"
tar -xzf refdata-gex-GRCh38-2020-A.tar.gz

### rename file to 10x format
cp scrna-seq/data/SRR* scrna-seq/copy/
cd scrna-seq/data/
for i in *_1*
do
    new=$(cut -d "_" -f1 <<< "$i")
    # echo "${new}_S0_L001_I1_001.fastq.gz"

    mv "$i" "${new}_S0_L001_I1_001.fastq.gz"
done

for i in *_2* ; do
    new=$(cut -d "_" -f1 <<< "$i")
    mv "$i" "${new}_S0_L001_R1_001.fastq.gz"
done

for i in *_3* ; do
    new=$(cut -d "_" -f1 <<< "$i")
    mv "$i" "${new}_S0_L001_R2_001.fastq.gz"
done

### cellranger count
mapfile -t files <  scrna-seq/srr.txt 
echo "${files[@]}"
for i in ${files[@]}
do
    module/cellranger/bin/cellranger count --id $i \
    --transcriptome reference/hg38_10x/ \
    --fastqs scrna-seq/data/ \
    --sample $i \
    --output-dir scrna-seq/output/$i \
    --localcores 36
done