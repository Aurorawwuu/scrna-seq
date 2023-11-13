#! bin/sh 

# download SRA toolkit
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -vxzf sratoolkit.tar.gz

mkdir module  ### move SRA toolkit to module file

prefetch /dir-to-txt-contains-all-SRR-name/.txt  --max-size ...
fastq-dump --split-files --gzip --outdir /output-dir /input-dir/SRR*

# download cellranger
curl -o cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1699944730&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=gyN3dZvw1tx3l1JS2-FcLrHRbIqOgGC-WOpqdaTyhH1K7db1AOUEggiMPOLCjuz1lBppGcJvWjouCPwt9oKxCEhX6Acbtuh8V1YUQ6nLr8mRJPHYleUIDy4AcCRF0LT0-~spWF8f2Fb8L5v~H7q-CoKD2e0im-~Jwf7SbcJMuo2OFiiaWY9Y4yLz~kDVmIkWX44N2UbEIV-9JvRNPBuEl0OV4YaHv6i0r~GdodSxlopXGEcR0e5woWGCrDyM~BNWjOpp7kqWCXfYsV8U9xIEiWof1bBdmquk58fsK0IDXgR~v71tZ~jw15RKdABZTL8qOd9TNpq1PdAcfyF4VLeYRg__"
tar -vxzf cellranger-7.2.0.tar.gz
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"
tar -vxzf

cellranger count