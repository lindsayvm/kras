#!/usr/bin/env bash
#Download TCGA data with token key

# Go to https://portal.gdc.cancer.gov/repository

#Select inclusion criteria:
  #Colon adenomas and adenocarcinomas TCGA
  #RNAseq counts
  #Simple nucleotide variation annotated somatic mutation WXS vcf

#Download manifest

#From terminal in own laptop:
rsync -avhzP -e 'ssh -o "ProxyCommand ssh -A l.leek@rhpc.nki.nl -W %h:%p"' ~/Downloads/gdc* l.leek@darwin:/home/l.leek/data/triflu/portal_GDC

#From terminal in darwin:
cd /home/l.leek/data/triflu/portal_GDC

#Secure the token (only readable and writeable for owner)
chmod 600 /home/l.leek/data/triflu/portal_GDC/gdc-user-token.txt

cd ..

# go to https://gdc.cancer.gov/access-data/gdc-data-transfer-tool and use link of Ubunu 64
wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip
unzip gdc-client_v1.6.1_Ubuntu_x64.zip
rm gdc-client*.zip
#Check whether it is working
./gdc-client

#Download secured data
./gdc-client download -m /home/l.leek/data/triflu/portal_GDC/gdc_manifest_WES.txt -t /home/l.leek/data/triflu/portal_GDC/gdc-user-token.txt

rm [!gdjgfdshtrszjtsxj]
