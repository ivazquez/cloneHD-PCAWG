# Pan-Cancer Analysis of Whole Genomes (PCAWG) - cloneHD workflow

To clone this repository, run the following command in a local directory:

    $ git clone https://github.com/ivazquez/cloneHD-PCAWG.git

[![Docker Repository on Quay](https://quay.io/repository/ivazquez/clonehd-pcawg/status "Docker Repository on Quay")](https://quay.io/repository/ivazquez/clonehd-pcawg)

## Manual

To set up the workflow in a container using `docker`:

    $ docker build -t quay.io/ivazquez/clonehd-pcawg:v1.0.0 .

This will build and compile cloneHD and cloneHD-tools, plus all dependencies.

    $ wget https://console.developers.google.com/m/cloudstorage/b/galaxyproject_images/o/Tumour2.tar.gz
    $ tar -xvfz Tumour2.tar.gz && cd Tumour2
    $ docker run -it -v `pwd`/Tumour2.mutect.vcf:/Tumour2.mutect.vcf `pwd`/Tumour2.battenberg.txt:/Tumour2.battenberg.txt ivazquez/pcawg_clonehd_workflow:v1.0-0
  
Now that you are within the `docker` container, you can execute:

    $ /usr/local/bin/cloneHD /Tumour2.vcf /Tumour2.battenberg.txt

## Automatic

Fetch a descriptor file in CWL format that tells Dockstore what are the cloneHD inputs and outputs:

    $ dockstore tool cwl --entry quay.io/ivazquez/clonehd-pcawg:v1.0-0 > Dockstore.cwl

You can create a runtime JSON template and edit it (or use the content of sample_config.json above).

    $ dockstore tool convert cwl2json --cwl Dockstore.cwl > Dockstore.json

You can now run it locally with the Dockstore CLI:

    $ dockstore tool launch --entry quay.io/ivazquez/clonehd-pcawg:v1.0-0 --json Dockstore.json

## How to cite

The cloneHD and filterHD software is free under the GNU General Public License v3. If you use this software in your work, please cite the accompanying publication:

Andrej Fischer, Ignacio Vázquez-García, Christopher J.R. Illingworth and Ville Mustonen. High-definition reconstruction of clonal composition in cancer. Cell Reports **7 (5)**, 1740-1752 (2014). [DOI: 10.1016/j.celrep.2014.04.055](http://dx.doi.org/10.1016/j.celrep.2014.04.055).