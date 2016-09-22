#!/usr/bin/env cwl-runner

class: CommandLineTool
description: "A Docker container for cloneHD SNV clustering. See the [cloneHD](http://www.sanger.ac.uk/science/tools/clonehd) website for more information."
id: "cloneHD"
label: "cloneHD tool"
cwlVersion: v1.0

description: |
    The cloneHD subclonal reconstruction workflow for the ICGC PanCancer Analysis of Whole Genomes (PCAWG) 
		project. For more information see the PCAWG project [page](https://dcc.icgc.org/pcawg) and our 
		GitHub [page](http://www.sanger.ac.uk/science/tools/clonehd) for our code including the source for 
		[this workflow](https://github.com/ivazquez/pcawg-clonehd-workflow).
    ```
    Usage:
    # fetch CWL
    $> dockstore tool cwl --entry quay.io/ivazquez/clonehd-pcawg:1.0-0 > Dockstore.cwl
    # make a runtime JSON template and edit it
    $> dockstore convert cwl2json --cwl Dockstore.cwl > Dockstore.json
    # run it locally with the Dockstore CLI
    $> dockstore launch --entry quay.io/ivazquez/clonehd-pcawg:1.0-0 \
        --json Dockstore.json
    ```

dct:creator:
  "@id": "http://orcid.org/0000-0003-0427-2639"
  foaf:name: Ignacio Vazquez-Garcia
  foaf:mbox: "mailto:ivg@sanger.ac.uk"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/ivazquez/clonehd-pcawg:1.0-0"
  - { import: node-engine.cwl }

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 4092
    outdirMin: 512000
    description: "the process requires at least 4G of RAM"

inputs:
  - id: "#vcf_input"
    type: File
    description: "The VCF file used as input, it must be sorted."
    inputBinding:
      position: 1
  - id: "#cna_input"
    type: File
    description: "The CNA file used as input, it must be sorted."
    inputBinding:
      position: 2

outputs:
  - id: "#clonehd_report"
    type: File
    outputBinding:
      glob: clonehd_report.zip
    description: "A zip file that contains the cloneHD report and various graphics."

baseCommand: ["bash", "/usr/local/bin/cloneHD_workflow"]