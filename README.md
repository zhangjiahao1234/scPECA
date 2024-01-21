# scPECA
## Introduction
This is a python version of PECA2 gene regulatory 
network construction software designed for single-cell data.
It has a faster running speed and a lower memory footprint.

## Installation
```commandline
pip install scPECA
```
Or install via Github
```commandline
git clone https://github.com/zhangjiahao1234/scPECA.git
cd scPECA
python setup.py install
```

## Run scPECA
### Input
scPECA requires to input the paired (sc)RNA-seq and (sc)ATAC-seq data,
and it provides data pre-processing in some formats. 
(Pre-processing can also be done manually by the user)
#### Format 1
Paired bulk RNA-seq and ATAC-seq count data

sample_name_RNA.txt

<table>
    <tr>
  		<td>gene1</td> 
        <td>10</td>
    </tr>
    <tr>
        <td>gene2</td> 
        <td>3</td>
    </tr>
</table>

sample_name_ATAC.txt


<table>
    <tr>
  		<td>chr1_10000_10100</td> 
        <td>1</td>
    </tr>
    <tr>
        <td>chr1_20000_20100 </td> 
        <td>0</td>
    </tr>
</table>

#### Format 2
scRNA-seq count data and scATAC-seq count data within the same cluster
without meta information

sample_name_scRNA.csv

|        | barcode1 | barcode2 |
|--------|----------|----------|
| gene1  | 1        | 19       |
| gene2  | 3        | 6        |

sample_name_scATAC.csv

|                  | barcode1 | barcode2 |
|------------------|----------|----------|
| chr1_10000_10100 | 1        | 0        |
| chr1_20000_20100 | 0        | 1        |

#### Format 3
scRNA-seq count data and scATAC-seq count data 
with meta information

sample_name_scRNA.csv

|        | barcode1 | barcode2 |
|--------|----------|----------|
| gene1  | 1        | 19       |
| gene2  | 3        | 6        |

sample_name_scRNA_meta.csv

<table>
    <tr>
  		<td>barcode1</td> 
        <td>celltype1</td>
    </tr>
    <tr>
        <td>barcode2</td> 
        <td>celltype2</td>
    </tr>
</table>

sample_name_scATAC.csv

|                  | barcode1 | barcode2 |
|------------------|----------|----------|
| chr1_10000_10100 | 1        | 0        |
| chr1_20000_20100 | 0        | 1        |

sample_name_scATAC_meta.csv

<table>
    <tr>
  		<td>barcode1</td> 
        <td>celltype1</td>
    </tr>
    <tr>
        <td>barcode2</td> 
        <td>celltype2</td>
    </tr>
</table>

### Run example

```commandline
import scPECA
import os

data_path = os.path.join(os.path.dirname(scPECA.__file__), 'Cones') # demo data path
demo = scPECA.scPECA(data_path, 'Cones', 'hg38') # Create a scPECA class
demo.RNA_process(2) # Format 2 RNA data processing
demo.ATAC_process(2) # Format 2 ATAC data processing
demo.network('Cones', 20, 0) # PECA2 GRN construction
```

The details of each function can be viewed in python.

### Output 
sample_name_network.txt

## System & Software Requirements
System: linux

Linux software: Homer

Python package: pybedtools, ismember, scipy, numpy_groupies, 

## Considerations

1. RNA data preprocessing currently supports only hg38 and mm10 genomes.
2. pybedtools may not support the latest version of 
python and will require the creation of a new environment.


