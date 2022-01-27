## Table of contents
* [General info](#general-info)
* [System requirements](#technologies)
* [Installation guide](#setup)
* [Demo](#demo)
* [License](#license)

## General info
'scarv' is the Python package that contains functionality required to compute the Selective Constraint Against Rare SNVs (SCARV) score. This includes 1. querying variants from a VCF file, 2. training a model of mutation probabilities from the local nucleotide context, 3. genome-wide comparison of the observed and predicted genetic variation to yield the SCARV score, 4. the integration of the SCARV score into a classifier (SCARV-clf) that predicts the pathogenicity status of a non-coding variant, and 5. the querying of non-coding pathogenic variants as well as their respective SCARV and SCARV-clf scores.

## System requirements
### OS Requirements
This package is supported for macOS and Linux.
The package has been tested on the following systems:
* Linux: Scientific 7.9

### Python Dependencies
Project is created with:
* numpy
* pandas
* pyranges
* pybedtools
* xgboost
* sklearn
* itertools
* keras
* keras_genomics
* pysam
* gzip
* pyBigWig


## Installation guide
### Install from PyPi
The typical install time on a "normal" desktop computer is less than 1 minute.
```
$ pip install scarv
```

## Demo
The package contains a wide range of functions. Below we highlight and demonstrate a number of these.
### Example 1: Querying local nucleotide contexts
Pick random coordinates on chromosome 1, and random nucleotide context width (2).
Example requires fasta file for the human reference genome hg38.
Expect function to return one-hot encoded "AXGXTXTXC", "XGXTXTXCX", and "GXTXTXCXA".
Expected run time is less than 10 seconds.
```
$ import pandas as pd
$ import pyranges as pr
$ import pybedtools
$ from scarv import *

$ region_df = pd.DataFrame({'Chromosome':["chr1"], 'Start':[1e8], 'End':[1e8+2]})
$ region_gr = pr.PyRanges(region_df)

$ flank = 2
$ genome = pybedtools.genome_registry.hg38
$ reference_fasta = "hg38.fa"

$ sequence = scarv_queries.query_sequence(region_gr, flank, genome, reference_fasta)
$ print(sequence)
```

### Example 2: Extracting variants from a VCF file
Example query from gnomAD v3.0 VCF file.
To demonstrate functionality within reasonable time, focus on first 10,000 rows.
Example requires download of the relevant VCF file (https://gnomad.broadinstitute.org/downloads).
Expected run time is less than 1 second for both snippets.

First, in bash:
```
$ zcat gnomad.genomes.r3.0.sites.vcf.bgz | head -10000 | bgzip - > gnomad_vcf_first_10k_rows.vcf.bgz
```

Second, in Python:
```
$ import os
$ from scarv import *

$ variant_type = "snv"
$ ancestry = "nfe"
$ pass_snvs_query = {'variant_type': variant_type, 'ancestry': ancestry, 'PASS': True}
$ query_results = scarv_queries.query_vcf("gnomad_vcf_first_10k_rows.vcf.bgz", [pass_snvs_query])
$ print(query_results[0][:10])
```

Expected output:
[('chr1', 10113, 10114, 'T', 'C', '3', '11754'), ('chr1', 10131, 10132, 'T', 'C', '6', '20808'), ('chr1', 10137, 10138, 'T', 'C', '3', '14916'), ('chr1', 10139, 10140, 'A', 'C', '2', '33744'), ('chr1', 10139, 10140, 'A', 'G', '2', '33744'), ('chr1', 10140, 10141, 'C', 'G', '1', '25742'), ('chr1', 10143, 10144, 'T', 'C', '11', '4530'), ('chr1', 10146, 10147, 'C', 'A', '9', '242'), ('chr1', 10146, 10147, 'C', 'G', '1', '246'), ('chr1', 10148, 10149, 'C', 'A', '1', '6742')]


### Example 3: Fitting the convolutional neural network
Fitting the CNN on 1000 random training samples.
Expected run time is less than 10 seconds.
Predictions from the reverse complementary sequence should be the reverse of the vector of predictions derived from the original sequence.
```
$ import numpy as np
$ from scarv import *

$ X = np.random.random(1000*25*5).reshape(1000,25,5)
$ y = np.random.multinomial(1, [1/4.]*5, size=1000)
$ cnn = scarv_fit.fit_cnn(X, y)
$ print(cnn.summary())

$ X_test = np.random.random(2*25*5).reshape(2,25,5)
$ print('Original sequence predictions:\n', cnn.predict(X_test))
$ print('Reverse complementary sequence predictions:\n', cnn.predict(X_test[:,::-1,::-1]))
```

## License
This project is covered under the MIT License.

