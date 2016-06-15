## rMATS-STAT: The Stand-Alone rMATS Statistical Model

Requirements
------------
1. Install Python 2.6.x or Python 2.7.x and corresponding versions of NumPy and
SciPy.
2. Add the Python directory to the $PATH environment variable.

Installation:
------------
The source code can be directly called from Python.

Usage:
--------------------------------
Use on replicates without pairs:

$ python rMATS_unpaired.py input_read_count_file output_folder number_processor diff_cutoff

Use on replicates with pairs (Each replicate is paired with another between the two sample groups):

$ python rMATS_paired.py input_read_count_file output_folder number_processor diff_cutoff

The examples of input files are available in the depositry. The input format is detailed in the Example section.

Example:
--------------------------------
Test rMATS statistical model on sample input:
Run rMATS_unpaired.py as below to test the script.

    $ python rMATS_unpaired.py inc.txt ./ 1 0.1

Input: One input file is required for the script. inc.txt: Each row of
this input file contains an alternative splicing event. The 5 columns of this
input file contain the read counts for the two isoforms of the alternative
splicing event in the two sample groups (2-3 columns are for the first group, 4-5 columns are for the second group). The read counts for different patients are separated by commas
in the column. As an example, for exon skipping events, each row defines a
skipped exon and the columns contain the read counts for inclusion and skipping
isoforms:
- ID: User defined ID for the alternative splicing event.
- IJC1: inclusion junction counts for group 1, patients are separated by comma.
- SJC1: skipping junction counts for group 1, patients are separated by comma.
- IJC2: inclusion junction counts for group 2, patients are separated by comma.
- SJC2: skipping junction counts for group 2, patients are separated by comma.
- IncFormLen: length of inclusion form, used for normalization.
- SkipFormLen: length of skipping form, used for normalization.

Output: The output folder contains the rMATS_Result_P.txt file. For each alternative splicing event, rMATS outputs the P-values that evaluate the associations between alternative splicing and patient survival.
- ID: User defined ID for the alternative splicing event.
- IJC1: inclusion junction counts for group 1, patients are separated by comma.
- SJC1: skipping junction counts for group 1, patients are separated by comma.
- IJC2: inclusion junction counts for group 2, patients are separated by comma.
- SJC2: skipping junction counts for group 2, patients are separated by comma.
- IncFormLen: length of inclusion form, used for normalization.
- SkipFormLen: length of skipping form, used for normalization.
- PValue: P-values of the alternative splicing event.


Contacts and bug reports
------------------------
Shihao Shen
shihao@ucla.edu

If you found a bug or mistake in this project, we would like to know about it.
Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
   fixed.
2. Check that your input is in the correct format and you have selected the
   correct options.
3. Please reduce your input to the smallest possible size that still produces
   the bug; we will need your input data to reproduce the problem, and the
   smaller you can make it, the easier it will be.

Publication
------------

Copyright and License Information
---------------------------------
Copyright (C) 2015 University of California, Los Angeles (UCLA)
Shihao Shen and Yi Xing

Authors: Shihao Shen and Yi Xing

This program is licensed with commercial restriction use license. Please see the attached LICENSE file for details.

