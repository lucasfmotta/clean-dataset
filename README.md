# clean-dataset
Clean a dataset

## Dependencies

Biopython & tqdm

```
pip install biopython
pip install tqdm
```

## Installation

```
sudo wget https://raw.githubusercontent.com/lucasfmotta/clean-dataset/main/clean-dataset.py -O /usr/local/bin/clean-dataset
sudo chmod a+rx /usr/local/bin/clean-dataset
```

## Usage
```
clean-dataset [Blast-parameters]
  [optional]
```
-h, -help: prints help
-U: updates to newest version
Blast parameters: you can use every BlastN parameters except -db -query or -out (for more inforomation read [Blast documentation](https://www.ncbi.nlm.nih.gov/books/NBK279684/#_appendices_Options_for_the_commandline_a_))
