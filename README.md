# protein-gd
## Data analysis code for [Yu et al. (2023)](https://doi.org/10.1038/s41587-022-01598-3)
### Requirements:

Python >= 3.10

numpy

scipy

matplotlib

seaborn

tslearn 

scikit-learn 

### Installation:
Simply clone the repository, create and activate a new environment (recommended), run ```pip install -r requirements.txt```, download the data files, and run through the jupyter notebooks.

### Usage: 
[DTW_Analysis.ipynb](https://github.com/wanunulab/protein-gd/blob/master/DTW_Analysis.ipynb) contains the script for the dynamic time warping analysis of protein translocation unidirectionality.

[GBC_Analysis.ipynb](https://github.com/wanunulab/protein-gd/blob/master/GBC_Analysis.ipynb) contains the script for machine learning classification and discrimination of protein translocation events via a Gradient Boosting Classifier (GBC) model.

## Authors & Credits:
This analysis was written and performed by Ali Fallahi and Amr Makhamreh (Department of Bioengineering, Northeastern University).
The [pypore](https://github.com/jmschrei/PyPore) package was originally written by Jacob Schreiber and Kevin Karplus. We slightly modified it to make it compatible with python 3 and our datafile types.

## Downloading Data:
To replicate the analysis, download the [labeled data](https://figshare.com/s/5cd39ee415c62a316a6f?file=36587610) and the [unlabeled data](https://figshare.com/s/5cd39ee415c62a316a6f?file=36587613). All data for this work is located [here](https://figshare.com/s/5cd39ee415c62a316a6f). 

## Cite us
Please cite this work as:

Yu, L., Kang, X., Li, F. et al. Unidirectional single-file transport of full-length proteins through a nanopore. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-022-01598-3
  
