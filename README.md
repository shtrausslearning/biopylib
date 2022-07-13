![](https://i.imgur.com/x1chlWl.png)

#### **Project Information**
![](https://camo.githubusercontent.com/d38e6cc39779250a2835bf8ed3a72d10dbe3b05fa6527baa3f6f1e8e8bd056bf/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f436f64652d507974686f6e2d696e666f726d6174696f6e616c3f7374796c653d666c6174266c6f676f3d707974686f6e266c6f676f436f6c6f723d776869746526636f6c6f723d326262633861) ![](https://badgen.net/badge/status/WIP/orange) 

#### **This package aims to:**
- Minimise/remove the need for coding that required from the user end for a specific biological analysis process
- This is achieved by interpreting & analysing the input (**desired series of operations**) using **NPL** methods

### Change Log:
- <code>0.0.7</code> : Added **Deterministic Motif Discovery** Classes
- <code>0.0.6</code> : Added **GenBank** reader

#### Kaggle version: 
- **<code>[biopyseq](https://www.kaggle.com/datasets/shtrausslearning/biopylib)</code>**

#### pypi version
[![PyPI version](https://badge.fury.io/py/biopylib.svg)](https://badge.fury.io/py/biopylib)

#### **src** contents:
- <code>blast</code> - BLAST sequence query operations
- <code>read_sequence</code> - Read biological sequence file formats
- <code>sequence</code> - Biological sequence operations
- <code>sequence_alignment</code> - Biological sequence alignment operations
- <code>sce</code> - Single Cell Experiment Analysis (SCE) operations
- <code>dmotif</code> - Deterministic Motif Discovery

#### **forkable testing notebook:**
- **[biopylib playground](https://www.kaggle.com/code/shtrausslearning/biopylib-playground)**

#### **Prepare whl file using:**

- python setup.py bdist_wheel --universal

#### **Install in notebook via:**

- !pip install /path/

#### **Upload to pip:**

- !pip install twine
- twine upload biopylib-0.0.5-py2.py3-none-any.whl
