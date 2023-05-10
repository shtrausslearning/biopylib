![](https://i.imgur.com/zX8Qn0i.png)

#### **Project Information**
![](https://camo.githubusercontent.com/d38e6cc39779250a2835bf8ed3a72d10dbe3b05fa6527baa3f6f1e8e8bd056bf/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f436f64652d507974686f6e2d696e666f726d6174696f6e616c3f7374796c653d666c6174266c6f676f3d707974686f6e266c6f676f436f6c6f723d776869746526636f6c6f723d326262633861) ![](https://badgen.net/badge/status/WIP/orange) 

- The project is currently at a stage where the bio modules are being added

#### **This package aims to:**
- Minimise/remove the need for coding that required from the users end for a specific biological analysis process
- This is achieved by interpreting & analysing the input (**desired series of operations**) using **NPL** methods

#### Change Log:
- <code>0.0.8</code> : Added **Deterministic Motif Discovery** Classes
- <code>0.0.7</code> : Added **GenBank** reader

#### **Kaggle** | **Github** version: 

- **<code>[biopylib](https://www.kaggle.com/datasets/shtrausslearning/biopylib)</code>** **0.0.7** | **<code>[biopylib](https://github.com/shtrausslearning/biopylib)</code>** **0.0.8**

#### pypi version
[![PyPI version](https://badge.fury.io/py/biopylib.svg)](https://badge.fury.io/py/biopylib)

#### **src** contents:

biopylib

- <code>blast</code> - BLAST sequence query operations
- <code>read_sequence</code> - Biological sequence file read/write operations
- <code>sequence</code> - Biological sequence operations
- <code>sequence_alignment</code> - Biological sequence alignment operations
- <code>sce</code> - Single Cell Experiment Analysis (SCE) operations
- <code>motif</code> - Deterministic & Probabilistic Motif Discovery operations

#### **biopylib notebooks:**

You can fork and test the latest version on Kaggle:

- **[biopylib playground](https://www.kaggle.com/code/shtrausslearning/biopylib-playground)**

Notebooks relating to **biopylib**:

- **[Biological Sequence Operations](https://www.kaggle.com/code/shtrausslearning/biological-sequence-operations)** (operations usising **blast**, **read_sequence**, **sequence** classes)
- **[Biological Sequence Alignment](https://www.kaggle.com/code/shtrausslearning/biological-sequence-alignment)** (operations using **read_sequence**, **sequence_alignment** classes)

#### **Prepare whl file using:**

```
python setup.py bdist_wheel --universal
```

#### **Install in notebook via:**

```
!pip install /path/
```

#### **Upload to pip:**

```
!pip install twine
twine upload biopylib-0.0.5-py2.py3-none-any.whl
```
