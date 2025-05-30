# graphalign

This repository contains the code for the re-implementation of the code from the article [ProGrAlign](#1)  for progressive multiple graph alignment. The code is located in the GitHub repository found [here](https://github.com/MarcosLaffitte/Progralign).

It is a part of the Masters thesis by Tobias Klink Lehn and Kasper Halkjær Beider titled Graph Alignment in Systems Chemistry. 

Graph alignment is the problem of aligning graphs in such a manner to optimise a scoring function. Our focus is on aligning multiple graphs, by computing pairwise alignments that are again graphs. As the title indicate, our primary area of interest is in chemistry, specifically aligning molecule reactions sharing the same reaction core.


## Structure of this GitHub

The over all structure of this GitHub is the following: 

* **src/** contains the source code for the implementation. 

* **test/** provides code to test the program and instances to test the program with. 

* **results/** contains pdfs displaying the results of calling the program on some of the test instances provided in the test folder. 

* **report/** contains the report for the Masters thesis. 

See the READMEs in the separate folders for their individual structures. 


## Getting started

* Clone this repository

* run ```cd graphalign/test/```

* run ```make```

Now, the executable ```align``` should be located inside ```graphalign/test/```. The file ```config.json``` contains the parameters for the progressive graph alignment. See the report for details.

Run ```align``` inside ```graphalign/test/``` as: 

```./align -{s, st, g} INSTANCE_PATH.json CONFIG_PATH.json OPTIONAL_GUIDE_TREE_STRING  ```

* ```-s``` runs the given instance without a timeout duration.
* ```-st``` runs the given instance with the timeout duration specified in ```CONFIG_PATH.json```.
* ```-g``` runs the given instance without a timeout duration and performs the alignment based on the guide tree specified by ```OPTIONAL_GUIDE_TREE_STRING```.

As an example, one can execute the following commands after having invoked ```make``` inside ```graphalign/test/```:

* ```./align -st instances/formose/Aldol\ Addition\ backward\(0\,5\).json config.json```

* ```./align -g instances/ppp/Aldolase\(17\,22\).json config.json "((((2,4),0),(3,5)),1)"```

Using the parameter ```-st``` or ```-g``` also writes the found alignment to ```INSTANCE_PATH_report.json``` which can be used for manual inspection and visualisation purposes.

## Visualisation
[ProGrAlign](#1) introduced a visualisation tool for their found alignments. We have adapted their Python code in ```test/visualisation.py``` to fit the format of our alignments. To visualise a report, one can invoke the following inside graphalign/test/:

* ```python3 visualisation.py Aldolase\(17\,22\)_report.json```

after having executed ```align``` as described in Getting Started. The visualisation creates a 1D and 2D plot of the alignment and saves them as separate PDFs.

* The 1D plot shows the alignment graph and the corresponding alignment table.

* The 2D plot shows the individual input graphs and their embeddings in the final alignment graph. Moreover, the common subgraph depicted by the set of match columns is shown on the last page.


## References
<a id="1">[1]</a> González Laffitte, Marcos E., and Peter F. Stadler. ‘Progressive Multiple Alignment of Graphs’. Algorithms, vol. 17, no. 3, 2024, https://doi.org/10.3390/a17030116.
