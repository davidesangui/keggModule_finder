# keggModule_finder
Small script to parse eggNOG-mapper annotations and KEGG database to identify KEGG modules encoded by microbial genomes.
## Installation
No installation is needed. Simply download the file keggModule_finder.py from this repository and launch it with Python3.  
If the KEGG database is to be parsed to obtain module definitions, then the python package Bio is required, which can be obtained via pip:
```
pip install bio
```
## How it works
keggModule_finder uses the rules of KEGG module definitions and the KEGG orthologues identified by eggNOG-mapper to define the capability of the corresponding microbial genomes to encode KEGG modules.  
Each definition is divided in step, and each step requires one or more KEGG orthologues depending on a set of logical relationship between them. keggModule_finder defines a KEGG module as encoded by a microbial genome if at least N-1 steps can be performed by the genome, where N is the total number of steps of the module.  
A binary matrix defining presence and absence of KEGG modules in a genome is given as output, as well as the count of encoded steps for each genome for each module.
Starting from the count output, any user can define the encoding potential in the way he finds most suitable.
## Usage
keggModule_finder.py requires eggNOG-mapper annotations. In particular, the file format ".emapper.annotations" is required. keggModule_finder.py takes as input a folder containing these annotation files for any number of genomes. Only the files ending with the suffix ".annotations" will be evaluated.  
keggModule_finder.py requires a "module_definition" file, as the one in this repository. Only the modules listed in this file will be evaluated. This file can be passed with the flag ```-definition_file```. A definition file listing all KEGG modules (July 2024) can be obtained from this repository. A more updated version can be obtained by parsing the KEGG database by setting the flag ```--get_definitions``` (which takes five-ten minutes). Such a file can be modified by inserting custom modules, as long as the custom definitions observe KEGG definitions' logic and format.
### Example
```
python3 keggModule_finder.py annot_folder output -definition_file module_definitions_JULY24.tsv
```
Two output files are generated. One (output_count.tsv) represents the count of steps for each genome for each module:
```
genome	M00001,length=9	M00002,length=5	M00003,length=7 [...]
genome1	9	5	7   [...]
genome2	8	4	6   [...]
genome3	9	5	6   [...]
```
The other (output_minusOne.tsv) is a binary matrix where 1 means that the corresponding genome can encode at least N-1 steps of the corresponding module, where N is the total number of steps of the module (as shown in the previous output).
```
genome	M00001 M00002	M00003 [...]
genome1	1	1	1   [...]
genome2	1	1	1   [...]
genome3	1	1	1   [...]
```
