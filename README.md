# MutantBuilder

Modeller and python3 are needed to run. 

PARAMETERS:

    -pdb: PDB file 
    -c: Chain to mutate
    -f: Config file containing which positions to mutate (defualt 'config.txt')
    -o: Output folder (defualt 'mutants')
    -r: add flag to not use homology modeling and use restricted model instead 

Example use: 

python mutantBuilder.py -pdb 6waq.pdb -c A

CONFIG FILE:

- The config.txt is used to guide the mutagensis each line contains a residue postion and amino acid to mutate to. 
- A '*' can be used to mutate to all resdiues. 
- Multiple mutants can be done at the same time using a comma as a seperator. (on the same line)

        EX 1: 2 seperate single mutants. 
        1 ALA
        4 GLN 

        EX 2: 1 double mutant.
        1 ALA, 4 GLN 

        EX 3: All mutants at postion 5
        5 *

