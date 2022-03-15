# 4vHelix_GUI
GUI for 4vHelix

## Introduction
The program 4vHelix allows to transform any mesh structure in a .ply format to a DNA structure at the nanoscale. Using 4vHelix, you can select specific edges that you want to reinforce using a 4-helix dna bundle, which adds additional rigidity to the structure. You can then generate a cadnano .json file of the DNA structure, where your selected edges are reinforced. Finally, a sequence generating routine can be run to generate lists of nucleotide sequences for the corresponding scaffold and staple strands, so they can be ordered quickly and efficiently. 

## Usage
The program can be run using the following command:

```python
python3 4vHelix.py
```
## Input
Three files of a structure are required for edge reinforcement, these are:
- .ply
- .rpoly
- .ntrail

Note: the .rpoly and .ntrail files can be generated from a .ply file, using [bscor](https://github.com/mohamma1/bscor).

## Output
The edge reinforcement routine will generate two output files:
- <filename>.json - contains the cadnano json file of the structure including reinforcement of the selected edges.
- <filename>._virt_scaff.xlsx - contains sequences for the edges (might be deprecated).

The sequence designer/generator will generate three output files:
- scaffolds_<filename>.txt - contains the sequences of the scaffold strands. Moreover, it contains the start and end location, and the length of each scaffold.
- staples_<filename>.txt - contains the sequences of the staple strands. Moreover, it contains the start and end location, and the length of each staple.
- visualized_sequence_<filename>.txt - contains a nicely formatted visualization of the scaffold and staple sequence data, analogous to the visual representation in cadnano. This might be useful for checking the final results.

More information about the sequence generating routine can be found [here](https://github.com/ItsTheSebbe/SequenceDesigner).
  
## Example Workflow
The following is an example for workflow of reinforcing edges.
- Launch 4vHelix.py.
- Open .rpoly file.
- Open .ply file.
- Open .ntrail file.
- Select edges that need to be reinforced.
- Press "Reinforce selected edges".
- Press "Run sequence designer" and select which scaffold you want to use.
  
![Screenshot from 2022-03-14 12-23-42](https://user-images.githubusercontent.com/28595211/158162789-7e09d253-d2cb-4f33-b7b9-7173e6946193.png)
