# 4vHelix_GUI
GUI for 4vHelix

## Overview
This program allows you to load any wireframe structure in ply format and select edges that you want to reinforce using a 4-helix dna bundle. It will then generate a cadnano json file where the are selected edges reinforced. A sequence generating routine can then be run to generate a list sequences for the scaffold and staple strands. 

## Usage
The program can be run using the following command:

```python
python3 GUI.py
```
## Input
Three files of a structure are required for edge reinforcement, these are:
- .ply
- .rpoly
- .ntrail

## Output
The edge reinforcement routine will generate two output files:
- <filename>.json - contains the cadnano json file of the structure including reinforcement of the selected edges.
- <filename>._virt_scaff.xlsx - contains sequences for the edges (might be deprecated).

The sequence designer/generator will generate three output files:
- scaffolds.txt - contains the sequences of the scaffold strands. Moreover, it contains the start and end location, and the length of each scaffold.
- staples.txt - contains the sequences of the staple strands. Moreover, it contains the start and end location, and the length of each staple.
- visualized_sequence.txt - contains a nicely formatted visualization of the scaffold and staple sequence data, analogous to the visual representation in cadnano. This might be useful for checking the final results.

More information about the sequence generating routine can be found [here](https://github.com/ItsTheSebbe/SequenceDesigner).
  
## Example Workflow
Here is an example for workflow of reinforcing edges.
- Launch GUI.py.
- Open rpoly file.
- Open ply file.
- Open ntrail file.
- Select edges that need to be reinforced.
- Press "Reinforce selected edges".
- Press "Run sequence designer".
  
![Screenshot from 2022-03-14 12-23-42](https://user-images.githubusercontent.com/28595211/158162789-7e09d253-d2cb-4f33-b7b9-7173e6946193.png)
