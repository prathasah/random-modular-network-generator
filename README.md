Random Modular Network Generator
================================
This is the source code used for the following paper: 

Sah, Pratha*, Lisa O. Singh, Aaron Clauset, and Shweta Bansal. "Exploring community structure in biological networks with random graphs." BMC Bioinformatics 2014, 15:220. doi:10.1186/1471-2105-15-220

This Python script generates undirected, simple, connected graphs with a specified degrees and pattern of communities, while maintaining a graph structure that is as random as possible.

email: ps875@georgetown.edu, sb753@georgetown.edu

Please cite the paper above, if you use our code in any form or create a derivative work.

Sample Output Graph
================================

![alt tag](https://github.com/prathasah/random-modular-network-generator/blob/master/Figure2.jpg)

Modular random graphs with network size = 150, 375 edges, 3 modules (or communities) of size, *s* = 50 and degree distribution is power law with modularity values of: a) *Q* = 0.1; b) *Q* = 0.3; and c) *Q* = 0.6. 

From  Figure 2, Sah *et al.* (2014)

Dependencies
================================
Networkx Python package is required. Link to install Networkx can be found [here](https://networkx.github.io/). The generator also requires sequence_generator.py script which is provided.

Usage
================================

For a quick demo of the code, open the script "random_modular_generator_variable_modules.py" and scroll down to the *main* function. Adjust the model parameters, degree distribution function and module size distribution function and run in terminal using the command:

`$ python random_modular_generator_variable_modules.py`



Sample Code
================================
For users unfamiliar to Python, I have uploadeded a sample code file (mock_code.py) demonstrating how the graph generator can be imported and used in a script. The mock code can be run using the command

`$ python mock_code.py`





