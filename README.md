# DSC180B_Q2
Data Science Senior Capstone Project

To set up environment:

1. Make sure you have access to at least one GPU. If not, this code will take days to run.
2. Install requirements.txt.
3. If torch is already installed, uninstall it.
4. Now (re)install torch with: ```pip install torch==2.1.2 --index-url https://download.pytorch.org/whl/cu118```
5. Run ```pip show torch```
Note the torch version in your terminal output after 'Version: '. It should be something like 2.1.2+cu118
6. (Please read and review the following command carefully) Now install: ```pip install torch-scatter torch-sparse torch-cluster -f https://data.pyg.org/whl/torch-{HERE}.html``` where the {HERE} is the version from Step 4 i.e. 2.1.2+cu118.
Ideally, the command should read: ```pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-2.1.2+cu118.html```
7. Lastly install torch-geometric: ```pip install torch-geometric```


To run the code:

Navigate to the directory Modeling, all paths used later will refer to Modeling as the root directory. This is where all the code lies. All the notebooks/scripts in this directory (not sub-directories) are all you need to replicate our analysis. 1 and 2 are the most vital deliverables of this project - you can run this code through without any additional setup.

1. node2vec.ipynb contains the notebook that does our whole clustering pipeline. Run this notebook through to replicate the clustering analysis including node2vec embeddings, clustering, statistical tests, and HI-C/gene analysis.
2. graphsage.py contains the code needed to run our classification model. You will be able to achieve an accuracy of around ~81%. Run this script with ```python graphsage.py```. 
3. structure.ipynb has the code for the 3d visualizations of ecDNA and HSR. All the different iterations of this plot, including gene labels, cluster labels, graphs, are in this notebook.
4. eda.ipynb has some initial visualizations and the preprocessing associated with our graph data.
5. data_prep.ipynb is the notebook where we prune our adjacency matrix based on Euclidean distances.
