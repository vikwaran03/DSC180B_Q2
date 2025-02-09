# DSC180B_Q2
Data Science Senior Capstone Project

To set up environment:

1. Make sure you have access to at least one GPU. If not, this code will take days to run.
2. Install requirements.txt.
3. If torch is already installed, uninstall it.
4. Now (re)install torch with: ```pip install torch==2.1.2 --index-url https://download.pytorch.org/whl/cu118``
5. Run ```pip show torch```
Note the torch version in your terminal output after 'Version: '. It should be something like 2.1.2+cu118
6. (Please read and review the following command carefully) Now install: ```pip install torch-scatter torch-sparse torch-cluster -f https://data.pyg.org/whl/torch-{HERE}.html``` where the {HERE} is the version from Step 4 i.e. 2.1.2+cu118.
Ideally, the command should read: ```pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-2.1.2+cu118.html```
7. Lastly install torch-geometric: ```pip install torch-geometric```


To run the code:

For the midterm checkpoint, we just have notebooks since they are easy to do experiments in them. For the final, we will consolidate these into python scripts. All the code for modeling in under Q2_prep. Navigate to that directory, all paths used later will refer to Q2_prep as the root directory.

1. graphsagenotebook.ipynb is the noteboook that has our classification with graphsage. You should be able to train a classifier that trains to 89-92% test accuracy.
2. node2vec.ipynb is the notebook that has the node2vec embeddings and then further clustering using PCA and DBSCAN.
3. eda.ipynb has some initial visualizations with our graph data.
4. VIK ADD HOW TO RUN YOUR APP