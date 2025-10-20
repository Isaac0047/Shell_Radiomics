# Shell_Radiomics
This repository stores the code for shell radiomics extraction and analysis

Step 1: shell_worldmap_visualize_model.py generates shells from tumor and visulize the shell contours.

tumor_visualize.m is a Matlab code for shell contour visualization.

Step 2: feature_extraction_shell.py extracts radiomics features from shells and save to excel format.

Step 3: radiogenomic_2d_slice.py, radiogenomic_3D.py and shell_radiomic_NN.py conduct machine learning predictions with different types of radiomic extraction methods.

Step 4: shell_radiomic_NN_extend.py includes the additional analysis after the machine learning model is fully trained.


