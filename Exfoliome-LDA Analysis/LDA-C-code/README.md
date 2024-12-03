# **LDA code for running 1, 2 and 3 feature classification & calculating bresub**

**Compiling Code on the command line:**

`gcc -o ../OneFeature One_Feature/*.c -lm`

`gcc -o ../TwoFeature Two_Feature/*.c -lm`

`gcc -o ../ThreeFeature Three_Feature/*.c -lm`

**To run the LDA Script:**
1. Copy your input data to the Data folder.
2. If you want the results to be saved to a location other than the **Results** folder, you will need to edit *outdir* variable in the *run_LDA.sh* script before running.
3. `bash run_LDA.sh`
4. After the run is completed, the 1, 2 and 3 feature classification results will be found the the **Results**.

**If you would like to change the Input or Results Directory within the LDA code**
Edit the main.c file within each folder line 38 "fp = fopen("Data/Input_Data.txt", "r");" with the folder/filename you would like to use for input data.
Edit the main.c file within each folder line 125 "fpout = fopen("Results/LDA_output_*_feature.txt", "w");;" with the folder/filename you would like to use to save results.