import os
import pandas as pd
from final import createbins

#Read Hartwigs metadata with cancer patient info
meta = pd.read_csv("C:\Research\Research_project\Colcc_metadata.csv", sep=';',header= 0)

#Create list with all patient sample IDs
list_of_files = meta['sampleId'].to_list()

#Create an empty list
new_list = []

#Add suffix in every sample ID name to make it compatible for further reading
for i in list_of_files:
    new_list.append(i+".purple.cnv.somatic.tsv")

#Loop which walks in every directory and creates a list of file and its filepath
for root, dirs, files in os.walk("C:\Research\Research_project\rawfiles"):
    print(files)
    for file in files:
        filepath = os.path.join(root, file)
        print(file)
        for i in new_list:
            if i == file:
                createbins(file,filepath)
print(new_list)