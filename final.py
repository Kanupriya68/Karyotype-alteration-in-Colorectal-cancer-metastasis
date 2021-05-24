##Import pandas
import pandas as pd
import math
import os
import openpyxl
import numpy as np

#Create function to calculate 1MB bins
def createbins(file,filepath):
    #Read files one by one and save in a variable
    somatic_cnv_original = pd.read_csv(filepath, sep="\t", header= (0))
    #Name columns in each file
    somatic_cnv_original= somatic_cnv_original[['chromosome','start','end','copyNumber']]

#Create empty list to append chromosome numbers
    chromosomeList = []
    #Append autosome chromosomes one by one in the loop
    for i in range(1,23):
        chromosomeList.append(str(i))

    #Append sex chromosomes as well
    chromosomeList.append("X")
    chromosomeList.append("Y")

    #Create empty list
    appendedData = []

    #Create a loop to read every chromosome from chromosome list and create bins
    for w in chromosomeList:

        somatic_cnv_chr = somatic_cnv_original[somatic_cnv_original.chromosome==w]
        #Create condition for bin size
        a = somatic_cnv_chr.shape[0]
        if  a <1:
            break
        else:

            bin_size = 1000000
            chromosomeName = w
            chr_stop = somatic_cnv_chr["end"]
            max_chr_stop = chr_stop.max()
            chr_bins = max_chr_stop/bin_size
            chr_bins = math.ceil(chr_bins)


            #Create an Empty database
            binList = []
            chromosome = []
            start = []
            end = []

            for i in range(1,chr_bins+1):
                binList.append(i)

            for a in range(1,chr_bins+1):
                chromosome.append(chromosomeName)

            for n in range(1,chr_bins+1):
                start.append((n-1)*bin_size+1)

            for n in range(1,chr_bins+1):
                end.append(n*bin_size)
            #Create dictionary
            data = {'binsNumber': binList, 'chromosome': chromosome, 'start': start, 'stop': end}
            #Insert new dataframe into blank dataframe
            blankDf = pd.DataFrame(data, columns= ['binsNumber','chromosome','start','stop'])

            #Create function to read copy number values bin wise
            def getCnNumber(start,end):
                try:
                    jkl = somatic_cnv_chr[((somatic_cnv_chr["start"] <= start) & (somatic_cnv_chr["end"] >= end))]

                    #print(len(jkl.index))
                    if len(jkl.index) < 2:
                        return jkl.head(1)["copyNumber"].iloc[0]

                except:
                    #print("Start: "+str(start)+" | End: "+str(end))
                    #print(somatic_cnv_chr[((somatic_cnv_chr["start"]>=start)& (somatic_cnv_chr["start"]<=end)) | ((somatic_cnv_chr["end"]>=start)& (somatic_cnv_chr["end"]<=end))])
                    b = somatic_cnv_chr[((somatic_cnv_chr["start"]>=start)& (somatic_cnv_chr["start"]<=end)) | ((somatic_cnv_chr["end"]>=start)& (somatic_cnv_chr["end"]<=end))]
                    #print(b)
                    b['diff'] = b['end']-b['start']
                    listG = b['diff'].tolist()
                    #print("FailA")
                    #print(listG)
                    #print("FailB")


                    List2 = listG[1:-1]

                    #print(b.head(1)["end"].iloc[0])
                    c = b.head(1)["end"].iloc[0]
                    #print(b.tail(1)["start"].iloc[0])
                    d = b.tail(1)["start"].iloc[0]
                    y = c-start
                    z = end-d
                    #print(y)
                    #print(z)


                    j = len(b.index)
                    lo=[]
                    lo.append(y)
                    lo = lo + List2
                    lo.append(z)
                    #print(lo)

                    #print(b)

                    weights=[]
                    for i in lo:
                        weights.append(i/bin_size)

                    #b["weights"]=weights
                    #print(b)

                    listCN = b['copyNumber'].tolist()
                    #print(listCN)

                    CN = 0
                    for i in range(0,len(listCN)):
                        CN += weights[i]*listCN[i]

                    #print(CN)
                    return CN

            #getCnNumber(5000001,6000000)


            cNList = []

            for n in range(1,chr_bins+1):
                cNList.append(getCnNumber((n-1)*bin_size+1,n*bin_size))

            #print(cNList)
            #print(len(cNList))
            blankDf['copyNumber'] = cNList
            #print(blankDf.head)

            appendedData.append(blankDf)


    appended_data = pd.concat(appendedData)
    output_folder = "C:\Research\Research_project\indexCN"
    output_path = os.path.join(output_folder,file)
    #appended_data.to_excel('C://Research/Research_project/CPCT02170013T-binned.xlsx')
    appended_data.to_csv(output_path, sep=',', index=True,mode='a')

    #np.savetxt('CPCT02170013T-binned.txt', appended_data.values, fmt='%d', delimiter="\t")


createbins("WIDE01010942T.purple.cnv.somatic.tsv","C:/Research/Research_project/missingcnv/200901_HMFregWIDE_FR16616764_FR16612260_WIDE01010942/WIDE01010942T.purple.cnv.somatic.tsv")