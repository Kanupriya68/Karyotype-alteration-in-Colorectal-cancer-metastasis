
import pandas as pd
import os

import warnings

warnings.filterwarnings("ignore")

tcgaData = pd.read_csv("C:\\Research\\Research_project\\tcga\\Lung.csv", sep=',')


tcgaData["Chromosome"] = tcgaData["Chromosome"].astype(str).str.split('.', expand = True)[0]

baseHartwigDf = pd.read_csv("C:\\Research\\Research_project\\tcga\\dataCN.csv", sep=';', header=0)
baseHartwigDf = baseHartwigDf[baseHartwigDf.columns[0:5]]
baseHartwigDf["chromosome"] = baseHartwigDf["chromosome"].astype(str)

tSampleIds = tcgaData.Sample.unique()

#tcgaSampleIds = ["TCGA-WD-A7RX"]

t1 = tSampleIds[0:950]

j=0

for i in t1:
        currentId = i
        j=j+1

        tcgabyID = tcgaData[tcgaData.Sample == currentId]

        print("processing: "+ currentId)
        print("processing file no:"+ str(j) )


        def getCnNumber(start,end,chr):
                    tcgaData_chr = tcgabyID[(tcgaData.Chromosome == chr)]
                    baseInsideTcgaBin = tcgaData_chr[((tcgaData_chr["Start"] <= start) & (tcgaData_chr["End"] >= end))]

                    if len(baseInsideTcgaBin.index) > 0:
                        return baseInsideTcgaBin.head(1)["Modal_Total_CN"].iloc[0]
                    else:
                            baseLeftTcgaBin = tcgaData_chr[
                                ((tcgaData_chr["Start"] >= start) & (tcgaData_chr["Start"] <= end) &
                                 (tcgaData_chr["End"] >= start) & (tcgaData_chr["End"] >= end))]

                            if len(baseLeftTcgaBin.index)>0:
                                lengthBinInsideTcgaBin = end - baseLeftTcgaBin.head(1)["Start"].iloc[0]
                                if lengthBinInsideTcgaBin >= 500000:
                                    return baseLeftTcgaBin.head(1)["Modal_Total_CN"].iloc[0]
                                else:
                                    return "NA"

                            else:
                                    baseRightTcgaBin = tcgaData_chr[
                                        ((tcgaData_chr["Start"] <= start) & (tcgaData_chr["Start"] >= end) &
                                         (tcgaData_chr["End"] >= start) & (tcgaData_chr["End"] >= end))]



                                    if len(baseRightTcgaBin.index)>0:



                                        lengthBinInsideTcgaBinR = baseRightTcgaBin.head(1)["End"].iloc[0] - start



                                        if lengthBinInsideTcgaBinR >= 500000:
                                            return baseRightTcgaBin.head(1)["Modal_Total_CN"].iloc[0]
                                        else:
                                            return "NA"

                                    else:

                                        baseOutSideTcgaBin = tcgaData_chr[
                                            ((tcgaData_chr["Start"] >= start) & (tcgaData_chr["End"] <= end))]

                                        if len(baseOutSideTcgaBin.index)>0:


                                            lengthBinInsideTcgaBinI = baseOutSideTcgaBin.head(1)["End"].iloc[0] - \
                                                                      baseOutSideTcgaBin.head(1)["Start"].iloc[0]

                                            if lengthBinInsideTcgaBinI >= 500000:
                                                return baseOutSideTcgaBin.head(1)["Modal_Total_CN"].iloc[0]
                                            else:
                                                return "NA"
                                        pass


        cNList = []

        for n in range(0, len(baseHartwigDf)):
            cNList.append(getCnNumber(baseHartwigDf.start[n],baseHartwigDf.stop[n], baseHartwigDf.chromosome[n]))

        print(cNList)
        print(len(cNList))
        baseHartwigDf[i] = cNList



#print(baseHartwigDf)




output_folder = "C:\\Research\\Research_project\\TCGAbinned"
output_path = os.path.join(output_folder,"Lungoutput.csv")

baseHartwigDf.to_csv(output_path, sep=',', index=False,mode='a')




