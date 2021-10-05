#Author: Abrar Al-Shaer
#Date: 8/29/19
#This program filters unwanted PCR smaple results and averages technical replicates within the file from the mean and standard deviation (blank cells were removed in excel).
#After filteration the program calculates the Ct(also known as Cq) value for each sample normalized to the reference, fold change, SEMs, and averaged standard deviations.

#!/usr/bin/env python

#Initializing imported packages
import pandas as pd
from pandas import DataFrame
import numpy as np
from scipy import stats
import os, sys

reference = input("What is your reference? Make sure to include quotes: ")
print "Your Reference: ", reference
one_target = input("What is your target? Make sure to include quotes: ")
print "Your target: ", one_target

#FILE INPUTS
print "File Inputs:"
PCR_file1 = 'Blk1_4.AvUCP.PreHeat.Sex.Diet.Inc_CLEAN_Female.csv'
print "FILE 1:", PCR_file1
PCR_file2 = 'Blk5_8.AvUCP.PreHeat.Sex.Diet.Inc_CLEAN_Female.csv'
print "FILE 2:", PCR_file2
PCR_file3 = 'Blk9_12.AvUCP.PrePostHeat.Sex.Diet.Inc_CLEAN_Female.csv'
print "FILE 3:", PCR_file3
PCR_file4 = 'Blk13-16.AvUCP.PostHeat.Sex.Diet.Inc_CLEAN_Female.csv'
print "FILE 4:", PCR_file4
PCR_file5 = 'Blk16-20.AvUCP.PostHeat.Sex.Diet.Inc_CLEAN_Female.csv'
print "FILE 5:", PCR_file5

def csv_init(file_PCR):
    print "---------Initiation of CSV file & Filteration---------"
    columns = [] #flag = columns
    from_csv = pd.read_csv(file_PCR, sep=None, header=0, engine='python') #reading in the CSV file
    columns = from_csv[['Target','Content','Sample','Biological Set Name','Cq','Cq Mean','Cq Std. Dev']] #reads header columns
    #filteration step
    #blank Cq rows (or NTC) were removed in excel via F5 command
    #df = dataframe
    df = columns[columns.Cq != 0] #removes zeroes if there is a zero in the Cq column
    df = df.drop_duplicates('Content', keep='first') #removing duplicates based on content column value
    print df #print first 5 lines of dataframe
    return df #function returns dataframe as output

def rows_init_store(dataframe):
    print "---------Defining and Storing Rows From Dataframe---------"
    #Initializing blank lists
    Cq = []
    CqDev = []
    CqAvg = []
    content = []
    target = []
    sample = []
    set_name = []
    #loop through all the rows in the dataframe and add each element of each row to the proper column list
    for index, row in dataframe.iterrows():
        print row[1], row[0], row[5], row[6] #print content variable, target, mean, and standard deviation
        target.append(row[0]) #append = add to list
        content.append(row[1])
        Cq.append(row[4])
        CqAvg.append(row[5])
        CqDev.append(row[6])
        sample.append(row[2])
        set_name.append(row[3])
    return target, content, CqAvg, CqDev, sample, set_name #function returns lists as output

def Ct_calculations(target, content, CqAvg, CqDev, sample, set_name):
    print "--------------Ct Calculation Inputs---------------"
    CqAvgAll = []
    CqDevAll = []
    sampleAll = []
    set_nameAll = []
    reference_flag = 0
    target_flag = 0
    if reference in target[0]:
        reference_flag = 1
    if one_target in target[0]:
        target_flag = 1
    print "Reference FLAG:", reference_flag
    print "Target FLAG:", target_flag
    for i in range(0, len(content), 2):
        if reference_flag == 1:
            print sample[i], target[i+1], "-", target[i], "Means:" , CqAvg[i+1], "-", CqAvg[i], "Standard Dev:", CqDev[i+1], "&", CqDev[i]
            sampleAll.append(sample[i])
            set_nameAll.append(set_name[i])
            CqAvgAll.append(float(CqAvg[i+1])-float(CqAvg[i])) #Cq calculation wither averages
            CqDevAll.append((float(CqDev[i+1])**2.0+float(CqDev[i])**2.0)**0.5) #Cq calculation with standard deviations
        if target_flag == 1:
            print sample[i], target[i], "-", target[i+1], "Means:" , CqAvg[i], "-", CqAvg[i+1], "Standard Dev:", CqDev[i], "&", CqDev[i+1]
            sampleAll.append(sample[i])
            set_nameAll.append(set_name[i])
            CqAvgAll.append(float(CqAvg[i])-float(CqAvg[i+1])) #Cq calculation wither averages
            CqDevAll.append((float(CqDev[i])**2.0+float(CqDev[i+1])**2.0)**0.5) #Cq calculation with standard deviations
    return sampleAll, set_nameAll, CqAvgAll, CqDevAll

def Ct_calculations_print(sampleAll, set_nameAll, CqAvgAll, CqDevAll, fileNum):
    print "------------Ct Calculation Results-----------------"
    for i in range(0, len(CqAvgAll)): #loop through and print all the contents of the lists
        print sampleAll[i], set_nameAll[i], CqAvgAll[i], CqDevAll[i]
    #Ct Calculations DataFrame
    print "\n********Ct Calc File*************\n"
    df_Ct = pd.DataFrame({'Sample IDs': sampleAll, 'Biological Sets': set_nameAll, 'Cq Averages': CqAvgAll, 'Cq Standard Deviations': CqDevAll}) #creates a dataframe
    print df_Ct #print all rows of dataframe
    print "\n********Ct Calc File SORTED*************\n"
    df_sort = df_Ct.sort_values('Biological Sets')
    print df_sort
    fileName = 'Blks1-20.AvUCP_Female_Cq_calculations_'+str(fileNum)+'.csv' #make a CSV file out of the dataframe unique to the file number (fileNum)
    df_Ct.to_csv(fileName) #print CSV contents to a file
    return df_sort

def Ct_calculations_merge(df_sort1, df_sort2, df_sort3, df_sort4, df_sort5):
    means = []
    set_names_merge = []
    #File1
    print "CALC MERGE DATAFRAME 1\n"
    print df_sort1
    set_nameAll_1 = df_sort1['Biological Sets'].tolist()
    CqAvgAll_1 = df_sort1['Cq Averages'].tolist()
    print "SET NAME ALL 1:\n", set_nameAll_1

    #File2
    print "CALC MERGE DATAFRAME 2\n"
    print df_sort2
    set_nameAll_2 = df_sort2['Biological Sets'].tolist()
    CqAvgAll_2 = df_sort2['Cq Averages'].tolist()
    print "SET NAME ALL 2:\n", set_nameAll_2

    #File3
    print "CALC MERGE DATAFRAME 3\n"
    print df_sort3
    set_nameAll_3 = df_sort3['Biological Sets'].tolist()
    CqAvgAll_3 = df_sort3['Cq Averages'].tolist()
    print "SET NAME ALL 3:\n", set_nameAll_3

    #File4
    print "CALC MERGE DATAFRAME 4\n"
    print df_sort4
    set_nameAll_4 = df_sort4['Biological Sets'].tolist()
    CqAvgAll_4 = df_sort4['Cq Averages'].tolist()
    print "SET NAME ALL 4:\n", set_nameAll_4

    #File5
    print "CALC MERGE DATAFRAME 5\n"
    print df_sort5
    set_nameAll_5 = df_sort5['Biological Sets'].tolist()
    CqAvgAll_5 = df_sort5['Cq Averages'].tolist()
    print "SET NAME ALL 5:\n", set_nameAll_5

    for i in range(0, len(set_nameAll_1)): #loop through all biological sets
        set_names_merge.append(set_nameAll_1[i])
        means.append((CqAvgAll_1[i]+CqAvgAll_2[i]+CqAvgAll_3[i]+CqAvgAll_4[i]+CqAvgAll_5[i]/5.0)) # +CqAvgAll_2[i]/2 dividing by 6 because there are 6 replicates - we are getting the average of the Ct calculations across all the files
        print set_nameAll_1[i], set_nameAll_2[i], set_nameAll_3[i], set_nameAll_4[i], set_nameAll_5[i]
    for j in range(0, len(set_names_merge)): #loop and print results
        print set_names_merge[j], means[j]
    return set_names_merge, means

def sem_calculation(set_names_merge, means):
    #Initialize the list of biological sets
    sets = ['FpreA', 'FpreB', 'FpreC', 'FpreD', 'FpreE', 'FpreF', 'FpostA', 'FpostB', 'FpostC', 'FpostD', 'FpostE', 'FpostF'] #MAKE SURE CALIBRATOR IS FIRST VALUE IN THE LIST
    FpreA = []
    FpreB = []
    FpreC = []
    FpreD = []
    FpreE = []
    FpreF = []
    FpostA = []
    FpostB = []
    Fpostc = []
    FpostD = []
    FpostE = []
    FpostF = []
    SEMs = []
    set_averages = []
    for i in range(0, len(set_names_merge)): #loop through and add to each set it's corresponding values
        if set_names_merge[i] == 'FpreA':
            FpreA.append(means[i])
        if set_names_merge[i] == 'FpreB':
            FpreB.append(means[i])
        if set_names_merge[i] == 'FpreC':
            FpreC.append(means[i])
        if set_names_merge[i] == 'FpreD':
            FpreD.append(means[i])
        if set_names_merge[i] == 'FpreE':
            FpreE.append(means[i])
        if set_names_merge[i] == 'FpreF':
            FpreF.append(means[i])
        if set_names_merge[i] == 'FpostA':
            FpostA.append(means[i])
        if set_names_merge[i] == 'FpostB':
            FpostB.append(means[i])
        if set_names_merge[i] == 'FpostC':
            FpostC.append(means[i])
        if set_names_merge[i] == 'FpostD':
            FpostD.append(means[i])
        if set_names_merge[i] == 'FpostE':
            FpostE.append(means[i])
        if set_names_merge[i] == 'FpostF':
            FpostF.append(means[i])
    #taking the mean of each set
    FpreA_avg = sum(FpreA)/len(FpreA)
    FpreB_avg = sum(FpreB)/len(FpreB)
    FpreC_avg = sum(FpreC)/len(FpreC)
    FpreD_avg = sum(FpreD)/len(FpreD)
    FpreE_avg = sum(FpreE)/len(FpreE)
    FpreF_avg = sum(FpreF)/len(FpreF)
    FpostA_avg = sum(FpostA)/len(FpostA)
    FpostB_avg = sum(FpostB)/len(FpostB)
    FpostC_avg = sum(FpostC)/len(FpostC)
    FpostD_avg = sum(FpostD)/len(FpostD)
    FpostE_avg = sum(FpostE)/len(FpostE)
    FpostF_avg = sum(FpostF)/len(FpostF)
    #averages combined into list
    set_averages.append(FpreA_avg)
    set_averages.append(FpreB_avg)
    set_averages.append(FpreC_avg)
    set_averages.append(FpreD_avg)
    set_averages.append(FpreE_avg)
    set_averages.append(FpreF_avg)
    set_averages.append(FpostA_avg)
    set_averages.append(FpostB_avg)
    set_averages.append(FpostC_avg)
    set_averages.append(FpostD_avg)
    set_averages.append(FpostE_avg)
    set_averages.append(FpostF_avg)
    #SEMs calculated and appended to list
    SEMs.append(stats.sem(FpreA))
    SEMs.append(stats.sem(FpreB))
    SEMs.append(stats.sem(FpreC))
    SEMs.append(stats.sem(FpreD))
    SEMs.append(stats.sem(FpreE))
    SEMs.append(stats.sem(FpreF))
    SEMs.append(stats.sem(FpostA))
    SEMs.append(stats.sem(FpostB))
    SEMs.append(stats.sem(FpostC))
    SEMs.append(stats.sem(FpostD))
    SEMs.append(stats.sem(FpostE))
    SEMs.append(stats.sem(FpostF))
    for j in range(0, len(set_averages)):
        print "Biological Set:", sets[j], "Average:", set_averages[j], "SEM:", SEMs[j]
    return sets, set_averages, SEMs

def delta_delta_Ct(sets, set_averages):
    calibrator_sample = set_averages[0] #corresponds to IS_0 since it's the first value in the list [index 0]
    dCt_calc = []
    dCt_sets = []
    for i in range(0, len(sets)):
        dCt_sets.append(sets[i])
        dCt_calc.append(set_averages[i]-calibrator_sample)
        print sets[i], set_averages[i]-calibrator_sample
    return dCt_sets, dCt_calc

def fold_change(dCt_sets, dCt_calc, SEMs):
    FC = []
    FC_range_1 = []
    FC_range_2 = []
    FC_sets = dCt_sets
    for i in range(0, len(dCt_calc)):
        print "Biological Set:", FC_sets[i], "Fold Change:", 2**(dCt_calc[i]*-1), "FC Range:", 2**((dCt_calc[i]+SEMs[i])*-1), "&", 2**((dCt_calc[i]-SEMs[i])*-1)
        FC.append(2**(dCt_calc[i]*-1)) #fold change calculation
        FC_range_1.append(2**((dCt_calc[i]+SEMs[i])*-1)) #fold change range 1
        FC_range_2.append(2**((dCt_calc[i]-SEMs[i])*-1)) #fold change range 2
    return FC_sets, FC, FC_range_1, FC_range_2

#Calling functions for PCR_file1
print "\n********PCR FILE 1*************\n"
csv1 = csv_init(PCR_file1)
target1, content1, CqAvg1, CqDev1, sample1, set_name1 = rows_init_store(csv1)
sampleAll1, set_nameAll1, CqAvgAll1, CqDevAll1 = Ct_calculations(target1, content1, CqAvg1, CqDev1, sample1, set_name1)
print_calc1 = Ct_calculations_print(sampleAll1, set_nameAll1, CqAvgAll1, CqDevAll1, 'file1')

#Calling functions for PCR_file2
print "\n********PCR FILE 2*************\n"
csv2 = csv_init(PCR_file2)
target2, content2, CqAvg2, CqDev2, sample2, set_name2 = rows_init_store(csv2)
sampleAll2, set_nameAll2, CqAvgAll2, CqDevAll2 = Ct_calculations(target2, content2, CqAvg2, CqDev2, sample2, set_name2)
print_calc2 = Ct_calculations_print(sampleAll2, set_nameAll2, CqAvgAll2, CqDevAll2, 'file2')

#Calling functions for PCR_file3
print "\n********PCR FILE 3*************\n"
csv3 = csv_init(PCR_file3)
target3, content3, CqAvg3, CqDev3, sample3, set_name3 = rows_init_store(csv3)
sampleAll3, set_nameAll3, CqAvgAll3, CqDevAll3 = Ct_calculations(target3, content3, CqAvg3, CqDev3, sample3, set_name3)
print_calc3 = Ct_calculations_print(sampleAll3, set_nameAll3, CqAvgAll3, CqDevAll3, 'file3')

#Calling functions for PCR_file4
print "\n********PCR FILE 4*************\n"
csv4 = csv_init(PCR_file4)
target4, content4, CqAvg4, CqDev4, sample4, set_name4 = rows_init_store(csv4)
sampleAll4, set_nameAll4, CqAvgAll4, CqDevAll4 = Ct_calculations(target4, content4, CqAvg4, CqDev4, sample4, set_name4)
print_calc4 = Ct_calculations_print(sampleAll4, set_nameAll4, CqAvgAll4, CqDevAll4, 'file4')

#Calling functions for PCR_file5
print "\n********PCR FILE 5*************\n"
csv5 = csv_init(PCR_file5)
target5, content5, CqAvg5, CqDev5, sample5, set_name5 = rows_init_store(csv5)
sampleAll5, set_nameAll5, CqAvgAll5, CqDevAll5 = Ct_calculations(target5, content5, CqAvg5, CqDev5, sample5, set_name5)
print_calc5 = Ct_calculations_print(sampleAll5, set_nameAll5, CqAvgAll5, CqDevAll5, 'file5')

#Calling function for merger of PCR Files
print "\n********PCR Calculations Merge*************\n"
set_names_merge_f, means_f = Ct_calculations_merge(print_calc1, print_calc2, print_calc3, print_calc4, print_calc5) # print_calc2 add print calc for file 4

#Calling function for calculations of mean and SEM of PCR Files
print "\n********PCR Averages & SEM*************\n"
bio_sets, avg, std_err = sem_calculation(set_names_merge_f, means_f)

#Calling function for calculations delta delta Ct
print "\n********PCR Delta Delta Ct*************\n"
delta_sets, delta_calc = delta_delta_Ct(bio_sets, avg)

#Calling function for calculating fold change
print "\n********PCR Fold Change*************\n"
fold_c_sets, fold_c_calc, fold_c_r1, fold_c_r2 = fold_change(delta_sets, delta_calc, std_err)

#Final DataFrame
print "\n********Fold Change File*************\n"
df_final = pd.DataFrame({'Biological Sets': fold_c_sets, 'Fold Change': fold_c_calc, 'Fold Change Range 1': fold_c_r1, 'Fold Change Range 2': fold_c_r2, 'SEM': std_err}) #creates a dataframe
print df_final.head() #print first 5 rows of dataframe

df_final.to_csv('Female_Chicken_FoldChanges.csv') #send dataframe to a file
