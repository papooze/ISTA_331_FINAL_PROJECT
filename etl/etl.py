"""
Filename: etl.py
Author: Michael Burman
Group collaberation with: Ash Rau
"""
import pandas as pd
import numpy as np

def csv_reader(filenames):
    '''This function takes in a list of filenames and returns a list of dataframes.'''
    dataframes = []
    for file in filenames:
        dataframe = pd.read_csv(file)
        dataframes.append(dataframe)
    return dataframes

def get_descriptive_stats(dataframes):
    '''This function takes in a list of dataframes and returns a list of descriptive stats of each csv. (unique values, column types)'''
    for df in dataframes:
        print(df.describe())
        print(df.nunique())

def main():
    dataframes = csv_reader(['/Users/micha/Documents/ISTA 331 DOCS/Homework/Final Project/data/eukaryotes.csv','/Users/micha/Documents/ISTA 331 DOCS/Homework/Final Project/data/prokaryotes.csv','/Users/micha/Documents/ISTA 331 DOCS/Homework/Final Project/data/viruses.csv'])
    get_descriptive_stats(dataframes)
    
main()