'''
Basically until Michael gets back, I want to figure out how to group the damn data here
'''

# Suggestions: 
# Seaborn for generating colors versus picking 44 colors by hand
# Make a classifier that can classify new bacteria based off of GC content and perhaps genome size?
import pandas as pd, numpy as np
import matplotlib as plt
import seaborn as sns

def df_clean_up(df):
    df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")
    df.columns = df.columns.str.replace("%", "")
    df.columns = df.columns.str.replace("#", "")

def group_em(df):
    bacteria_dict = {}
    for ind in df.index:
        group = df.loc[ind, "organism_groups"]
        group = group[9:]
        gc = df.loc[ind, "gc"]
        pos = group.find(";")
        family = group[:pos]
        sub_group = group[pos+1:]
        if family == "Proteobacteria":
                if sub_group == "delta/epsilon subdivisions":
                    sub_group = "Detla/Epsilonproteobacteria"
        if family not in bacteria_dict:
            bacteria_dict[family] = {sub_group: [gc]}
        elif sub_group not in bacteria_dict[family]:
            bacteria_dict[family][sub_group] = [gc]
        else:
            bacteria_dict[family][sub_group].append(gc)
    return bacteria_dict

def clean_bacteria_dict(bacteria_dict):
    pop_em = []
    for family in bacteria_dict.keys():
        for sub_family in bacteria_dict[family].keys():
            length = len(bacteria_dict[family][sub_family])
            if length == 1:
                pop_em.append([family, sub_family])
            new_val = sum(bacteria_dict[family][sub_family])/length
            bacteria_dict[family][sub_family] = new_val
    for fails in pop_em:
        if fails[0] in bacteria_dict.keys():
            bacteria_dict[fails[0]].pop(fails[1])
            if len(bacteria_dict[fails[0]]) == 0:
                bacteria_dict.pop(fails[0])

def main():
    df = pd.read_csv("data/prokaryotes.csv")
    df_clean_up(df)
    bacteria_dict = group_em(df)
    bacteria_dict.pop("unclassified Bacteria")
    clean_bacteria_dict(bacteria_dict)
    sum_this = 0
    df = pd.DataFrame(bacteria_dict)
    for family in bacteria_dict.keys():
        for sub_fam in bacteria_dict[family].keys():
            sum_this += 1
    df.to_csv()

if __name__ == '__main__':
    main()