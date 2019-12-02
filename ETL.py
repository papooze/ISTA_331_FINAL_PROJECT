'''Filename: ETL.py
This program creates visuals in order to explore our data for the final project.\

These visuals are: 
    
    A scatter plot
    A k-means clustering plot
    ???
    
    data epsilon, proteobacteria
    makes it uniform'''
    
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler


def read_csv():
    '''Creates a dataframe from the prokaryotes .csv, cleans the dataframe, then returns the dataframe.'''
    gc_df = pd.read_csv("data/prokaryotes.csv")
    gc_df.drop(gc_df.columns.difference(['Organism Groups','GC%', 'Scaffolds', 'Size(Mb)']), 1, inplace=True)
    for index, row in gc_df.iterrows():
        org = row['Organism Groups'].split(';')
        if org[-1] == 'delta/epsilon subdivisions':
            gc_df.loc[index, 'Organism Groups'] = 'Delta/Epsilonproteobacteria'
            continue
        gc_df.loc[index, 'Organism Groups'] = org[-1]
    gc_df = gc_df.rename(columns= {'Organism Groups': 'Type', 'GC%' : 'GC Content', 'Size(Mb)' : 'Genome Size'})
    gc_df.loc[gc_df['Type'].str.contains("unclassified"), 'Type'] = "Unclassified"
    gc_df.loc[gc_df['Type'].str.contains("environmental"), 'Type'] = "Unclassified"
    gc_df.loc[gc_df['Type'].str.contains("candidate phyla"), 'Type'] = "Unclassified"
    gc_df = gc_df[gc_df.Type != 'Unclassified']
    return gc_df


def euclidian_distance(ser1, ser2):
    '''
    Params:
        ser1 = first series
        ser2 = second series
    Returns:
        Euclidian distance as a float.
    Takes in two series objects and returns the Euclidian distance between the two.
    This function is instrumental in utilizing the k-means algorithm'''
    common_labels = set(ser1.index) & set(ser2.index)
    euc_dist = 0
    for label in common_labels:
        euc_dist += (ser1[label] - ser2[label])**2
    euc_dist = euc_dist ** .5
    if common_labels == set():
        return np.nan
    return euc_dist 

def scale(gc_df):
    '''
    Params:
        gc_df = The GC content dataframe.
    This function takes the GC content of the dataframe and scales it so they run between 0 and 1. '''
    scaler = MinMaxScaler(copy=False)
    scaler.fit_transform(gc_df)

def classify(centroidframe, fv):
    '''
    Params:
        centroidframe = Centroid dataframe.
        fv = Feature vector (a row from the climate frame)
    Returns:
        The label of the cluster whose centroid is closest to the feature vector.'''
    min_dist = euclidian_distance(centroidframe.loc[0], fv) #sets first centroid as minimum distance
    closest = 0 
    for cen in centroidframe.index:
        distance = euclidian_distance(centroidframe.loc[cen], fv) 
        if distance < min_dist:
            closest = cen
            min_dist = distance
    return closest

def get_initial_centroids(gc_df, k):
    '''
    Params:
        gc_df  = Scaled GC dataframe.
        k = Number of centroids
    Returns:
        Returns a dataframe containing starting values for centroids.'''
    centroid_frame = pd.DataFrame(data=[list(gc_df.iloc[i * int(len(gc_df) / k)]) for i in range(k)], columns= gc_df.columns)
    return centroid_frame

def get_labels(gc_df, centroidframe):
    '''
    Params:
        climateframe = Scaled climate dataframe.
        centroidframe = Centroid dataframe.
    Returns: 
        Returns a series that maps GC content of the bacteria to the labels of clusters that data belongs to.'''
    label_list = []
    for bac_type in gc_df.index:
        closest_cen = classify(centroidframe, gc_df.loc[bac_type])
        label_list.append(closest_cen)
    return pd.Series(data=label_list, index = gc_df.index)

def update_centroids(gc_df, centroidframe, labels):
    '''
    Params:
        gc_df = Scaled GC dataframe.
        centroidframe = Centroid dataframe.
        labels = Labels series.
    Goes through the labels series, and updates the centroid rows in centroidframe with averages of its rows in the labels.
    '''
    for col in centroidframe.columns:
        centroidframe[col].values[:] = 0
    for day in gc_df.index:
        centroid = labels[day]
        centroidframe.loc[centroid] += gc_df.loc[day]
    counts = labels.value_counts()
    for centroid in counts.index:
        centroidframe.loc[centroid] = centroidframe.loc[centroid]/counts[centroid]

def k_means(gc_df, k):
    '''
    Params:
        gc_df = Scaled GC dataframe.
        k = A value for how many centroids are desired.
    Returns:
        centroidframe = The final sorted centroid dataframe.
        labels = The labels of the data for the centroids.
    This function performs the k-means algorithm to adjust the centroids of the data 
    until their optimal spot to sort all of the data have been found.'''
    
    centroidframe = get_initial_centroids(gc_df, k)
    labels = get_labels(gc_df, centroidframe)
    while True:
        prev_labels = labels.copy()    
        update_centroids(gc_df, centroidframe, labels)
        labels = get_labels(gc_df, centroidframe)
        if labels.eq(prev_labels).all() == True: #If the label frames match-- ie the centroids have been found.
            break
    return centroidframe, labels

def distortion(df, labels, centroids):  
    '''This function measures how well clustering fits the data'''
    result = 0     
    for index in df.index:         
        label = labels[index]        
        result += euclidian_distance(df.loc[index], centroids.loc[label])**2    
    return result 
 
def kmeans_list(df, max_k):     
    list_of_kmeans = []     
    for k in range(1, max_k + 1):         
        centroids, labels = k_means(df, k)         
        sum_distances = distortion(df, labels, centroids)         
        list_of_kmeans.append({'centroids': centroids, 'labels': labels, 'k': k, 'distortion': sum_distances})     
    return list_of_kmeans     
    
def extract_distortion_dict(list_of_kmeans):     
    distortion_dict = {}     
    for d in list_of_kmeans:         
        distortion_dict[d['k']] = d['distortion']     
    return distortion_dict 

def make_distortion_graph(gc_df, k):
    '''Creates distortion graph for our data for use with the k-means algorithm.'''
    list_of_kmeans = kmeans_list(gc_df, k)
    dist_dict = extract_distortion_dict(list_of_kmeans)
    dist_ser = pd.Series(dist_dict)
    plt.figure(figsize = (15, 10))
    ax = dist_ser.plot(color = 'green')
    ax.set_xlabel('Centroid Number', size=24)
    ax.set_ylabel('Distortion', size=24)
    plt.figure()
    
def make_centroid_graph(gc_df, k):
    '''Creates scatter plot with centroids utilizing the k-means algorithm. Takes in a dataframe and k as an integer.'''
    cf, labels = k_means(gc_df, k)
    gc_df = pd.concat([gc_df,pd.DataFrame(labels)], axis=1)
    gc_df = gc_df.rename(columns= {0 : "Centroid"})
    plt.figure(figsize= (15, 10))
    ax = sns.scatterplot(x= 'GC Content', y= 'Genome Size', hue='Centroid', s=100, data= gc_df, palette= sns.husl_palette(k))
    plt.show()

def make_scatter(gc_df):
    '''Makes two scatter plots out of the data. One has a regressor.'''
    gc_df = gc_df.reset_index()
    plt.figure(figsize = (15, 10))
    ax = sns.scatterplot(x= 'GC Content', y= 'Genome Size', hue='Type', s=100, data= gc_df, palette= sns.husl_palette(len(gc_df)))
    ax.legend_.remove()
    plt.show()
    plt.figure(figsize = (15, 10))
    ax = sns.regplot(x= 'GC Content', y= 'Genome Size', data= gc_df, color='red')
    plt.show()
    
def main():
    gc_df = read_csv()
    gc_df = gc_df.groupby('Type').mean()
    scale(gc_df)
    make_centroid_graph(gc_df, 3)
    make_centroid_graph(gc_df, 5)
    make_centroid_graph(gc_df, 15)
    make_scatter(gc_df)
    make_distortion_graph(gc_df, 20)



    
    
main()