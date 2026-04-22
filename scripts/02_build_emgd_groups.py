from pathlib import Path
import numpy as np
import pandas as pd
import json
import argparse
from sklearn.cluster import KMeans
RANDOM_STATE = 2
EXTREMES_FRACTION = 0.10

#proteomics_train_df = pd.read_csv('', index_col=0) #Label-free proteomics data with samples on the columns and protein markers on the index
#genedependency_prob_df = pd.read_csv('', index_col=0) #Depmap probability data with samples on the columns and genes on the index
#n_train_size = len(proteomics_train_df.index)
#common_samples = proteomics_train_df.columns.intersection(genedependency_prob_df.columns)

def build_emgd_median(dependency_df, samples):
    '''Build EMGD groups by splitting samples into essential (ESS) and non-essential (NES) groups based on the median gene dependency probability for each gene.
    Args:
        dependency_df (pd.DataFrame): DataFrame containing gene dependency probabilities with samples on the columns and genes on the index.
        samples (list): List of sample names to consider for grouping. 
        Returns: A tuple of two dictionaries containing the NES and ESS groups for each gene.'''

    nes_group, ess_group = {},{}

    for gene in dependency_df.index:

        row = dependency_df.loc[gene, samples]
        median = np.median(row)

        nes_group[gene] = row[row < median].index.tolist()
        ess_group[gene] = row[row >= median].index.tolist()

    return nes_group, ess_group

#emgd_median_groups= build_emgd_median(genedependency_prob_df, common_samples)

def build_emgd_quartile(dependency_df, samples):
    '''Build EMGD groups by splitting samples into essential (ESS) and non-essential (NES) groups based on the first and third quartiles of dependency probabilities for each gene.
    Args:  
        dependency_df (pd.DataFrame): DataFrame containing gene dependency probabilities with samples on the columns and genes on the index.
        samples (list): List of sample names to consider for grouping.
        Returns: A tuple of two dictionaries containing the NES and ESS groups for each gene.'''
    
    nes_group, ess_group = {},{}

    for gene in dependency_df.index:
    
        row = dependency_df.loc[gene, samples]

        q1 = np.quantile(row, 0.25)
        q3 = np.quantile(row, 0.75)

        nes_group[gene] = row[row < q1].index.tolist()
        ess_group[gene] = row[row > q3].index.tolist()

    return nes_group, ess_group

#emgd_quartile_groups =  build_emgd_quartile(genedependency_prob_df, common_samples)

def build_emgd_high_low(dependency_df, samples):
    '''Build EMGD groups by splitting samples into essential (ESS) and non-essential (NES) groups based on fixed thresholds of dependency probabilities for each gene.
    Args:
        dependency_df (pd.DataFrame): DataFrame containing gene dependency probabilities with samples on the columns and genes on the index.
        samples (list): List of sample names to consider for grouping.      
        Returns: A tuple of two dictionaries containing the NES and ESS groups for each gene.'''

    nes_group, ess_group = {},{}
    for gene in dependency_df.index:
        row = dependency_df.loc[gene, samples]
        nes_group[gene] = row[row < 0.33].index.tolist()
        ess_group[gene] = row[row >= 0.66].index.tolist()
    return nes_group, ess_group

#emgd_high_low_groups = build_emgd_high_low(genedependency_prob_df, common_samples)


def build_emgd_extremes(dependency_df, fraction, samples):
    '''Build EMGD groups by selecting the top and bottom (10%) fractions of samples based on dependency probabilities for each gene.
    Args:
        dependency_df (pd.DataFrame): DataFrame containing gene dependency probabilities with samples on the columns and genes on the index.
        fraction (float): Fraction of samples to include in the essential (ESS) and non-essential (NES) groups for each gene.
        samples (list): List of sample names to consider for grouping.
        Returns: A tuple of two dictionaries containing the NES and ESS groups for each gene.'''
    
    nes_group, ess_group = {},{}

    for gene in dependency_df.index:
    
        ccl_dependency_zip = sorted(list(zip( samples,dependency_df.loc[gene,samples])), key = lambda x: x[1])
        names_sorted_by_dependency = [x[0] for x in ccl_dependency_zip]

        if len(names_sorted_by_dependency) < 3:
            raise ValueError("Not enough samples to build extremes groups for gene: {}".format(gene))  
        
        n_extremes = max(1, round(len(names_sorted_by_dependency)*fraction))  # Ensure at least one sample in each group

        nes_group[gene] = names_sorted_by_dependency[0:n_extremes]
        ess_group[gene] = names_sorted_by_dependency[-n_extremes:]

    return nes_group, ess_group

#emgd_ext_groups = build_emgd_extremes(genedependency_prob_df, fraction=0.10, common_samples)

def build_emgd_cluster(dependency_df, samples, random_state):
    '''Build EMGD clusters (n=2) using KMeans clustering on the dependency probabilities for each gene across samples.
    Args:
        dependency_df (pd.DataFrame): DataFrame containing gene dependency probabilities with samples on the columns and genes on the index.
        samples (list): List of sample names to consider for clustering.
        random_state (int): Random state for reproducibility of KMeans clustering.
        Returns: A tuple of two dictionaries containing the NES and ESS groups for each gene'''
    nes_group, ess_group = {},{}
    for gene in dependency_df.index:
        row = np.array(dependency_df.loc[gene, samples])

        kmeans = KMeans(n_clusters=2, n_init=10, random_state=random_state).fit(row.reshape(-1,1))
        labels = kmeans.labels_

        centroids = kmeans.cluster_centers_.flatten()
        sample_series = pd.Series(samples, index=np.arange(len(samples)))

        lower_cluster = int(np.argmin(centroids))
        upper_cluster = int(np.argmax(centroids))

        nes_group[gene] = sample_series[labels == lower_cluster].tolist()
        ess_group[gene] = sample_series[labels == upper_cluster].tolist()

    return nes_group, ess_group

#emgd_cluster_groups = build_emgd_cluster(genedependency_prob_df, common_samples)
def save_groups(output_path, nes_group, ess_group):
    payload = {
        "nes_group": nes_group,
        "ess_group": ess_group,
    }
    with open(output_path, "w") as f:
        json.dump(payload, f, indent=2)
        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dependency", required=True, help="Path to genedependency.csv file with samples on the columns and genes on the rows")
    parser.add_argument("--proteomics", required=True, help="Path to label-free proteomics.csv file with samples on the rows and protein names on the columns")
    parser.add_argument("--output-dir", required=True, help="Directory for output JSON files")
    args = parser.parse_args()


  
    dependency_df = pd.read_csv(args.dependency, index_col=0) #Gene dependency probability data with samples on the columns and genes on the rows
    proteomics_df = pd.read_csv(args.proteomics, index_col=0) #Label-free proteomics data with samples on the rows and protein names on the columns


    common_samples = proteomics_df.index.intersection(dependency_df.columns)

    if len(common_samples) == 0:
        raise ValueError("No common samples found between proteomics and dependency data.")

    dependency_df = dependency_df.loc[:,common_samples]
    proteomics_df = proteomics_df.loc[common_samples,:]

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    nes_group, ess_group = build_emgd_median(dependency_df, common_samples)
    save_groups(output_dir / "emgd_median_groups.json", nes_group, ess_group)

    nes_group, ess_group = build_emgd_quartile(dependency_df, common_samples)
    save_groups(output_dir / "emgd_quartile_groups.json", nes_group, ess_group)

    nes_group, ess_group = build_emgd_high_low(dependency_df, common_samples)
    save_groups(output_dir / "emgd_high_low_groups.json", nes_group, ess_group)

    nes_group, ess_group = build_emgd_extremes(dependency_df, fraction=EXTREMES_FRACTION, samples=common_samples)
    save_groups(output_dir / "emgd_extremes_groups.json", nes_group, ess_group)

    nes_group, ess_group = build_emgd_cluster(dependency_df, common_samples, RANDOM_STATE)
    save_groups(output_dir / "emgd_cluster_groups.json", nes_group, ess_group)

if __name__ == "__main__":
    main() 