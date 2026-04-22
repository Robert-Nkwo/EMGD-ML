 #%% main loop for identifying emperical markers of gene dependence
from pathlib import Path
from collections import Counter
import argparse
import json
import multiprocessing
import time
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import RepeatedKFold


RANDOM_STATE = 2    

warnings.filterwarnings(
    "ignore",
    message="Precision loss occurred in moment calculation",
    category=RuntimeWarning
)

def load_emgd_groups(group_json_path):
    with open(group_json_path, 'r') as f:
        group_dict = json.load(f)
    return group_dict["nes_group"], group_dict["ess_group"]

def student_ttest(group1_df, group2_df, proteins):
    group1_mean = np.mean(group1_df[proteins], axis=0)
    group2_mean = np.mean(group2_df[proteins], axis=0)
    mean_delta = group1_mean - group2_mean

    pvals = stats.ttest_ind(
        group1_df[proteins],
        group2_df[proteins],
        axis=0,
        nan_policy="omit",
        equal_var=False,
    ).pvalue

    pval_dict = dict(zip(proteins, pvals))
    delta_dict = dict(zip(proteins, mean_delta))

    return pval_dict, delta_dict

def run_ttests(nes_comp_df, ess_comp_df):
    proteins = nes_comp_df.columns.tolist()
    return student_ttest(nes_comp_df, ess_comp_df, proteins)



def generate_emgd(proteomics_df,
                   gene_list, 
                   nes_group, 
                   ess_group, 
                   output_csv,
                   min_group_size = 9,
                   n_splits = 3,
                   n_repeats = 20,
                   n_processes = 2,
                   min_markers = 100,
                   max_markers = 300
):
   
    """Generate EMGD markers for a list of genes by comparing proteomics data between precomputed NES and ESS groups.
    Args:
        proteomics_df (pd.DataFrame): DataFrame containing proteomics data with samples on the index and protein markers on the columns.
        gene_list (list): List of gene names for which to generate EMGD markers.
        nes_group (dict): Dictionary mapping each gene to a list of sample names in the non-essential (NES) group.
        ess_group (dict): Dictionary mapping each gene to a list of sample names in the essential (ESS) group.
        output_csv (str): Path to save the output CSV file containing EMGD markers for each gene.
        min_group_size (int): Minimum number of samples required in each group to perform comparisons (default: 9).
        n_splits (int): Number of splits for repeated K-fold cross-validation (default: 3).
        n_repeats (int): Number of repeats for repeated K-fold cross-validation (default: 100).
        n_processes (int): Number of parallel processes to use for t-tests (default: 2).
        min_markers (int): Minimum number of ESS or NES EMGD markers for each gene (default: 100).
        max_markers (int): Maximum number of ESS or NES EMGD markers for each gene (default: 300).
    Returns:
        pd.DataFrame: DataFrame containing EMGD markers for each gene, with columns for NES markers, ESS markers, and counts of each."""
    
    results_df = pd.DataFrame(
    columns=["marker_nes", "marker_ess", "n_nes", "n_ess"],
    index=gene_list,
    ) 
    total_time = []
    print(f'Processing a total of {len(gene_list)} genes.')
    with multiprocessing.Pool(processes=n_processes) as pool:
        for i, gene in enumerate(gene_list):
            print(f'Processing gene {i+1}/{len(gene_list)}: {gene}')
            if gene not in nes_group or gene not in ess_group:
                print(f'Warning: {gene} is missing in the sample dependency dict, skipping.')
                continue

            nes_samples = [s for s in nes_group[gene] if s in proteomics_df.index]
            ess_samples = [s for s in ess_group[gene] if s in proteomics_df.index]
                    
            if len(nes_samples) < min_group_size or len(ess_samples) < min_group_size:
                print(
                    f"Warning: Comparisons for {gene} has fewer than {min_group_size} samples in one or both groups; skipping"
                )
                continue

            start = time.time()

            X_nes = proteomics_df.loc[nes_samples,:]
            X_ess = proteomics_df.loc[ess_samples,:]

            kf = RepeatedKFold(
                n_splits=n_splits,
                n_repeats=n_repeats,
                random_state=RANDOM_STATE,
            )

            nes_kfold_bank = [split[1] for split in kf.split(X_nes)]
            ess_kfold_bank = [split[1] for split in kf.split(X_ess)]

            fold_args = []
            for nes_idx, ess_idx in zip(nes_kfold_bank, ess_kfold_bank):
                nes_comp_df = X_nes.iloc[nes_idx, :]
                ess_comp_df = X_ess.iloc[ess_idx, :]
                fold_args.append((nes_comp_df, ess_comp_df))

        
            fold_results = pool.starmap(run_ttests, fold_args) #list of tuples of pval_dict and delta_dict for each fold - see student_ttest function

            pval_bank = [x[0] for x in fold_results]
            delta_bank = [x[1] for x in fold_results]

            dir_down = []
            dir_up = []

            for delta_dict in delta_bank:
                for protein, delta in delta_dict.items():
                    if delta <= 0:
                        dir_down.append(protein)
                    else:
                        dir_up.append(protein)

            # Count how many times each protein is downregulated or upregulated across comparisons
            # downregulated = higher expression in ESS samples, upregulated = higher expression in NES samples
            tally_dir_down = Counter(dir_down)
            tally_dir_up = Counter(dir_up)

            pval_dict = {}
            for protein in proteomics_df.columns:
                combined_p = stats.combine_pvalues(
                    [fold_pvals[protein] for fold_pvals in pval_bank]
                )[1]

                pval_dict[protein] = 2 if np.isnan(combined_p) else combined_p

            top_proteins = sorted(pval_dict, key=lambda x: float(pval_dict[x]))

            nes_markers = []
            ess_markers = [] 
            n_randomised = len(pval_bank)

            # NES selection
            for protein in top_proteins[:max_markers]:
                if tally_dir_up[protein] >= round(n_randomised / 2) and pval_dict[protein] < 0.05:
                    nes_markers.append(
                        (protein, float(tally_dir_up[protein]) / n_randomised, float(pval_dict[protein]))
                    )

            # ESS selection
            for protein in top_proteins[:max_markers]:
                if tally_dir_down[protein] >= round(n_randomised / 2) and pval_dict[protein] < 0.05:
                    ess_markers.append(
                        (protein, float(tally_dir_down[protein]) / n_randomised, float(pval_dict[protein]))
                    )

            if len(nes_markers) < min_markers:
                to_add = min_markers - len(nes_markers)
                n_added = 0
                for protein in top_proteins[max_markers:]:
                    if (
                        tally_dir_up[protein] >= round(n_randomised / 2)
                        and pval_dict[protein] < 0.05
                        and n_added < to_add
                    ):
                        nes_markers.append(
                            (protein, float(tally_dir_up[protein]) / n_randomised, float(pval_dict[protein]))
                        )
                        n_added += 1

            if len(ess_markers) < min_markers:
                to_add = min_markers - len(ess_markers)
                n_added = 0
                for protein in top_proteins[max_markers:]:
                    if (
                        tally_dir_down[protein] >= round(n_randomised / 2)
                        and pval_dict[protein] < 0.05
                        and n_added < to_add
                    ):
                        ess_markers.append(
                            (protein, float(tally_dir_down[protein]) / n_randomised, float(pval_dict[protein]))
                        )
                        n_added += 1

            results_df.loc[gene, "marker_nes"] = nes_markers
            results_df.loc[gene, "marker_ess"] = ess_markers
            results_df.loc[gene, "n_nes"] = len(nes_markers)
            results_df.loc[gene, "n_ess"] = len(ess_markers)

            elapsed = time.time() - start
            total_time.append(elapsed)
            print(f"{gene} has {len(nes_markers)} NES samples")
            print(f"{gene} has {len(ess_markers)} ESS samples")
            print(f"{gene} completed in {elapsed:.2f} seconds")

            results_df.to_csv(output_csv)

    return results_df
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dependency", required=True, help="Path to dependency CSV")
    parser.add_argument("--proteomics", required=True, help="Path to proteomics CSV")
    parser.add_argument("--groups", required=True, help="Path to EMGD group JSON")
    parser.add_argument("--output", required=True, help="Path to output marker CSV")
    parser.add_argument("--min-group-size", type=int, default=9)
    args = parser.parse_args()

    dependency_df = pd.read_csv(args.dependency, index_col=0)
    proteomics_df = pd.read_csv(args.proteomics, index_col=0)

    # dependency_df: genes x samples
    # proteomics_df: samples x proteins
    shared_samples = proteomics_df.index.intersection(dependency_df.columns)

    if len(shared_samples) == 0:
        raise ValueError("No common samples found between dependency and proteomics data.")

    dependency_df = dependency_df.loc[:, shared_samples]
    proteomics_df = proteomics_df.loc[shared_samples, :]

    nes_group, ess_group = load_emgd_groups(args.groups)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    generate_emgd(
        proteomics_df=proteomics_df,
        gene_list=dependency_df.index.tolist(),
        nes_group=nes_group,
        ess_group=ess_group,
        output_csv=output_path,
        min_group_size=args.min_group_size,
    )

if __name__ == "__main__":
    main() 