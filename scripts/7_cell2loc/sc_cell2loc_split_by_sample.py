
# AnnData tutorial
# https://anndata-tutorials.readthedocs.io/en/latest/getting-started.html

import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import cell2location

from cell2location.utils.filtering import filter_genes
from cell2location.models import Cell2location, RegressionModel

print(ad.__version__)

import torch
torch.cuda.is_available()
torch.cuda.current_device()
torch.cuda.get_device_name(torch.cuda.current_device())

import jax
print(jax.devices())


# results_folder = "/mnt/data0/projects/ASAP/ST/st_sc_deconvolution/results_V2_08032023"
results_folder = "/mnt/data0/projects/ASAP/ST/st_sc_deconvolution/results_V4_02202024"
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f"{results_folder}/reference_signatures"
run_name = f"{results_folder}/cell2location_map"
epoch_plot_name = f"{results_folder}/epoch_convergence_plots"

# sc_dir = "/mnt/data0/projects/ASAP/ST/data/sc_for_deconvolution/sc_for_deconvolution_V2_08032023"
sc_dir = "/mnt/data0/projects/ASAP/ST/data/sc_for_deconvolution/sc_for_deconvolution_V4_02202024"
# adata_vis_all = sc.read_h5ad("/mnt/data0/projects/ASAP/ST/data/comb_seurat_smoothed_final_01102023.h5ad")
adata_vis_all = sc.read_h5ad("/mnt/data0/projects/ASAP/ST/data/comb_seurat_post_smoothed_label_s5_diet_11302023.h5ad")
# adata_vis_all is converted from Seurat, so X contains the lognormalized "data" slot -
# substitute the counts in raw.X for X
adata_vis_all = adata_vis_all.raw.to_adata()

# drop_cols = [x for x in adata_vis_all.obs.columns if x.startswith("Spatial_snn_res")] + \
#             [x for x in adata_vis_all.obs.columns if x.startswith("smoothed_label_")] + \
#             ["seurat_clusters", "is_layer1", "final_smoothed_label_minusL1", "layer_label"]
drop_cols = [x for x in adata_vis_all.obs.columns if x.startswith("Spatial_snn_res")] + \
            ["seurat_clusters", "smoothed_label_s1", "s2_subclust", "s2_subclust_score",
             "s2_nbr_score", "smoothed_label_s2", "smoothed_label_s3", "smoothed_label_s4"]
adata_vis_all.obs = adata_vis_all.obs.drop(columns=drop_cols)
# add shortened sample_id, i.e. BN0347_rep5 --> BN0347
adata_vis_all.obs["sample_id"] = adata_vis_all.obs["orig.ident"].str.split("_rep").str[0]

# get list of spatial samples with matching single cell file
paired_samples = []

for single_sample in sorted(set(adata_vis_all.obs["sample_id"])):
    # sc_meta_file = os.path.join(sc_dir, f"sc_count_{single_sample}_meta_02062023.csv")
    sc_meta_file = os.path.join(sc_dir, f"{single_sample}_Metadata.csv")
    if not os.path.exists(sc_meta_file):
        continue
    print(f"Running {single_sample}")
    

    # filter for a single sample_id
    adata_vis = adata_vis_all[adata_vis_all.obs["sample_id"]==single_sample, :]



    # https://github.com/BayraktarLab/cell2location/blob/master/docs/notebooks/cell2location_short_demo.ipynb
    # Calculate QC metrics
    adata_vis.X = adata_vis.X.toarray()
    sc.pp.calculate_qc_metrics(adata_vis, inplace=True)
    adata_vis.X = csr_matrix(adata_vis.X)
    adata_vis.var['mt'] = [gene.startswith('MT-') for gene in adata_vis.var.index]
    adata_vis.obs['mt_frac'] = adata_vis[:, adata_vis.var['mt'].tolist()].X.sum(1).A.squeeze()/adata_vis.obs['total_counts']

    # add sample name to obs names
    # adata_vis.obs["new_spot_ids"] = [str(i) for i in adata_vis.obs["new_spot_ids"]]
    # adata_vis.obs_names = adata_vis.obs["new_spot_ids"]
    # adata_vis.obs.index.name = "spot_id"

    # mitochondria-encoded (MT) genes should be removed for spatial mapping
    adata_vis.obsm['mt'] = adata_vis[:, adata_vis.var['mt'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['mt'].values]


    # remove genes with very low expression
    # sc.pp.filter_genes(adata_vis, min_counts=10, inplace=True)
    # parameters from https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
    # On a few samples, filter_genes triggers the following warning:
    # /home/jy1008/anaconda3/envs/cell2loc/lib/python3.9/site-packages/pandas/core/arraylike.py:396: RuntimeWarning: divide by zero encountered in log10
    selected = filter_genes(
        adata_vis, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
    )
    adata_vis = adata_vis[:, selected].copy()




    # import the sc data
    # https://chanzuckerberg.github.io/scRNA-python-workshop/preprocessing/00-tabula-muris.html
    # https://docs.scvi-tools.org/en/stable/tutorials/notebooks/cell2location_lymph_node_spatial_tutorial.html
    # counts_file = sc_meta_file.replace("_meta", "")
    counts_file = sc_meta_file.replace("Metadata", "Counts")
    print(sc_meta_file)
    print(counts_file)
    
    # put quotes around brackets (9,10]
    # sed -i 's/(/"(/g' BN*_Metadata.csv
    # sed -i 's/\[/"\[/g' BN*_Metadata.csv
    # sed -i 's/)/)"/g' BN*_Metadata.csv
    # sed -i 's/\]/\]"/g' BN*_Metadata.csv
    count_df = pd.read_csv(counts_file, index_col=0)
    count_df = count_df.transpose()
    meta_df = pd.read_csv(sc_meta_file, index_col=0, quotechar='"')
    # meta_df = meta_df.drop(columns=["orig.ident", "RNA_snn_res.1.2", "seurat_clusters", "class", "dblscore"])
    meta_df = meta_df.drop(columns=["orig.ident", "RNA_snn_res.1.5", "seurat_clusters", "class", "dblscore"])

    sc_adata = sc.AnnData(X = count_df, obs = meta_df)

    # Calculate QC metrics
    # sc_adata.X = sc_adata.X.toarray() # didn't use csr_matrix earlier, so already an array
    sc.pp.calculate_qc_metrics(sc_adata, inplace=True)
    sc_adata.X = csr_matrix(sc_adata.X)
    sc_adata.var['mt'] = [gene.startswith('MT-') for gene in sc_adata.var.index]
    sc_adata.obs['mt_frac'] = sc_adata[:, sc_adata.var['mt'].tolist()].X.sum(1).A.squeeze()/sc_adata.obs['total_counts']

    # mitochondria-encoded (MT) genes should be removed for spatial mapping
    sc_adata.obsm['mt'] = sc_adata[:, sc_adata.var['mt'].values].X.toarray()
    sc_adata = sc_adata[:, ~sc_adata.var['mt'].values]

    # parameters from https://docs.scvi-tools.org/en/stable/tutorials/notebooks/cell2location_lymph_node_spatial_tutorial.html
    # This results in 16111 genes kept
    selected = filter_genes(
        sc_adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
    )
    sc_adata = sc_adata[:, selected].copy()

    # in SC data, remove very sparse categories
    label_freq_thresh = 20
    # cell_type_column = "CellSubtypes"
    cell_type_column = "Complex_Assignment"
    # adata_vis.obs["final_smoothed_label_plusL1"].value_counts()
    cell_counts = sc_adata.obs[cell_type_column].value_counts()
    cell_types = cell_counts[cell_counts > label_freq_thresh].index.tolist()
    sc_adata = sc_adata[sc_adata.obs[cell_type_column].isin(cell_types), :]


    # prepare anndata for the regression model
    RegressionModel.setup_anndata(
        adata=sc_adata,
        # 10X reaction / sample / batch
        batch_key=None,
        # cell type, covariate used for constructing signatures
        labels_key=cell_type_column,
        # multiplicative technical effects (platform, 3' vs 5', donor effect)
        categorical_covariate_keys=None,
        continuous_covariate_keys=None
    )
    mod = RegressionModel(sc_adata)
    mod.view_anndata_setup()
    mod.train(
        max_epochs=2000,
        batch_size=None,
        train_size=1,
        lr=0.002,
    )

    for skip_epoch in [0, 20, 100]:
        plt.figure()
        mod.plot_history(skip_epoch)
        epoch_plot = "cell2location_sc_nbregression_skipepoch%s_%s_02232023.pdf" % (skip_epoch, single_sample)
        print(epoch_plot)
        # plt.savefig(os.path.join(results_folder, epoch_plot))
        plt.savefig(os.path.join(epoch_plot_name, epoch_plot))
        plt.close()


    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    sc_adata = mod.export_posterior(
        sc_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
    )

    # Save model
    mod.save(dir_path=f"{ref_run_name}", prefix=f"sc_nbregression_{single_sample}_", overwrite=True)

    # Save anndata object with results
    adata_file = f"{ref_run_name}/sc_adata_{single_sample}.h5ad"
    sc_adata.write(adata_file)
    adata_file

    # mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", sc_adata)
    # adata_file = f"{ref_run_name}/sc.h5ad"
    # sc_adata = sc.read_h5ad(adata_file)

    # export estimated expression in each cluster
    if "means_per_cluster_mu_fg" in sc_adata.varm.keys():
        print("means_per_cluster_mu_fg in sc_adata.varm.keys()")
        inf_aver = sc_adata.varm["means_per_cluster_mu_fg"][
            [f"means_per_cluster_mu_fg_{i}" for i in sc_adata.uns["mod"]["factor_names"]]
        ].copy()
    else:
        inf_aver = sc_adata.var[
            [f"means_per_cluster_mu_fg_{i}" for i in sc_adata.uns["mod"]["factor_names"]]
        ].copy()
    
    inf_aver.columns = sc_adata.uns["mod"]["factor_names"]
    inf_aver.iloc[0:5, 0:5]

    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()


    Cell2location.setup_anndata(adata=adata_vis, batch_key="sample_id")

    mod = Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=10,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection (using default here):
        detection_alpha=20,
    )
    mod.view_anndata_setup()

    mod.train(
        max_epochs=30000,
        # train using full data (batch_size=None)
        # Batch size of 1024 is used in Supp of Kleshchevnikov et al., 2020
        batch_size=1024,
        # use all data points in training because
        # we need to estimate cell abundance at all locations
        train_size=1,
    )

    # plot ELBO loss history during training, removing first 100 epochs from the plot

    for skip_epoch in [0, 100, 1000]:
        plt.figure()
        mod.plot_history(skip_epoch)
        epoch_plot = "cell2location_sc_cell2location_skipepoch%s_%s_02232023.pdf" % (skip_epoch, single_sample)
        print(epoch_plot)
        plt.legend(labels=["full data training"]);
        plt.savefig(os.path.join(epoch_plot_name, epoch_plot))
        plt.close()


    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs, "use_gpu": True},
    )

    # Save model
    mod.save(dir_path=f"{run_name}", prefix=f"st_cell2loc_map_{single_sample}_", overwrite=True)

    # mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

    # Save anndata object with results
    adata_file = f"{run_name}/st_adata_{single_sample}.h5ad"
    adata_vis.write(adata_file)
    adata_file


    # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
    # to adata.obs with nice names for plotting
    # (this is reportedly more accurate than using mean)
    adata_vis.obs[adata_vis.uns["mod"]["factor_names"]] = adata_vis.obsm[
        "q05_cell_abundance_w_sf"
    ]

    # write just the obs to file, then join with metadata in R
    adata_vis.obs.to_csv(f"{run_name}/st_cellprop_table_{single_sample}.csv")
