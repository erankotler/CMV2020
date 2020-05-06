
import pandas as pd
import numpy as np
from scipy import stats 
import os
# from Utils import shell_command
# from queue.qp import qp
# from addloglevels import sethandlers
from matplotlib import pyplot as plt
from graphics.Heatmap import Heatmap
import pickle
# from numpy.lib.scimath import log10
from matplotlib import pyplot as plt


def create_data_window(data_df, center_pos, bin_size, half_win_size):
    '''
    Extract a region (window) around a center point (e.g. TSS) in the data matrix.
    data_df should be a deeptools data matrix
    center_pos - window center position in bp from data start point
    half_win_size - number of bp to include on each side of the center_pos
    Output is a pandas dataframe with two columns of gene info and in-window data values
    '''
    win_start_pos = (center_pos - half_win_size ) / bin_size
    win_end_pos = (center_pos + half_win_size) / bin_size
    win_df = data_df.iloc[:,win_start_pos:win_end_pos+1] # extract data in window's range
    
    return win_df

def sortby_win_means(df, center_pos, bin_size=10, half_win_size=500):    
    '''Examine average across small window (usually around TSS). Sort genes by window average and save files:
    1) CSV file with gene names (sorted list)
    2) CSV file with genes and their binding scores (sorted)
    3) pickle file (.dat) of sorted means
    '''
    win_df = create_data_window(df, center_pos, bin_size, half_win_size) #extract window
    means = win_df.mean(axis=1) # calculate mean across window
    
    # Sort genes by average in window
    means.sort(axis=1, inplace=True, ascending=False) # sort means by descending order    
    sorted_indexes = means.index # list of gene names, sorted by window average
    pd.to_pickle(means, os.path.join(output_dir, out_base+"_win_means.dat"))
    means.to_csv(os.path.join(output_dir, out_base+"_sorted_win_means.csv"))
    sorted_means_df = means.to_frame("Mean_win_signal")
    gene_names = sorted_means_df.drop("Mean_win_signal", axis=1)
    gene_names.to_csv(os.path.join(output_dir, out_base + "_gene_list_sorted_by_win.csv"))
    print ("sorted genes (by window_mean signals) saved to pickle and to csv files")
    
    return means, win_df, sorted_indexes

def plt_win_avg_vs_genes(means, output_file_path):
    ''' Create dot plot of mean signal in examined window (usually TSS) per gene in data-set.
    (Input is a series)'''
    plt.plot(means, "ro")
    plt.title("Examined window mean across genes")
    plt.xlabel("Genes (sorted by window mean signal)")
    plt.ylabel("Mean signal in examined window")
    plt.savefig(output_file_path)
    plt.close()
    return

def get_expression_data():
    # get gene expression data (sorted low to high)    
    genes_sorted_by_exp_file = "/home/eranko//Data/mRNA/Expression/Human/MRC5/GSM579888_expression_noNA_sorted.txt"
    exp_df = pd.DataFrame.from_csv(genes_sorted_by_exp_file, sep='\t', index_col = 0, header = None)
    exp_df.columns=["Expression"]
 
    return exp_df

def plt_signal_vs_exp_level(sorted_means_df, exp_df, exp_output_file_path, logexp_output_file_path):
    ''''Plot signal in TSS region vs. expression level of the gene'''
#     exp_df = get_expression_data() # Get expression data:
#     sorted_means_df = sorted_means.to_frame("Mean_win_signal") # Convert window_means series to dataframe

    # plot signal vs expression:
    signal_vs_exp_df = pd.merge(exp_df, sorted_means_df, how='left', left_index=True, right_index=True, sort=False)
    signal_vs_exp_df = signal_vs_exp_df.dropna().sort_values(by="Expression")

    plt.plot(signal_vs_exp_df["Expression"], signal_vs_exp_df["Mean_win_signal"], "ro")
    plt.xlabel("Expression")
    plt.ylabel("Mean signal across window")
    r, p_val = stats.pearsonr(signal_vs_exp_df["Expression"], signal_vs_exp_df["Mean_win_signal"])
    plt.title("pearson r = {r}, p-val = {p_val}".format(r=r, p_val=p_val))
    plt.savefig(exp_output_file_path)
    plt.close()    
        
    # plot signal vs log10(expression):
    log_exp_df = np.log10(exp_df)
    signal_vs_logexp_df = pd.merge(log_exp_df, sorted_means_df, how='left', left_index=True, right_index=True, sort=False)
    signal_vs_logexp_df = signal_vs_logexp_df.dropna()

    plt.plot(signal_vs_logexp_df["Expression"], signal_vs_logexp_df["Mean_win_signal"], "ro")
    plt.xlabel("Expression (log10)")
    plt.ylabel("Mean signal across window")
    r, p_val = stats.pearsonr(signal_vs_logexp_df["Expression"], signal_vs_logexp_df["Mean_win_signal"])
    plt.title("pearson r = {r}, p-val = {p_val}".format(r=r, p_val=p_val))
    plt.savefig(logexp_output_file_path)
    plt.close()
    
    return        


#===================================================================================================
#===================================================================================================

def get_paths_and_parameters():
#   ### Directories:
    data_dir = "Z:/Runs/CMV/Deeptools_analyses/TSS_plots/"
    output_dir =  "Z:/Runs/CMV/Deeptools_analyses/TSS_plots/py_analyses/"

#   ### Matrix pickle file (created by /Runs/Erank/CMV/deepToolsMatrix_to_df_pickle.py):
#     matrix_pickle_file = "IE1_norm_dCTD_All_TSSs_RefPoint_matrix.dat"
    matrix_pickle_file = "IE1_output_All_TSSs_RefPoint_matrix.dat"
#     matrix_pickle_file = "IE1_dCTD_output_All_TSSs_RefPoint_matrix.dat"

#   ### Output file names base:
#     out_base = "HA_IE1_hg19_treat_tss_refGene_after5000"
#     out_base = "IE1_norm_dCTD_High2000"
#     out_base = "IE1_norm_dCTD_All_TSS_500bp"
#     out_base = "IE1_norm_dCTD_All_TSS_2000bp"
#     out_base = "IE1_output_All_TSS_2000bp"
#     out_base = "IE1_dCTD_output_All_TSS_2000bp"
    out_base = "IE1_output_All_TSS_500bp"
#     out_base = "IE1_dCTD_output_All_TSS_500bp"

    half_win_size = 500
#     half_win_size = 2000
    
    return data_dir, output_dir, matrix_pickle_file , out_base, half_win_size 


if __name__ == "__main__":
    # Get file paths and parameters
    data_dir, output_dir, matrix_pickle_file , out_base, half_win_size = get_paths_and_parameters()
    
    # Get data into pandas df:
#     data_file = os.path.join(data_dir,gz_f)
    df = pd.read_pickle(os.path.join(output_dir,matrix_pickle_file))
    
    # Extract per-gene mean binding in a window centered around the feature (TSS) 
    sorted_means, win_df, sorted_indexes = sortby_win_means(df, center_pos=5000, bin_size=10, half_win_size=half_win_size) # Sort genes by win average -> to CSV 
    sorted_means_df = sorted_means.to_frame("Mean_win_signal") # Convert window_means series to dataframe

    # Plot window mean per gene -> save figure file:
    plt_win_avg_vs_genes(sorted_means, os.path.join(output_dir, out_base+"_win_means.png"))
    
    
    print ("Done")
    
    
    #===================================================================================================
    #===================================================================================================
    
    ###################################################################################################3
##     ### gz_file (deeptools compute matrix output): 
# #     gz_f = "HA_IE1_hg19_treat_tss_refGene_after5000_matrix.gz"
# #     gz_f = "IE1_norm_dCTD_High2000_RefPoint_matrix.gz"
# #     gz_f = "IE1_norm_dCTD_All_TSSs_RefPoint_matrix.gz"
# #     gz_f = "IE1_norm_dCTD_All_TSSs_RefPoint_matrix.tab"    
#  
#### Expression analysis using MRC5 microarray data:  
#     exp_df = get_expression_data() # get expression data in MRC5
# 
#     # Plot window mean signal versus expression level in MRC5:    
#     plt_signal_vs_exp_level(sorted_means_df, exp_df, os.path.join(output_dir, out_base+"_sig_vs_expression.png"), \
#                             os.path.join(output_dir, out_base+"_sig_vs_logexpression.png"))
#     
#     # Look only at highly expressed genes:
#     plt.plot(range(len(exp_df)), exp_df) # plot expressions to set threshold
#     plt.show()
#     plt.close()
#
#     exp_threshold = 6 # set threshold to define high gene expression
#     fract_high = len(exp_df[exp_df["Expression"]>=6]) / float(len(exp_df)) # fract genes above thresh
#     high_exp_df = exp_df[exp_df["Expression"]>=6]
#     plt_signal_vs_exp_level(sorted_means_df, high_exp_df, os.path.join(output_dir, out_base+"_sig_vs_expression_high_exp_genes.png"), \
#                             os.path.join(output_dir, out_base+"_sig_vs_logexpression_high_exp_genes.png"))
#     
#     signal_vs_highexp_df = pd.merge(sorted_means_df,high_exp_df, how='left', left_index=True, right_index=True, sort=False)
#     signal_vs_highexp_df = signal_vs_highexp_df.dropna().sort_values(by="Mean_win_signal", ascending=False)    
#     signal_vs_highexp_df.to_csv(os.path.join(output_dir, out_base+"_sig_vs_expression_high_exp_genes.csv"))
#     
#     # sort dataframes:
#     sorted_win_df = win_df.ix[sorted_indexes] #sort window df
#     sorted_df = df.ix[sorted_means.sort_values().index] # sort entire data by window's average
#     
#     #create heatmap:
#     print("creating heatmap")
#     # create heatmap:
#     heatmap = Heatmap(sorted_df)
#     heatmap.savefig(os.path.join(output_dir, out_base+"_sorted_by_win.png"))

    



