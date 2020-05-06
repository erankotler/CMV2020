import pandas as pd
import os 
from matplotlib import pyplot as plt
import numpy as np
# from graphics.Heatmap import *
from scipy import stats
from data_manipulations import standardize_df
from graphic_funcs import two_var_bar_graph, bar_graph, norm_dist_fit_hist
# from graphic_funcs import kernel_density_scatter_plt as kern_sct


def create_log_ratio_columns(df):   
    ''' Create df with TSS data and different ratios of expression data'''
    log_ratios_df=pd.DataFrame()
    
    # Add small constant instead of zeros for exp values:
    exp_uninfected = df["uninfected_ (mrna)"]
    exp_5hr_UV = df["5hr_UV_(mrna)"]
    exp_24hr = df["24hr_(mrna)"]
    exp_72hr = df["72hr_ (mrna)"]
    for d in [exp_uninfected, exp_5hr_UV, exp_24hr, exp_72hr]:
#         d.replace(to_replace=0, value=0.02, inplace=True) # alternaitve - replace zeros with 0.2
        if d.min()==0: # replace zeros (if exist) with next lowest value 
            lowest_non_zero = d.sort_values().unique()[1]
            d.replace(to_replace=0, value=lowest_non_zero, inplace=True)    

    log_ratios_df["TSS_signal"]= df['TSS_signal']
    log_ratios_df["24hr/uninfected"]=np.log2(exp_24hr/exp_uninfected)
    log_ratios_df["72hr/uninfected"]=np.log2(exp_72hr/exp_uninfected)
    log_ratios_df["24hr/5hrUV"]=np.log2(exp_24hr/exp_5hr_UV)
    log_ratios_df["72hr/5hrUV"]=np.log2(exp_72hr/exp_5hr_UV) 
    
    return log_ratios_df

def get_file_paths(data_dir, exp_file, TSS_means_file, output_analysis_dir, output_file_base):
    TSS_path = os.path.join(data_dir, TSS_means_file) 
    exp_path = os.path.join(data_dir, exp_file)
    output_base_path = os.path.join(output_analysis_dir, output_file_base)
    return TSS_path, exp_path, output_base_path

def get_exp_of_diff_bound_genes(df, num_of_top_bottom_genes):
    '''
    Extract expression levels of genes which are on top/bottom of TSS binding levels.
    Remove duplicate lines (same gene name) and keep the lines with the more extreme TSS binding value (high or low, depending on group) 
    '''
    high_bound_genes = df.sort_values('TSS_signal') # sort by ascending binding
    high_bound_genes['GeneName'] = high_bound_genes.index # add column of gene names
    high_bound_genes= high_bound_genes.drop_duplicates(subset='GeneName', keep='last') # drop duplicated genes (keep high)
    high_bound_genes= high_bound_genes.drop('GeneName', axis =1).tail(num_of_top_bottom_genes)
    
    low_bound_genes = df.sort_values('TSS_signal', ascending = False) # sort by descending binding
    low_bound_genes['GeneName'] = low_bound_genes.index
    low_bound_genes= low_bound_genes.drop_duplicates(subset='GeneName', keep='last') # drop duplicated genes (keep low)
    low_bound_genes= low_bound_genes.drop('GeneName', axis =1).tail(num_of_top_bottom_genes)
    
#     # Extraction of highly vs lowly bound genes without removing duplicates:
#     low_bound_genes = df.sort_values('TSS_signal').head(num_of_top_bottom_genes)
#     high_bound_genes = df.sort_values('TSS_signal').tail(num_of_top_bottom_genes)
    
    # Discard column of TSS binding from dfs:
    exp_of_low_bound_genes = low_bound_genes.iloc[:,1:]
    exp_of_high_bound_genes = high_bound_genes.iloc[:,1:]
    
    return exp_of_low_bound_genes, exp_of_high_bound_genes
 
def compare_exp_and_binding_signals(data_dir, exp_file, TSS_means_file, output_analysis_dir, output_file_base, num_of_top_bottom_genes=2000):
    '''
    Compare ChIP binding signals (for each gene, mean signal in a window around a feature e.g. TSS) to expression levels taken
    from RNA seq data (Noam Stern-Ginossar's paper).
    1) Load both datasets
    2) compare top highly expressed genes to lowly expressed genes in their mean binding levels.
    3) Compare highly bound genes to lowly bound genes in their expression levels.
    4) Save figures as output p.png file
    '''
    # Get paths for files
    TSS_path, exp_path, output_base_path = get_file_paths(data_dir, exp_file, TSS_means_file, output_analysis_dir, output_file_base)
    
    TSS_df = pd.read_csv(TSS_path, index_col =0, header=None)
    exp_df = pd.read_excel(exp_path, header=1, index_col = 0)
    RNA_df = exp_df.iloc[:, 17:24]
    
    df = TSS_df.join(RNA_df)
    df.columns.values[0]=u"TSS_signal"

    red_df = df.iloc[:,(0,5)] # chose column for specific comparison (col 4 = 24hr expression, col 5 = 72hr)
    d = red_df.dropna() # reduced df (without NaNs)
    
    #Log-transform expression levels (replacing zeros with next lowest value in given condition)
        
    d_log_exp = d.iloc[:,0].to_frame() # insert TSS values to new log-exp df
    
    for col in range(1, len(d.columns)): # replace zeros with lowest non-zero value
        d_exp = d.iloc[:,col]
        lowest_val = d_exp.sort_values().unique()[0] 
        second_lowest_val = d_exp.sort_values().unique()[1]
        lowest_non_zero = lowest_val if lowest_val !=0 else second_lowest_val
        d_exp.replace(to_replace=0, value=lowest_non_zero, inplace=True)    
    
        d_log_exp[d.columns[col]] = np.log2(d_exp) # log2 transform expression data   


    # Assess normality of distributions:
    for i,c in enumerate(d_log_exp.columns):
        fig = norm_dist_fit_hist(d_log_exp[c])
        if c == 'TSS_signal':
            plt_title = output_file_base + ' normal dist fit ' + c
        else:
            plt_title = c + ' normal dist fit'                
        plt.title(plt_title)
        plt.xlabel("Level")
        plt.ylabel("Count")
        plt.savefig(output_base_path + '_normality_' + str(i+1))
        plt.close()
    

    ### Compare TSS binding of genes that are highly lowly expressed: 
    '''
    Remove duplicte gene binding values, leaving the row with the highest binding Z score.
    This is under the assumption that binding, rather than depletion is what would affect gene
    expression changes
    '''
    # Removing duplicates:
    no_dups_log_ratio_df = d_log_exp.sort_values('TSS_signal')
    no_dups_log_ratio_df['GeneName'] = no_dups_log_ratio_df.index # add column of gene names
    no_dups_log_ratio_df = no_dups_log_ratio_df.drop_duplicates(subset='GeneName', keep='last') # drop duplicated genes (keep high)
    
    # Standardizing (z-scores) of TSS binding data:
    std_d = standardize_df(no_dups_log_ratio_df.iloc[:,:-1], column='TSS_signal') 
    
    # take top and bottom expressed genes and examine IE1 binding:
    low_exp = std_d.sort(std_d.columns[-1]).head(num_of_top_bottom_genes)
    high_exp = std_d.sort(std_d.columns[-1]).tail(num_of_top_bottom_genes)
    Mean_binding_of_low_exp = low_exp["TSS_signal"].mean()
    Mean_binding_of_high_exp = high_exp["TSS_signal"].mean()
    binding_STEs = (np.std(low_exp["TSS_signal"])/np.sqrt(num_of_top_bottom_genes), np.std(high_exp["TSS_signal"])/np.sqrt(num_of_top_bottom_genes))
    MannWhit_pval = stats.mannwhitneyu(low_exp["TSS_signal"], high_exp["TSS_signal"])[1]

    # plot comparison of expression levels:
    # Bar Graph (error bars=STE):
    plt.figure()
    groups = ('Lowly expressed','Highly expressed')
    y_pos = np.arange(len(groups))+1
    y_vals = [Mean_binding_of_low_exp, Mean_binding_of_high_exp]
    plt.bar(y_pos, y_vals, align='center', color='black', alpha=0.3, yerr=binding_STEs, ecolor='black')
    plt.ylim([-0.3,0.5])
    plt.xticks(y_pos, groups)
    plt.ylabel('IE1 binding (Z-score)')
    plt.title("MannWhit p = " + str(MannWhit_pval))
    plt.savefig(output_base_path + "_binding_" + str(num_of_top_bottom_genes) +"_HiLow_exp_genes_bar_graph.png")    
    plt.close()
    
    # Box plot (whiskers = 1.5 IQR)
    box_plt_data = [low_exp["TSS_signal"].values, high_exp["TSS_signal"].values] 
    plt.figure()
    plt.boxplot(box_plt_data, widths=0.5)
    plt.xticks([1,2],['Lowly expressed','Highly expressed'])
    plt.ylabel("TSS binding (Z-score)")
    plt.xlabel("Expression at 72hr post infection")
    plt.title(output_file_base + ' - ' + str(num_of_top_bottom_genes) + ' genes per group')
    plt.text(1.5,3.5,'MattWhitney p = ' + str(MannWhit_pval), horizontalalignment='center', verticalalignment='center')   
    plt.savefig(output_base_path + '_Binding_of_' + str(num_of_top_bottom_genes) +'_HiLow_exp_genes_boxplot')
    plt.close()
    
    
    ### Compare expression levels of genes that are highly/lowly bound by IE1 
    # Extract highly vs lowly bound genes, without getting duplicate genes:
    '''
    Note that the TSS binding data has many replications of the same gene name due to different 
    transcripts. This may cause some bias. I am removing duplicates, leaving the more extreme binding value
    (high values in highly bound gene group and lowest values in the lowly bound gene group
    '''
    # Remove duplicates and take top and bottom bound genes and examine expressions:
    exp_low_bound_genes, exp_high_bound_genes = get_exp_of_diff_bound_genes(d_log_exp, num_of_top_bottom_genes)
    
    # Calculate statistics for plots:
    Mean_exp_of_low_bound = exp_low_bound_genes.iloc[:,-1].mean()
    Mean_exp_of_high_bound = exp_high_bound_genes.iloc[:,-1].mean()
    exp_STEs = (np.std(exp_low_bound_genes.iloc[:,-1])/np.sqrt(num_of_top_bottom_genes), np.std(exp_high_bound_genes.iloc[:,-1])/np.sqrt(num_of_top_bottom_genes))
    MannWhit_pval = stats.mannwhitneyu(exp_low_bound_genes.iloc[:,-1], exp_high_bound_genes.iloc[:,-1])[1]
    
    # plot comparison of expression under High/Low binding levels:
    # Bar Graph:
    plt.figure()
    groups = ('Lowly bound','Highly bound')
    y_pos = np.arange(len(groups))+1
    y_vals = [Mean_exp_of_low_bound, Mean_exp_of_high_bound]
    plt.bar(y_pos, y_vals, align='center', color='black', alpha=0.3, yerr=exp_STEs, ecolor='black')
    plt.ylim([0,5])
    plt.xticks(y_pos, groups)
    plt.ylabel('log2(expression)')
    plt.title("MannWhit p=" + str(MannWhit_pval))
    plt.savefig(output_base_path + "_exp_of_" + str(num_of_top_bottom_genes) + "_highlow_bound_genes_bar_graph.png")    
    plt.close()
    
    # Box plot (whiskers = 1.5 IQR)
    plt.figure()
    plt.boxplot([exp_low_bound_genes.iloc[:,-1].values, 
                 exp_high_bound_genes.iloc[:,-1].values], widths=0.5, whis=1.5)
    plt.xticks([1,2],['Low bound', 'High bound'])
    plt.ylabel("log2 expression at 24hr post-infection")
    plt.xlabel("Binding to TSS")
    plt.title(output_file_base + ': ' + "MW p = " + str(MannWhit_pval) + ' (' + str(num_of_top_bottom_genes) + ' genes per group)')
    plt.savefig(output_base_path + '_Exp_' + str(num_of_top_bottom_genes) +'_highlow_bound_genes_boxplot')
    plt.close()
    
    
    # Do the same but looking at expression levels in all time points together: 
    df_no_na = df.dropna() 
    d_log_exp_all_timepoints = df_no_na.iloc[:,0].to_frame() # insert TSS values to new df (will hold all log-exp vals)
    
    for col in range(1, len(df_no_na.columns)): # replace zeros with lowest non-zero value
        d_exp = df_no_na.iloc[:,col]
        lowest_val = d_exp.sort_values().unique()[0] 
        second_lowest_val = d_exp.sort_values().unique()[1]
        lowest_non_zero = lowest_val if lowest_val !=0 else second_lowest_val
        d_exp.replace(to_replace=0, value=lowest_non_zero, inplace=True)    
    
        d_log_exp_all_timepoints[df_no_na.columns[col]] = np.log2(d_exp) # log2 transform expression data  

    # Standardizing (z-scores) of TSS binding data:
    d_std_TSS_rae_log_exp_all_timepoints = standardize_df(d_log_exp_all_timepoints, column='TSS_signal')    

    # Remove duplicates and take top and bottom bound genes and examine expressions:
    all_exp_low_bound_genes, all_exp_high_bound_genes = get_exp_of_diff_bound_genes(d_std_TSS_rae_log_exp_all_timepoints, num_of_top_bottom_genes)

    # Bar graph comparing expression levels of lowly bound vs highly bound genes across all RNA time points / conditions (error bars = STE): 
    fig = two_var_bar_graph(all_exp_low_bound_genes, all_exp_high_bound_genes, width=0.35)
    fig.set_size_inches(20, 10, 1)
    plt.legend(["Low binding", "High binding"], loc='best')
    plt.ylabel("log2(expression)")
    plt.ylim([0,5])
    plt.title(output_file_base + ' - ' + str(num_of_top_bottom_genes) + ' genes per group')  
    plt.savefig(output_base_path + '_all_TP_exp_of_' + str(num_of_top_bottom_genes) +'_highlow_bound_genes_raw')
    plt.close()


    # Standardizing both TSS binding values and log expression levels (z-scores):
    d_std_log_exp_all_timepoints = standardize_df(d_log_exp_all_timepoints, column='all')    

    # Remove duplicates and take top and bottom bound genes and examine expressions:
    all_exp_low_bound_genes, all_exp_high_bound_genes = get_exp_of_diff_bound_genes(d_std_log_exp_all_timepoints, num_of_top_bottom_genes)

    # Bar graph comparing expression levels of lowly bound vs highly bound genes across all RNA time points / conditions (error bars = STE): 
    fig = two_var_bar_graph(all_exp_low_bound_genes, all_exp_high_bound_genes, width=0.35)
    fig.set_size_inches(20, 10, 1)
    plt.legend(["Low binding", "High binding"], loc='best')
    plt.ylabel("log2(expression) - Z scores")
    plt.ylim([-0.3,0.15])
    plt.title(output_file_base + ' - ' + str(num_of_top_bottom_genes) + ' genes per group')
    plt.savefig(output_base_path + '_all_TP_exp_of_' + str(num_of_top_bottom_genes) +'_highlow_bound_genes_Zscores')
    plt.close()

    return

def compare_exp_fc_to_binding_signals(data_dir, exp_file, TSS_means_file, output_analysis_dir, output_file_base, num_of_top_bottom_genes=2000):
    '''
    Compare the binding signal of IE1 to genes which are up/down regulated following viral infection (according to RNA seq data)
    '''
    # Get file paths:
    TSS_path, exp_path, output_base_path = get_file_paths(data_dir, exp_file, TSS_means_file, output_analysis_dir, output_file_base)
    
    # Import data into united df:
    TSS_df = pd.read_csv(TSS_path, index_col =0, header=None)
    exp_df = pd.read_excel(exp_path, header=1, index_col = 0)
    RNA_df = exp_df.iloc[:, 17:24]
    df = TSS_df.join(RNA_df)
    df.columns.values[0]=u"TSS_signal"
    df = df.dropna()

    # Create df with log-ratios of 24hr & 72hr post-infection over uninfected & 5hr-UV infected cells (0s in expression are replaced with 0.02) 
    log_ratios_df = create_log_ratio_columns(df) 
    
    # Assess normality of distributions:
    for i,c in enumerate(log_ratios_df.columns):
        fig = norm_dist_fit_hist(log_ratios_df[c])
        if c == 'TSS_signal':
            plt_title = output_file_base + ' normal dist fit ' + c
        else:
            plt_title = c + ' normal dist fit'                
        plt.title(plt_title)
        plt.xlabel("Level")
        plt.ylabel("Count")
        plt.savefig(output_base_path + '_fc_normality_' + str(i+1))
        plt.close()
    
    
    # Data standardization:
    std_log_ratio_df = standardize_df(log_ratios_df, column='TSS_signal') # standardize only TSS signal to Z-scores, leave expression in raw scores
    std_all_log_ratio_df = standardize_df(log_ratios_df, column='all') # standardize both TSS signals and expression logFC to Z-scores


    ### Compare expression fold changes (up/down regulation) during viral infection in highly/lowly bound genes 
    # Extract highly vs lowly bound genes, without getting duplicate genes:
    '''
    Note that the TSS binding data has many replications of the same gene name due to different 
    transcripts. This may cause some bias. I decided to remove duplicates, but leave the more extreme binding value
    '''
    # Expression fold changes of highly/lowly bound genes (raw expression logFC scores), removing duplicate genes
    exp_of_low_bound_genes, exp_of_high_bound_genes = get_exp_of_diff_bound_genes(std_log_ratio_df, num_of_top_bottom_genes)
     
    # Bar graph comparing expression fold changes (log FC) of lowly bound vs highly bound genes (error bars = STE): 
    fig = two_var_bar_graph(exp_of_low_bound_genes, exp_of_high_bound_genes, width=0.35)
    plt.legend(["Low binding", "High binding"], loc='best')
    plt.ylabel("log2(expression fold change)")
    plt.ylim([0,0.9])
    plt.title(output_file_base + ' - ' + str(num_of_top_bottom_genes) + ' genes per group')
    plt.savefig(output_base_path + '_Exp_FC_of_' + str(num_of_top_bottom_genes) +'_highlow_bound_genes_raw')
    plt.close()
    
    # Box plot comparing expression fold changes at 72hr/5hrUV (log FC) of lowly bound vs highly bound genes. Whiskers = 1.5 IQR
    plt.figure()
    plt.boxplot([exp_of_low_bound_genes.iloc[:,-1].values, 
                 exp_of_high_bound_genes.iloc[:,-1].values], widths=0.5, whis=1.5)
    plt.xticks([1,2],['Low bound', 'High bound'])
    plt.ylabel("log2 expression fold change (72hr/5hrUV)")
    plt.xlabel("Binding to TSS")
    plt.title(output_file_base + ' - ' + str(num_of_top_bottom_genes) + ' genes per group')
    plt.savefig(output_base_path + '_Exp_FC_of_' + str(num_of_top_bottom_genes) +'_highlow_bound_genes_boxplot')
    plt.close()
    

    # Expression fold changes of highly/lowly bound genes (using standardized expression logFC (Z-scores))
    exp_Z_scores_of_low_bound_genes, exp_Z_scores_of_high_bound_genes = get_exp_of_diff_bound_genes(std_all_log_ratio_df, num_of_top_bottom_genes)
      
    # Bar graph comparing expression fold changes (log FC) of lowly bound vs highly bound genes (error bars = STE): 
    fig = two_var_bar_graph(exp_Z_scores_of_low_bound_genes, exp_Z_scores_of_high_bound_genes, width=0.35)
    plt.legend(["Low binding", "High binding"], loc='best')
    plt.ylabel("log2(expression fold change) - Z scores")
    plt.ylim([-0.3,0.15])
    plt.title(output_file_base + ' - ' + str(num_of_top_bottom_genes) + ' genes per group')
    plt.savefig(output_base_path + '_Exp_FC_of_' + str(num_of_top_bottom_genes) +'_highlow_bound_genes_Zscores')
    plt.close()


    ### Compare TSS binding of genes that are highly up/down-regulated during viral infection: 
    # Extract top up/down-regulated genes (expression changes the most between 72hr and 5hr UV):   
    '''
    Remove duplicte gene binding values, leaving the row with the highest binding Z score.
    This is under the assumption that binding, rather than depletion is what would affect gene
    expression changes
    '''
    no_dups_log_ratio_df = log_ratios_df.sort_values('TSS_signal')
    no_dups_log_ratio_df['GeneName'] = no_dups_log_ratio_df.index # add column of gene names
    no_dups_log_ratio_df = no_dups_log_ratio_df.drop_duplicates(subset='GeneName', keep='last') # drop duplicated genes (keep high)
    std_no_dups_log_ratio_df = standardize_df(no_dups_log_ratio_df.iloc[:,:-1], column='TSS_signal')
    
    
    downregulated_genes = std_no_dups_log_ratio_df.sort_values('72hr/5hrUV').head(num_of_top_bottom_genes)
    upregulated_genes = std_no_dups_log_ratio_df.sort_values('72hr/5hrUV').tail(num_of_top_bottom_genes)   

    # Extract only TSS binding data from these dfs:
    TSS_binding_of_downregulated_genes = downregulated_genes.iloc[:,0]
    TSS_binding_of_upregulated_genes = upregulated_genes.iloc[:,0]

    # Calculate statistics:
    MannWhit_pval = stats.mannwhitneyu(TSS_binding_of_downregulated_genes, TSS_binding_of_upregulated_genes)[1]
    
    # Bar graph (error bars = STE):
    fig = bar_graph(groups_names=['Top downregulated', 'Top upregulated'], 
                    y_vals=[TSS_binding_of_downregulated_genes.mean(),TSS_binding_of_upregulated_genes.mean()], 
                    yerr_vals=[np.std(TSS_binding_of_downregulated_genes)/np.sqrt(500),np.std(TSS_binding_of_upregulated_genes)/np.sqrt(500)])
    plt.ylabel("TSS binding (Z-score)")
    plt.xlabel("Expression change (72hr/5hrUV)")
    plt.ylim([-0.1,0.5])
    plt.title(output_file_base + ' - ' + str(num_of_top_bottom_genes) + ' genes per group')
    plt.text(1.5,0.4,'MattWhitney p = ' + str(MannWhit_pval), horizontalalignment='center', verticalalignment='center')
    plt.savefig(output_base_path + '_Binding_of_' + str(num_of_top_bottom_genes) +'_UpDown_regulated_genes')
    plt.close()

    # Box plot (whiskers = 1.5 IQR)
    box_plt_data = [TSS_binding_of_downregulated_genes.values, TSS_binding_of_upregulated_genes.values]
    plt.figure()
    plt.boxplot(box_plt_data, widths=0.5)
    plt.xticks([1,2],['Top downregulated', 'Top upregulated'])
    plt.ylabel("TSS binding (Z-score)")
    plt.xlabel("Expression change (72hr/5hrUV)")
    plt.title(output_file_base + ' - ' + str(num_of_top_bottom_genes) + ' genes per group')
    plt.text(1.5,3.5,'MattWhitney p = ' + str(MannWhit_pval), horizontalalignment='center', verticalalignment='center')
    plt.savefig(output_base_path + '_Binding_of_' + str(num_of_top_bottom_genes) +'_UpDown_regulated_genes_boxplot')
    plt.close()

    print("Done") 

    return log_ratios_df, output_base_path

def intersect_lists(log_ratios_df, output_base_path, num_of_top_bottom_genes):
# data_dir, exp_file, TSS_means_file, output_analysis_dir, output_file_base, num_of_top_bottom_genes):
    '''
    Create list of genes that intersect between highly/lowly bound genes and up/down regulated during viral infection.
    The following intersections are returned (as dfs with binding Z-scores and raw log-FC expression changes):
        Highly bound AND Down-regulated
        Highly bound AND Upregulated 
        Lowly bound AND Upregulated 
        Lowly bound AND Down-regulated
    '''
    # standardize only TSS binding values to Z-scores, leave expression in raw scores (log-FC)
    std_log_ratio_df = standardize_df(log_ratios_df, column='TSS_signal') 
    
    # Organize data - drop dups (keep high TSS binding value) and sort by binding score:
    sorted_by_TSS = std_log_ratio_df.sort_values('TSS_signal') # sort by ascending binding
    sorted_by_TSS['GeneName'] = sorted_by_TSS.index # add column of gene names
    sorted_by_TSS= sorted_by_TSS.drop_duplicates(subset='GeneName', keep='last') # drop duplicated genes (keep highest binding value)
    
    # get list of highly bound TSSs:
    high_bound_genes= sorted_by_TSS.drop('GeneName', axis =1).tail(num_of_top_bottom_genes)
    
    # get list of lowly bound TSSs:
    low_bound_genes= sorted_by_TSS.drop('GeneName', axis =1).head(num_of_top_bottom_genes) 
    
    # Organize expression data - drop dups and sort by expression FC:
    no_dups_log_ratio_df = std_log_ratio_df.sort_values('TSS_signal') # sort by TSS bindinf for duplkicate removal
    no_dups_log_ratio_df['GeneName'] = no_dups_log_ratio_df.index # add column of gene names
    no_dups_log_ratio_df = no_dups_log_ratio_df.drop_duplicates(subset='GeneName', keep='last') # drop duplicated genes (keep high)
    no_dups_log_ratio_df = no_dups_log_ratio_df.iloc[:,:-1]
    
    # get list of up-regulated genes 72Hr/5hrUV:
    upregulated_genes = no_dups_log_ratio_df.sort_values('72hr/5hrUV').tail(num_of_top_bottom_genes)   
    
    # get list of up-regulated genes 72Hr/5hrUV:
    downregulated_genes = no_dups_log_ratio_df.sort_values('72hr/5hrUV').head(num_of_top_bottom_genes)

    # create 4 intersections
    # Highly bound AND Down-regulated
    high_bound_AND_downreg = pd.concat([high_bound_genes, downregulated_genes], axis=1, join='inner').iloc[:,:5]

    # Highly bound AND Up-regulated:
    high_bound_AND_upreg = pd.concat([high_bound_genes, upregulated_genes], axis=1, join='inner').iloc[:,:5]
 
    # Lowly bound AND Up-regulated:
    low_bound_AND_upreg = pd.concat([low_bound_genes, upregulated_genes], axis=1, join='inner').iloc[:,:5]
    
    # Lowly bound AND Down-regulated    
    low_bound_AND_downreg = pd.concat([low_bound_genes, downregulated_genes], axis=1, join='inner').iloc[:,:5]

    # save intersection gene lists with their corresponding values 
    high_bound_AND_downreg.to_csv(output_base_path + '_high_bound_AND_downreg_'+ str(num_of_top_bottom_genes) + '.csv')
    high_bound_AND_upreg.to_csv(output_base_path + '_high_bound_AND_upreg'+ str(num_of_top_bottom_genes) + '.csv')
    low_bound_AND_upreg.to_csv(output_base_path + '_low_bound_AND_upreg'+ str(num_of_top_bottom_genes) + '.csv')
    low_bound_AND_downreg.to_csv(output_base_path + '_low_bound_AND_downreg'+ str(num_of_top_bottom_genes) + '.csv')

    # save csv with all gene data (without dupplicates)
    no_dups_log_ratio_df.to_csv(output_base_path + '_all_genes_no_dups.csv')
    return

    


def main():
    data_dir = "C:/Users/Eran/Dropbox/Limudim/CMV project/TSS_Signal_vs_RNA_expression/Data"
    output_analysis_dir = "C:/Users/Eran/Dropbox/Limudim/CMV project/TSS_Signal_vs_RNA_expression/Analysis"
    exp_file = "S1_table_RPKM_hg19.xlsx" # data from Noam Stern-Ginnosar's paper
    
    TSS_means_file = "IE1_norm_dCTD_All_TSS_500bp_sorted_win_means.csv"
#     TSS_means_file = "IE1_output_All_TSS_500bp_sorted_win_means.csv"
#     TSS_means_file = "IE1_dCTD_output_All_TSS_500bp_sorted_win_means.csv"
    
    output_file_base = "IE1_norm_IE1dCTD"
#     output_file_base = "IE1_output"
#     output_file_base = "IE1_dCTD_output"

#     num_of_top_bottom_genes = 2000 # number of top/bottom-expressed or top/bottom-bound genes to examine
#     for i in [300, 500, 1000, 1500, 2000]:
    for i in [500]:



#   # Comparing ChIP signal of high/low expressed genes and comparing expression of high/low ChIP bound genes
#         compare_exp_and_binding_signals(data_dir, exp_file, TSS_means_file, output_analysis_dir, output_file_base, num_of_top_bottom_genes=i)
    
#   # Compare expression changes at late vs early infection (or before infection) to IE1 binding
        log_ratios_df, output_base_path = compare_exp_fc_to_binding_signals(data_dir, exp_file, TSS_means_file, output_analysis_dir, output_file_base, num_of_top_bottom_genes=i)
        
#   # Create lists for intersections of highly/lowly exoressed vs. highly/lowly bound genes:
        intersect_lists(log_ratios_df, output_base_path, num_of_top_bottom_genes=i)

    print ("Figures saved to: " + output_analysis_dir)
    print("Done Main") 

#####################################################################################################################
#####################################################################################################################

   

#     # Normalize by 5hr UV expression
#     RNA_df_norm_uv = RNA_df.div(RNA_df["5hr_UV_(mrna)"], axis='index')
#     RNA_df_norm_uv.drop("5hr_UV_(mrna)", axis=1, inplace=True)
#     RNA_df_norm_uv = np.log2(RNA_df_norm_uv)
#     
#     clr_range=[-6,6] # define heatmap value range
#     heatmap = Heatmap(RNA_df_norm_uv, color_range=clr_range) # create heatmap
#     
#     top = RNA_df_norm_uv.ix[:3000,:]
#     df = TSS_df.join(top)
#     df.columns.values[0]=u"TSS_signal"
#     df = df.dropna()
#     Top = df.groupby(df.index).first()
#     heatmap = Heatmap(Top, color_range=clr_range) # create heatmap
#     heatmap.cluster(axis='rows', row_method='weighted') # hierarchial clustering of rows
#     heatmap.show()
#     
# #     heatmap.cluster(row_method='weighted', column_method='weighted') # hierarchial clustering
#     heatmap.cluster(axis='rows', row_method='weighted') # hierarchial clustering of rows
# 
# #     heatmap.savefig(data_dir + "_tmp.png")
# #     heatmap.show()
#     


    
#     df.corr(method="spearman")
    




if __name__=='__main__':
    main()
    

    