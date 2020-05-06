''' 
Collection of functions for plotting
'''
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


def kernel_density_scatter_plt(x, y):
    ''' 
    Density scatter-plot using kernel density estimation (scatter plot colored by density of the data points)
    x and y are two vectors (e.g. pandas series or slices from df)
    '''
    from scipy.stats.kde import gaussian_kde

    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    
    # Plot and return figure
    fig, ax = plt.subplots()
    ax.scatter(x, y, c=z, s=50, edgecolor='')
    return fig


def density_scatter_plt(x,y):
    ''''Density scatter-plot using 2d histograms'''
    fig = plt.hist2d(x, y, (50, 50), cmap=plt.jet)
    plt.colorbar()
    return fig


def colored_scatter_plt(x, y, colorby):
    '''
    create scatter plot of y values vs x values, colored by a third variable (c).
    Returns a fig object
    '''
    fig = plt.scatter(x, y, c=colorby, s=50)
    plt.colorbar()
    return fig


def plot_hists_of_all_columns(df, hist_range):
    '''
    Create histograms of all columns in a pandas data frame (df).
    Currently built to plot up to 9 histograms.
    Histogram range (hist_range) is a list [min, max].
    '''
    fig = plt.figure()
    for subplt, col in enumerate(df.columns):
        if subplt>9:
            break
        sub_plt_str = int('33'+str(subplt))
        plt.subplot(sub_plt_str)
        plt.hist(df[col], bins=20, range = hist_range)
        plt.title(col)        
    return fig


def bar_graph(groups_names, y_vals, yerr_vals):
    '''
    Create simple bar graph. groups and y_vals should be lists.
    Returns a figure handle
    '''
    groups = ('Lowly expressed','Highly expressed')
    y_pos = np.arange(len(groups_names))+1
    fig = plt.bar(y_pos, y_vals, align='center', color='black', alpha=0.3, yerr=yerr_vals, ecolor='black')
    plt.xticks(y_pos, groups_names)
    return fig


def two_var_bar_graph(df1, df2, width = 0.35):
    '''
    Create a bar graph for two variables (two dfs) each containing the same number of columns.
    Data from each df will be plotted in different colored bars, grouped together by the column positions
    in the dfs.
    df_names is a list of names of the groups (e.g. ["df1", "df2"].
    Error bars represent 1 STE
    width = the width of the bars (default = 0.35)
    
    '''
    df1_means = df1.mean()
    df2_means = df2.mean()
    df1_STEs = np.std(df1)/np.sqrt(len(df1.index))
    df2_STEs = np.std(df2)/np.sqrt(len(df2.index))
    
    num_of_columns = len(df1_means)
    if num_of_columns != len(df2_means):
        print "Error - number of columns in the two dataframes must be equal"
    
    ind = np.arange(num_of_columns)  # the x locations for the groups   
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, df1_means, width, color='g', yerr=df1_STEs, ecolor="black")
    rects2 = ax.bar(ind + width, df2_means, width, color='b', yerr=df2_STEs, ecolor="black")
    
    # Text for labels, title and axes ticks
#     ax.set_ylabel('Mean expression (Z-score)')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(list(df1.columns))
#     ax.legend((rects1[0], rects2[0]), (df_names), loc='lower left')
    return fig



def norm_dist_fit_hist(data):
    '''
    Fit normal distribution to 1-D data array.
    Plot histogram and gaussian fit
    '''
    from scipy.stats import norm
    
    # Fit a normal distribution to the data:
    mu, std = norm.fit(data)

    # Plot the histogram.
    fig=plt.figure()
    plt.hist(data, bins=40, normed=True, alpha=0.5, color='b')
    
    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)
    text = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(text)
    return fig
#     plt.show()






    
def main():
    g=["sdf","dfd"]
    vals=[1234,532]
    yerr=[1,4]
    bar_graph(g, vals, yerr)
    
    # Generate random data
    x = np.random.normal(size=1000)
    y = x * 3 + np.random.normal(size=1000)
    
    fig1 = kernel_density_scatter_plt(x, y)
    
    fig2 = density_scatter_plt(x,y)

#     plt.show()

    
if __name__=="__main__":
    main()

