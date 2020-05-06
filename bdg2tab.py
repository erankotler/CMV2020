'''
UPDATED VERSION STORED ON CLUSTER
'''

def bdg2tab(bdg_file):
    '''
    Expand bdg coverage file of the format chr_name<tab>start<tab>end<tab>score to a .tab file
    of the format chr_name<tab>position<tab>score (at a single bp resolution
    '''
    tab_file = bdg_file[:-3] + "tab"
    
    bdg_f = open(bdg_file, 'r')
    tab_f = open(tab_file,'w')
    
    for line in bdg_f:
        bdg_line = line.split("\t")
        chrom = bdg_line[0]
        start = int(bdg_line [1])+1
        end = int(bdg_line [2])+1
        new_range = range(start,end,1)
        
        for bp in new_range:
            line_str = "\t".join([chrom, str(bp), bdg_line[3]])
            tab_f.write(line_str)
        
    bdg_f.close()
    tab_f.close()
    



def tab2bdg(tab_file):
    '''
    Compress tab coverage file of the format  chr_name<tab>position<tab>score (at a single bp 
    resolution) to a .bdg file of the format chr_name<tab>start<tab>end<tab>score 
    '''
    import pandas as pd
    
    bdg_file = tab_file[:-3] + "bdg"
       
    df = pd.read_table(tab_file, sep = "\t", names = ["Chr", "pos", "score"]) # import data from tab file
    df['score_grp'] = (df.score.diff(1) != 0).astype('int').cumsum() # score_grp will increment by one whenever Value changes
    
    grp = df.groupby('score_grp')
    bdg_df = pd.DataFrame({'Chr': df.Chr[df.score_grp.unique()].replace(" ", ""),
                  'Start' : grp.pos.first()-1,
                  'End' : df.groupby('score_grp').pos.last(),
                  'Score': grp.score.first()})

    bdg_df = bdg_df[['Chr', 'Start', 'End', 'Score']]
       
    bdg_df.to_csv(bdg_file, '\t', header = False, index = False)

    
    
    


if __name__=="__main__":
#     bdg_file = "C:/Users/Eran/Desktop/CMV Temp stuff/tmp2.bdg"
#     bdg2tab(bdg_file)
    tab_file = "C:/Users/Eran/Desktop/CMV Temp stuff/tmp2.tab"
    tab2bdg(tab_file)
    print ("Done")
