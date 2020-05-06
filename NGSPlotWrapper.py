'''
Run NGS Plot to create feature plots (wrapper for running through the command line, on local installation)
'''


import os
import subprocess

def run_ngsplot(in_bam, out_f, title = "", genome = "hg19", feature = "tss", flank_len = 3000, frag_len = 150):
    """
    Wrapper for running NGS-plot to create feature plots.
    Inputs:
        in_bam = string of path to input bam file for plotting (or two files separated by ":"
        out_f = prefix path for output (str)
        title = plot title string
        genome = refernce genome to use for feature extraction (e.g. "hg19")
        feature = tss / tes / genebody / exon / cgi / enhancer / dhs / bed

    See ngs.plot.r documentation for more details.
    """
    cmnd_dtr = "ngs.plot.r -G {G} -R {R} -C {in_f} -O {out_f} -T {T} -L {L} -FL {FL}".format(G=genome, R=feature, in_f=in_bam, out_f=out_f, T=title, L=flank_len, FL=frag_len)
    print ("Running...", cmnd_dtr)
    os.system(cmnd_dtr)
    # subprocess.Popen(cmnd_dtr) #***********
    print ("Complete")


def all_samples_feature_plots(bam_path, out_path, feature="tss"):
    """ Create simple feature plots for all 4 samples
    """
    samples = ["CMV1_HA-IE1_input_hg19_filtered.bam", "CMV2_HA-IE1_output_hg19_filtered.bam",
               "CMV3_HA-IE1_dCTD_input_hg19_filtered.bam", "CMV4_HA-IE1_dCTD_output_hg19_filtered.bam"]
    out_fs = ["IE1_input_hg19", "IE1_output_hg19", "IE1dCTD_input_hg19", "IE1dCTD_output_hg19"]
    titles = ["IE1_input", "IE1_output", "IE1dCTD_input", "IE1dCTD_output"]

    out_fs = [s + "_" + feature for s in out_fs]
    titles = [s + "_" + feature for s in titles]

    for i, s in enumerate(samples):
        in_bam = s
        out_f = out_fs[i]
        title = titles[i]

        in_bam = os.path.join(bam_path, in_bam)
        out_f = os.path.join(out_path, out_f)

        run_ngsplot(in_bam, out_f, title=title, genome="hg19", feature=feature, frag_len=150, flank_len=5000)

def input_normalized_feature_plt(samp1, samp2, bam_path, out_path, comparison_name, feature="tss"):
    out_f = os.path.join(out_path, comparison_name + "_" + feature)
    title = comparison_name + "_" + feature
    f1 = os.path.join(bam_path, samp1)
    f2 = os.path.join(bam_path, samp2)
    bam_str = f1 + ":" + f2
    run_ngsplot(in_bam = bam_str, out_f = out_f, title=title, genome="hg19", feature=feature, frag_len=150, flank_len=5000)


if __name__=="__main__":
    bam_path = "/Users/erankotler/Google\ Drive/workspace/CMV/Data/ChIPseq/bam"
    out_path = "/Users/erankotler/Dropbox/Limudim/CMV\ project/NGS_plot_output"

    # all_samples_feature_plots(bam_path, out_path, feature="tss")
    # all_samples_feature_plots(bam_path, out_path, feature="tes")
    # all_samples_feature_plots(bam_path, out_path, feature="genebody")
    # all_samples_feature_plots(bam_path, out_path, feature="exon")
    #
    # samp2 = "CMV1_HA-IE1_input_hg19_filtered.bam"
    # samp1 = "CMV2_HA-IE1_output_hg19_filtered.bam"
    # input_normalized_feature_plt(samp1, samp2, bam_path, out_path, comparison_name="IE1_vs_input", feature="tss")
    # input_normalized_feature_plt(samp1, samp2, bam_path, out_path, comparison_name="IE1_vs_input", feature="genebody")

    #
    # samp2 = "CMV3_HA-IE1_dCTD_input_hg19_filtered.bam"
    # samp1 = "CMV4_HA-IE1_dCTD_output_hg19_filtered.bam"
    # input_normalized_feature_plt(samp1, samp2, bam_path, out_path, comparison_name="IE1dCTD_vs_input", feature="tss")
    # input_normalized_feature_plt(samp1, samp2, bam_path, out_path, comparison_name="IE1dCTD_vs_input", feature="genebody")

    #
    # samp2 = "CMV1_HA-IE1_input_hg19_filtered.bam"
    # samp1 = "CMV3_HA-IE1_dCTD_input_hg19_filtered.bam"
    # input_normalized_feature_plt(samp1, samp2, bam_path, out_path, comparison_name="IE1_vs_IE1dCTD", feature="tss")
    # input_normalized_feature_plt(samp1, samp2, bam_path, out_path, comparison_name="IE1_vs_IE1dCTD", feature="genebody")


    configuration = "/Users/erankotler/Google\ Drive/workspace/CMV/ngsplot_config_files/IE1_dCTD_config.txt"
    # run_ngsplot(in_bam=configuration, out_f=os.path.join(out_path, "comparison_genebody"), title="Output - genebodies",
    #             genome="hg19", feature="genebody", frag_len=150, flank_len=5000)

    title = "Output_genebodies"
    out_f = os.path.join(out_path, "comparison_genebody")
    flank_len = 5000
    cmnd_dtr = "ngs.plot.r -G {G} -R {R} -C {in_f} -O {out_f} -T {T}".format(G="hg19", R="genebody", in_f=configuration, out_f=out_f, T=title)
    print ("Running...", cmnd_dtr)
    os.system(cmnd_dtr)


    # imr90_mnase_bam = "/Users/erankotler/Google\ Drive/workspace/CMV/Data/GSE21823_RAW/GSM543311_ResultCount_61DMKAAXX_s_7.map.bam"
    # # imr90_mnase_bam = "/Users/erankotler/Google\ Drive/workspace/CMV/Data/ChIPseq/bam/CMV1_HA-IE1_input_hg19_filtered.bam"
    #
    # run_ngsplot(in_bam = imr90_mnase_bam, out_f = os.path.join(out_path,"imr90_MNASe_genebody"), title="IMR90 MNAse seq", genome="hg19", feature="genebody", frag_len=150, flank_len=5000)
    #


    print ("Done")

