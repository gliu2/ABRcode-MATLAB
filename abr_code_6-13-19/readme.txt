Readme for abr_code MATLAB files
George Liu
6-14-19

Project: single trace ABR analysis with peak-peak amplitudes and template inner products

To use:
 - run ‘import_ABR_CSVfolder.m’ first to select folder containing ABR CSV data files, load ABR single trace data
- then run ‘jasaFigures_6-8-19.m’ to obtain figures, obtain summary stats (displayed and saved to variable for easy copying to Excel)
- run 'save_allfigsJASA.m' to select target folder, save figures. HARDCODE change dataset name. Run immediately to avoid altering figure order (which messes up figure naming).
- run 'jasaFigures_summary_6_12_19.m' to plot summary stats. Select Excel file containing summary stats, then select target folder to save plots.

Main files (2):
import_ABR_CSVfolder.m (10 dependencies)
jasaFigures_6-8-19.m  (1 dependency)
save_allfigsJASA.m
jasaFigures_summary_6_12_19.m (1 dependency)

Dependency files (12): (in directory ‘abr_utils’)
analyze_innerprod_ABR.m
import_ABRcsv.m
innerprod2pval.m
lags_xcov.m
plot_scatter2groups.m
PTDetect.m
same_xaxes.m
same_yaxes.m
vp2p_abr_sp_loc.m
vp2p_abr_sp.m
vp2p_abr.m
xgiveny.m