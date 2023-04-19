# Generate plots

The scripts in this directory can be used to generate the plots shown in Figures 4 and 5 of the article. Figure 4 (a) and (b) show the Parallel runtime breakdown of the receiver, sender (longest running) and total for GreeDIMM for the livejournal and wikipedia dataset respectively. Figure 5 (a) and (b) provides Parallel runtime breakdown for the receiver process of GreeDIMM, between its communicating thread and bucketing threads for the livejournal and wikipedia dataset respectively. 

The command to generate the plots in Figures 4 and 5 are:
`python3 plot_master.py`

This will result in the generation of the following:
* Figure4_livejournal_breakdown.PDF -- corresponding to the plot in Figure 4 (a) in the article.
* Figure4_wikipedia_breakdown.PDF -- corresponding to the plot in Figure 4 (b) in the article.
* Figure5_livejournal_receiver_breakdown.PDF -- corresponding to the plot in Figure 5 (a) in the article.
* Figure5_wikipedia_receiver_breakdown.PDF -- corresponding to the plot in Figure 5 (b) in the article.


