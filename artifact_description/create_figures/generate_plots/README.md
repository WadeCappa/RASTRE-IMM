# Generate plots

The scripts in this directory can be used to generate the plots shown in Figures 4 and 5 of the article. Figure 4 (a) and (b) shows the parallel runtime breakdown of the receiver, longest running sender, and the total runtime for GreeDIMM. Figure 4 (a) shows the livejournal breakdown and figure 4 (b) shows the wikipedia breakdown. Figure 5 (a) and (b) provides the parallel runtime breakdown for the receiver process, with seperate bars for the communicating thread and bucketing threads. Figure 5 (a) shows the livejournal receiver breakdown and figure 5 (b) shows the wikipedia receiver breakdown. 

The command to generate the plots in Figures 4 and 5 are: `python3 plot_master.py`

This will result in the generation of the following:
* Figure4_livejournal_breakdown.PDF -- corresponding to the plot in Figure 4 (a) in the article.
* Figure4_wikipedia_breakdown.PDF -- corresponding to the plot in Figure 4 (b) in the article.
* Figure5_livejournal_receiver_breakdown.PDF -- corresponding to the plot in Figure 5 (a) in the article.
* Figure5_wikipedia_receiver_breakdown.PDF -- corresponding to the plot in Figure 5 (b) in the article.


