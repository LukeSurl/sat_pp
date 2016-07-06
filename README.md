# sat_pp
Code to post-process and analyse satellite observations of trace gases and equivalent model output.

This utility first loads processed satellite observation and associated model data. Presently this data needs to be that in the final stage of processing (processed by ulitilties not *yet* contained in this package).
This utility launches a menu, through which various different ways of analysing the data are accessed. The ultility can create maps showing spatial data , scatter plots, time series, as well as other statisical analyses.

The language used is python, and various python packages are required. The utility runs on JASMIN if a virtual environment is used with the correct packages "pip install"ed.

I (Luke Surl) created this repositry on 2016-07-06. At initiation, this utility encompasses some, but not all, of the mostly independent final-stage visualisation tools I have previously developed locally. Work on incorporating the remaining visualisation tools will continue.

I will also incorporate into this repositry the tools developed for earlier-stage processing and analysis.
