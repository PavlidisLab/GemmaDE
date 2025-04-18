# Gemma DE
Gemma DE is an R Shiny app that allows users to search the Gemma database by gene (set) to uncover biological conditions which are commonly and robustly associated with the differential expression of that gene (set). This is done by ranking experiments in which one or several of the query genes have been found as differentially expressed and reporting the enrichment of the experiments' associated ontology terms against those contained in the entire corpus.

# Design of Gemma DE
At its core, Gemma DE is a two-file Shiny app. This means that the user interface is defined in [`ui.R`](main/ui.R) and the server logic is defined in [`server.R`](main/server.R). There are also three additional files for organization: [`process.R`](main/process.R) which holds the algorithms that drive the condition enrichment, [`renderTools.R`](main/renderTools.R) which contains scripts necessary for visualizations, and [`gemmaAPI.R`](main/gemmaAPI.R), which contains logic to interface with Gemma.

Users can pass options to Gemma DE, which are packed into a named list that gets passed around as appropriate, and as long as a user remains connected to Gemma DE via their web browser, they also maintain a server-side copy of their output in their session data. These data are not insignificantly small, so care should be taken to not allow too many simultaneous connections.

# Running Gemma DE
After installing the necessary dependencies (note that some dependencies are from my personal GitHub repository, not CRAN), sourcing [`start.R`](start.R) will sequentially load the required packages, data and then start the web server. Assuming the app is run on a dedicated server, users can connect on their local machine by SSH tunneling. This is as simple as executing the following command in a command prompt: `ssh -L <LOCAL PORT>:localhost:<GEMMA DE PORT> <USERNAME>@<GEMMA DE SERVER> -p <GEMMA DE PORT>`, replacing items in triangle braces with the correct options. Once completed, you can connect to the user interface by opening a web browser and navigating to `localhost:<LOCAL PORT>`.

# The algorithm
Gemma DE works by tabulating a number of per-experiment differential expression properties (magnitude of fold change, number of genes differentially expressed, etc.) and using these to calculate a ranking for experiment in the corpus. We calculate this ranking for 500,000 randomly selected genes, so that we obtain a distribution of rankings per experiment (roughly 5 billion total rankings). For later searches, we calculate so-called "fundamental relatedness scores" as the z-score of an experiment's ranking on this empirical distribution. These fundamental relatedness scores are then used as weights in an enrichment module which associates all experiments with their inferred biological contrasts.

The idea behind this is that if we can first predict how strongly and specifically the differential expression of a gene is to an experiment, we can then use this information across a pool of experiments to predict how this reflects on the biological conditions that they study.

## Regenerating data
Since Gemma DE works on a data freeze from Gemma, updates to the freeze will require updates to the calculated "prior" that this algorithm uses. Scripts for this purpose are located within the [`generate`](generate) directory. To regenerate the data, users should first delete old data:
1. The old `CACHE.BACKGROUND2.rds` and `DATA.HOLDER.rds`
2. The contents of the `updated_nulls2` directory

This ensures that new data will be created. Following this, complete regeneration (unfortunately) requires user interaction at multiple steps.
1. Run [start.R](start.R) to recreate tag caches and differential expression data
2. Once completed, run [bootstrap.R](generate/bootstrap.R) to create experiment-wise expected rankings (see the first paragraph of the above section)

This packages all necessary data into `CACHE.BACKGROUND.rds` (expanded ontology terms) and `DATA.HOLDER.rds` (experimental data). It also creates files for normalizing experiment scores in `updated_nulls2`. The experiment data takes up a lot of space (around 30 GB at the time of writing). So when `load.R` is called, the data is cold read into RAM momentarily and gene-major, read-only file-backed matrices of DE data are then used to replace them in RAM. These have the advantage of not requiring the entire matrix to be maintained in RAM and reducing the memory footprint of Gemma DE significantly.

# Using Gemma DE

## Through the User Interface
Genes can be input (either by entry into the textbox or by file upload) as NCBI or Ensembl gene IDs, official symbol, description, known aliases or by GO group. An optional differential expression signature can also be inputted (in the same order as the genes) to attempt to find biological conditions whose fold changes are correlated with the query. There are more options which can be expanded by clicking on the corresponding text to restrict which categories are reported, modify filtering options and etc. Of note in this category is the scoring method. For searches not involving a query differential expression signature, we recommend leaving this on the default. Searches involving a query signature should use the modified vector space model (M-VSM) or correlation (we describe these in our publication). Clicking search will begin an enrichment and update the URL, which can be saved for future reference.

Many visualization options beyond just the results table are available. These are accessible either via. the tabs (located above the results table) or by clicking on the "Visualize" button. Clicking the "Visualize" button will send a request to the Gemma database to fetch gene expression information (and may take some time to retrieve, depending on how much data is being fetched), while all tab options recalculate views based on already available data (and are thus likely faster to compute).

Data can also be downloaded by clicking the "Download" button (for the results table), or the download icon (in the "Visualize" view).

## Programmatically
We are currently working on a REST API for accessing Gemma DE.

# Evaluation of Gemma DE
The performance of Gemma DE was evaluated in large part by running enrichments on either entirely or partially simulated data. Other metrics were obtained on real data by searching gene sets with known expression differences in some biological processes (ie. biological sex).

# Future work
This can be used for a number of related analyses and some improvements to the codebase will help support them.

1. The reverse operation of this condition enrichment is relatively straightforward (query with a condition-comparison and get back a ranked set of genes). This might be interesting, and it is theoretically possible to construct a network of condition-comparisons in terms of how closely coincident their gene sets are.
2. It would be nice to have a way to threshold the output with some statistical certainty. This may be done by some sort of permutation test (ie. shuffling the contrast labels repeatedly and re-scoring).
