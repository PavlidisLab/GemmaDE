# Gemma DE
Gemma DE is an R Shiny app that allows users to search the Gemma database by gene (set) to uncover biological conditions which are commonly and robustly associated with the differential expression of that gene (set). This is done by ranking experiments in which one or several of the query genes have been found as differentially expressed and reporting the enrichment of the experiments' associated ontology terms against those contained in the entire corpus.

# Design of Gemma DE
At its core, Gemma DE is a two-file Shiny app. This means that the user interface is defined in `ui.R` and the server logic is defined in `server.R`. There are also three additional files for organization: `process.R` which holds the algorithms that drive the condition enrichment, `renderTools.R` which contains scripts necessary for visualizations, and `gemmaAPI.R`, which contains logic to interface with Gemma.

# Running Gemma DE
After installing the necessary dependencies (note that some dependencies are from my personal GitHub repository, not CRAN), sourcing `start.R` will sequentially load the required packages, data and then start the web server. Assuming the app is run on a dedicated server, users can connect on their local machine by SSH tunneling. This is as simple as executing the following command in a command prompt: `ssh -L <LOCAL PORT>:localhost:<GEMMA DE PORT> <USERNAME>@<GEMMA DE SERVER> -p <GEMMA DE PORT>`, replacing items in triangle braces with the correct options. Once completed, you can connect to the user interface by opening a web browser and navigating to `localhost:<LOCAL PORT>`.

# Using Gemma DE

## Through the User Interface
Genes can be inputted (either by entry into the textbox or by file upload) as NCBI or Ensembl gene IDs, official symbol, description, known aliases or by GO group. An optional differential expression signature can also be inputted (in the same order as the genes) to attempt to find biological conditions whose fold changes are correlated with the query. There are more options which can be expanded by clicking on the corresponding text to restrict which categories are reported, modify filtering options and etc. Of note in this category is the scoring method. For searches not involving a query differential expression signature, we recommend leaving this on the default. Searches involving a query signature should use the modified vector space model (M-VSM) or correlation (we describe these in our publication). Clicking search will begin an enrichment and update the URL, which can be saved for future reference.

Many visualization options beyond just the results table are available. These are accessible either via. the tabs (located above the results table) or by clicking on the "Visualize" button. Clicking the "Visualize" button will send a request to the Gemma database to fetch gene expression information (and may take some time to retrieve, depending on how much data is being fetched), while all tab options recalculate views based on already available data (and are thus likely faster to compute).

Data can also be downloaded by clicking the "Download" button (for the results table), or the download icon (in the "Visualize" view).

## Programmatically
We are currently working on a REST API for accessing Gemma DE.

# Evaluation of Gemma DE
The performance of Gemma DE was evaluated in large part by running enrichments on either entirely or partially simulated data. Other metrics were obtained on real data by searching gene sets with known expression differences in some biological processes (ie. biological sex).
