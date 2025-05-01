# Shiny_UNITE Tutorial

## Building a Shiny App to Match Taxa Against the UNITE Fungi Database

### Goal
This Shiny app checks which species from a list (e.g., UK lichens, NBN fungi) are represented in the UNITE database of ITS fungal sequences and optionally integrates GBIF/NBN synonym data.

### What You'll Need

- R & RStudio installed
- Required R packages (install if needed):

```
install.packages(c("shiny", "dplyr", "readr", "stringr"))
```

### Files you'll use
- UNITE FASTA header file (.txt) or processed .csv file (e.g. UNITE_public_19.02.2025_headers_rfmt.csv)
- Taxon spreadsheet (e.g. lichen-list-201803.csv or fungal species.csv)

**Optional:** Synonym file (e.g. 2025-04-30-uksi_species_synonyms(Sheet1).csv)


### Your shiny app will have the following directory structure

```
UNITE_TaxonMatcher_App/
├── app.R
└── data/
    ├── UNITE_public_19.02.2025.csv
    ├── UNITE_public_2024.04.11.csv
```

Where app.R is the Shiny App File and data are the header files for the different UNITE datasets that will be built in. The user will be able to select one of these in the form of a drop-down menu.

## Step 1: Create Your Shiny App File

At the top specify the libraries to load up
```
library(shiny)
library(dplyr)
library(readr)
library(stringr)
```

Then load up the available UNITE files from the `data/` folder 

```
unite_versions <- list.files("data", pattern = "\\.csv$", full.names = FALSE)
names(unite_versions) <- gsub("_", " ", tools::file_path_sans_ext(unite_versions))
```

Next define your UI: 

```
ui <- fluidPage(
  titlePanel("UNITE Database Taxon Matching"),
  sidebarLayout(
    sidebarPanel(
      selectInput("unite_version", "Select UNITE Version:", choices = unite_versions),    #selectInput specifies a drop-down menu option
      fileInput("taxa_file", "Upload Taxon List CSV:", accept = ".csv"),      #this will be the required input spreadsheet from the user
      fileInput("syn_file", "Upload GBIF/NBN Synonyms CSV (optional):", accept = ".csv"),      #this is the optional second file to check for any synonyms
      actionButton("run", "Run Matching"),        #This is the button to run the order of operations that we will specify next 
      downloadButton("download", "Download Matched Data")    #This will allow the user to download the output 
    ),

#The mainPanel then specifies the output as a table that is viewed within the app
    mainPanel(
      h4("Preview of Matching Results:"),
      tableOutput("preview"),
      verbatimTextOutput("summary")
    )
  )
)
```

The next step is to specify the functions and order of operations:

```
server <- function(input, output, session) {
  observeEvent(input$run, {
    req(input$taxa_file)

    # Load selected UNITE version
    unite_path <- file.path("data", input$unite_version)
    df <- read.csv(unite_path, stringsAsFactors = FALSE)

    # Load and clean taxon data
    taxa_df <- read.csv(input$taxa_file$datapath, stringsAsFactors = FALSE)
    taxa_df$Taxon.name <- gsub("var\\.|f\\.|s\\.|str\\.|lat\\.|subsp\\.", "", taxa_df$Taxon.name)
    taxa_df$Taxon.name <- trimws(taxa_df$Taxon.name)

    taxa_df$UNITE.Record <- taxa_df$Taxon.name %in% df$Species
    if ("Recent.Synonym" %in% colnames(taxa_df)) {
      taxa_df$UNITE.Record.Synonym <- taxa_df$Recent.Synonym %in% df$Species
      taxa_df[is.na(taxa_df$Recent.Synonym), "UNITE.Record.Synonym"] <- NA
    }

    if (!is.null(input$syn_file)) {
      syn_df <- read.csv(input$syn_file$datapath, stringsAsFactors = FALSE)
      item_counts <- sapply(unique(taxa_df$Taxon.name), function(x) sum(df$Species == x))
      df2 <- data.frame(item = unique(taxa_df$Taxon.name), count = item_counts)

      taxa_df <- taxa_df %>%
        left_join(df2, by = c("Taxon.name" = "item")) %>%
        left_join(syn_df, by = c("taxonID" = "taxon_id"))

      colnames(taxa_df)[ncol(taxa_df)] <- "synonym.name"
      taxa_df$Synonym.UNITE.Record <- taxa_df$synonym.name %in% df$Species
    }

    output$preview <- renderTable({
      head(taxa_df, 20)
    })

    output$summary <- renderPrint({
      list(
        UNITE_matches = sum(taxa_df$UNITE.Record, na.rm = TRUE),
        Synonym_matches = if ("UNITE.Record.Synonym" %in% colnames(taxa_df))
          sum(taxa_df$UNITE.Record.Synonym, na.rm = TRUE) else NA
      )
    })

    output$download <- downloadHandler(
      filename = function() {
        paste0("matched_taxa_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(taxa_df, file, row.names = FALSE)
      }
    )
  })
}
```

The final step is to put it all together: 

```
shinyApp(ui = ui, server = server)
```
Save the above as `app.R`. This is what will be needed to open the shiny app up in RStudio. Shiny Apps can be hosted from web servers as well. 


## Step 2: Run Your App

- Open the app.R file in RStudio.
- Click "Run App" in the top right of the script window.

The app will open in a new window or browser tab.

## Step 3: Upload Files & Explore
**Upload UNITE headers**
- As .txt from grep "^>" UNITE.fasta > headers.txt, OR
- As a reformatted .csv file from earlier scripts

**Upload your species list**

*Example:* sssi-lichen-list-201803.csv

**(Optional) Upload synonyms file**

*Example:* 2025-04-30-uksi_species_synonyms(Sheet1).csv

**Click "Run Matching"**


## Database options
|Db option|File source|Description|Reference|
|--|--|----|----|
|sh_general_release|sh_general_release_19.02.2025.tgz|UNITE general FASTA release for Fungi: the RepS/RefS of all SHs, adopting the dynamically use of clustering thresholds whenever available.|Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2025): UNITE general FASTA release for Fungi. UNITE Community. 10.15156/BIO/3301229|
|sh_general_release_s|sh_general_release_s_19.02.2025.tgz| UNITE general FASTA release for Fungi 2: the RepS/RefS of all SHs, adopting the dynamically use of clustering thresholds whenever available. Includes global and 97% singletons. | Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2025): UNITE general FASTA release for Fungi 2. UNITE Community. 10.15156/BIO/3301230|
|sh_general_release_all|sh_general_release_all_19.02.2025.tgz|	UNITE general FASTA release for eukaryotes: the RepS/RefS of all SHs, adopting the dynamically use of clustering thresholds whenever available.| Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2025): UNITE general FASTA release for eukaryotes. UNITE Community. 10.15156/BIO/3301231|
|sh_general_release_s_all|sh_general_release_s_all_19.02.2025.tgz| 	UNITE general FASTA release for eukaryotes 2 : the RepS/RefS of all SHs, adopting the dynamically use of clustering thresholds whenever available. Includes global and 97% singletons. |Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2025): UNITE general FASTA release for eukaryotes 2. UNITE Community. 10.15156/BIO/3301232Z|
|UNITE_public_19.02.2025 |UNITE_public_19.02.2025.fasta.gz| Full UNITE+INSD dataset for Fungi : comprises “all” fungal ITS sequences of the UNITE and INSDC databases present in UNITE Species Hypotheses, updated and released some four times a year. Locked UNITE sequences and low quality (and overly short) INSDC sequences are however excluded. | Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2025): Full UNITE+INSD dataset for Fungi. UNITE Community. 10.15156/BIO/3301227 |
|	UNITE_public_all_19.02.2025 | UNITE_public_all_19.02.2025.fasta.gz |Full UNITE+INSD dataset for eukaryotes: comprises “all” eukaryote ITS sequences of the UNITE and INSDC databases present in UNITE Species Hypotheses, updated and released some four times a year. Locked UNITE sequences and low quality (and overly short) INSDC sequences are however excluded.| 	Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2025): Full UNITE+INSD dataset for eukaryotes. UNITE Community. 10.15156/BIO/3301228|

## What You Get
A table preview of:
1. Your taxon list annotated with matches to UNITE species
2. A summary count of matched species
3. An option to download results as `.csv`

### Tips
- Clean species names before upload to improve match rates.
- The app works best if taxon files contain columns: Taxon.name, Recent.Synonym, and taxonID.
- Make sure column names match exactly; modify the code if your files use different names.
