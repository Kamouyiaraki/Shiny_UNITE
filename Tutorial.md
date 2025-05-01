# Shiny_UNITE Tutorial

## Building a Shiny App to Match Taxa Against the UNITE Fungi Database

### Goal
This Shiny app checks which species from a list (e.g., UK lichens, NBN fungi) are represented in the UNITE database of ITS fungal sequences and optionally integrates GBIF/NBN synonym data.

✅ Uses only one input file (user uploads a taxon list CSV).

✅ Allows user to select the taxon name and taxon ID columns from that file.

✅ Checks if taxon names and synonyms (from within that same file) are found in the selected UNITE database version.

✅ Does not allow the user to upload a second synonym file.

✅ Provides a clean UI, summary, and downloadable output.



### What You'll Need

- R & RStudio installed
- Required R packages (install if needed):

```
install.packages(c("shiny", "dplyr", "readr", "stringr"))
```

### Files you'll use
- The shiny app folder that includes the app R file and `data` directory
- Taxon spreadsheet (e.g. lichen-list-201803.csv or fungal species.csv)


### Your shiny app will have the following directory structure

```
UNITE_TaxonMatcher_App/
├── app.R
└── data/
    ├── UNITE_public_19.02.2025.csv
    ├── UNITE_public_all_19.02.2025.csv
...
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
      selectInput("unite_version", "Select UNITE Version:", choices = unite_versions),    #selectInput, "choices =" specifies a drop-down menu option
      fileInput("taxa_file", "Upload Taxon List CSV:", accept = ".csv"),      #this will be the required input spreadsheet from the user

      uiOutput("column_select_ui"), #once the file is uploaded the ui outputs to a column selection step

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

## input requirements: step 1 is the inputfile and then in reaction to that file being uploaded then the CSV is read 
   taxa_data <- reactive({
    req(input$taxa_file)
    read.csv(input$taxa_file$datapath, stringsAsFactors = FALSE)
  })

# Dynamically create column selectors once file is uploaded
  output$column_select_ui <- renderUI({
    req(taxa_data())
    colnames_df <- names(taxa_data())

## tagList() creates a simple list of tags to be seen in the ui ("choices =" specifies a drop-down option)
    tagList(
      selectInput("taxon_col", "Select Taxon Name Column:", choices = colnames_df),
      selectInput("synonym_col", "Select Synonym Column (optional):", choices = c("None", colnames_df))
    )
  })

observeEvent(input$run, {
    req(taxa_data(), input$taxon_col)

    # Load selected UNITE database version
    unite_df <- read.csv(file.path("data", input$unite_version), stringsAsFactors = FALSE)

    # Clean species names
    unite_df$Species <- str_trim(unite_df$Species)

    # Work with user data
    taxa_df <- taxa_data()

    # Clean taxon names
    taxa_df$CleanTaxon <- gsub("var\\.|f\\.|s\\.|str\\.|lat\\.|subsp\\.", "", taxa_df[[input$taxon_col]])
    taxa_df$CleanTaxon <- trimws(taxa_df$CleanTaxon)

    taxa_df$UNITE.Record <- taxa_df$CleanTaxon %in% unite_df$Species

# If synonym column selected and not "None"
    if (!is.null(input$synonym_col) && input$synonym_col != "None") {
      taxa_df$CleanSynonym <- gsub("var\\.|f\\.|s\\.|str\\.|lat\\.|subsp\\.", "", taxa_df[[input$synonym_col]])
      taxa_df$CleanSynonym <- trimws(taxa_df$CleanSynonym)
      taxa_df$UNITE.Record.Synonym <- taxa_df$CleanSynonym %in% unite_df$Species
      taxa_df[is.na(taxa_df[[input$synonym_col]]), "UNITE.Record.Synonym"] <- NA
    }

# Output preview
    output$preview <- renderTable({
      head(taxa_df, 20)
    })

    # Output summary stats
    output$summary <- renderPrint({
      num_direct <- sum(taxa_df$UNITE.Record, na.rm = TRUE)
      num_syn <- if ("UNITE.Record.Synonym" %in% names(taxa_df))               sum(taxa_df$UNITE.Record.Synonym, na.rm = TRUE) else 0
      num_any <- sum(taxa_df$UNITE.Match.Any, na.rm = TRUE)

  list(
    Total_Records = nrow(taxa_df),
    UNITE_Matches = num_direct,
    Synonym_Matches = if (num_syn > 0) num_syn else "N/A",
    Total_Matches_Combined = num_any
  )
})

    # Enable download
    output$download <- downloadHandler(
      filename = function() {
        paste0("UNITE_matches_", Sys.Date(), ".csv")
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
