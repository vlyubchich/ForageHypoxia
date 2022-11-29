# ForageHypoxia

## User memo
This GitHub repository is for collaborating on the analysis of Chesapeake Bay data, specifically the Maryland Sea Grant project on benthic organisms and hypoxia in the bay. 

1. If you are a collaborator, you have direct access to the repository **vlyubchich/ForageHypoxia**. Clone it to your computer to your desired location -- this action will create a folder **ForageHypoxia** with all the files inside.
If you are not a collaborator, you can fork, edit, then create a pull request for the repository.
2. Normal RStudio workflow:
    + Open in RStudio the project **ForageHypoxia**.
    + From the RStudio **Git** tab, "Pull" the project to get the most recent version of the files on your computer.
    + Edit, add or delete files, such as directly in the folder **ForageHypoxia**. Save your changes.
    + The changes will also appear in the **Git** tab -- check that sensitive, large, or non-plain text (DOCX / PDF / XLSX / etc.) files do not appear there. If they do, move the files to an ignored folder or add exceptions to the **.gitignore** file.
    + Put checkmarks to "Stage" the changes.
    + Click "Commit" and add a brief description ("Commit message") of the changes you are uploading. Commit the changes -- now they are saved in the version control on your computer.
    + Click "Push" to upload the changes to GitHub. Now all the collaborators can see and "Pull" the updates.
    
Note that file paths are relative to the project location on your computer, which is the default R working directory for this project. Hence, there should be no need to set a working directory. For example, to load data, use 
```{r}
D <- read.csv("./data_folder/FileName.csv")
```
to save processed data for future use
```{r}
write.csv(RObjectToSave, 
          "./data_folder/FileName.csv", 
          row.names = FALSE)
```

## Data

+ `data_rca/` is a folder with RCA model outputs. Due to their large size, the folder is ignored on Git. Please create a folder with this name and put the corresponding files to this folder on your computer manually. The latest versions of the files can be downloaded from Google Drive (`_2` in the file names denotes the version of the extraction):
    * [rca_cells_2.csv](https://drive.google.com/file/d/1fN1U_pKxkkZ9EqHIMVf5zoAuqHF22JNg/view?usp=share_link), 1 MB, contains information on all model cells, where FSM is the land mask variable (`1`=water, `0`=land, `-1`=river BC, `-2`=ocean BC)
    * [rca_ts_1986-2015_2.zip](https://drive.google.com/file/d/1olRbfDZeov8LvCFU4Yd6n7l4SlHSuNso/view?usp=sharing), 3 GB, time series for water cells, in separate files by year
+ `data_benthos/` folder contains benthic data, mostly extracted from raw Excel files using the code `benthos_extract.R`:
    * `benthos_strata.csv` information about 10 strata of the bay
    * `benthos_biomass.csv` biomass information, note that we need to use random sites `SITE_TYPE == "RANDOM"`    
    * `benthos_taxa.csv` taxonomic classifications of species from the CBTRUST project. Ignore this file, use the updated file with initials "rjw" listed below.
+ `data_benthos/Taxa_group_IDs_updated_rjw_v3_1.csv` an updated file (email 2022-11-22) with benthos-specific classifications by Ryan W. Use the "Aggregate_grp" as the ID column.
+ `data_benthos/stations_cells.csv` benthic stations with corresponding nearest RCA model cells, output of `code/stations_cells.R`
+ `outputs_Dan_2022-08/` folder with codes and outputs of the summer-2022 volunteer Dan McCrary.

## Code

The main files with code are described below, presumably in the order they enter the project workflow.

+ `code/rca_extract_0.R` contains first attempts to load RCA data and map cell centers over the Chesapeake Bay outline.
+ `code/rca_extract.R` processes RCA outputs and extracts them for analysis, see the contents of the `data_rca/` data folder described above.
+ `code/hynet_0.R` explores options, on subsampled or simulated data, for creating a hypoxia network.
+ `code/benthos_extract.R` combines the data from raw Excel files from the benthic project with CBTRUST and saves output in `data_benthos/` including `benthos_strata.csv`, `benthos_taxa.csv`, and `benthos_biomass.csv`.
+ `code/stations_cells.R` matches benthic stations with corresponding nearest RCA model cell. The code writes the results into `data_benthos/stations_cells.csv`.

