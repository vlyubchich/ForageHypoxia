# ForageHypoxia

## User memo
This GitHub repository is for collaborating on the analysis of Chesapeake Bay data, specifically the Maryland SeaGrant project on benthic organisms and hypoxia in the bay. The repository is set as private (not publicly accessible).

1. As a collaborator, you have direct access to the repository **vlyubchich/ForageHypoxia**. Clone it to your computer to your desired location -- this action will create a folder **ForageHypoxia** with all the files inside.
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

The main data files are described below.

+ `data_rca/` is a folder with RCA model outputs. Due to their large size, the folder is ignored on Git. Please put the corresponding files to this folder on your computer manually.


## Code

The main files with code are described below, presumably in the order they enter the project workflow.

+ `code/rca_extract.R` file to process RCA outputs and extract: 
    * `data_rca/rca_cells_[version].csv` information on all model cells, where FSM is the land mask variable (`1`=water, `0`=land, `-1`=river BC, `-2`=ocean BC);
    * `data_rca/rca_ts_YYYY_[version].csv` time series for water cells, separated by year, where `[version]` is the version of the extraction code.

