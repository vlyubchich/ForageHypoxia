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
D = read.csv("./data_folder/FileName.csv")
```
to save processed data for future use
```{r}
write.csv(RObjectToSave, "./data_folder/FileName.csv", row.names = FALSE)
``
