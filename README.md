# AKOYA-MIA
Data integration across different modalities.

For count data, refer to shared folder in CloudStor.

Run the Rscript commands inside the RStudio docker: 

docker run -d --name sigcm -p 8787:8787  -v "\<source folder>\:/home/rstudio/STvEA" -e USER=rstudio -e PASSWORD=singlecellmel123 adyprat/singlecellmel
