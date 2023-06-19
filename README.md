# AKOYA-MIA
Data integration across different modalities.

For count data, refer to shared folder in CloudStor.

Run the Rscript commands inside the RStudio docker: 

docker run -d --name sigcm -p 8787:8787  -v "\<source folder>\:/home/rstudio/singlecellmel" -e USER=rstudio -e PASSWORD=singlecellmel123 adyprat/singlecellmel

To view integrated .h5ad object, use cellxgene to launch:
--port 54890 <directory>/P1_Transferred_092022_1.h5ad
--host 0.0.0.0 --port 64501 <directory>/P3_Transferred_092022_1.h5ad
--host 0.0.0.0 --port 57790 <directory>/P2_Transferred_092022_1.h5ad
--host 0.0.0.0 --port 58456 <directory>/P5_Transferred_120522.h5ad
--host 0.0.0.0 --port 53759 <directory>/P6_Transferred_120522.h5ad

