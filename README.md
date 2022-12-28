# CONTINUED
## CONTINUED: Cluster, integration and annotation of *In situ* metabonomics data
CONTINUED is a algorithm to process *In situ* metabonomics data include reconstruction of histological spatial structure based on single sample clustering, integration of multiple samples and annotation with LC-MS data.
### Installation
description
### How does it work?
In the `example/` folder, we provide plots and files to perform a procession, taken from our data. Alternatively you can use your own data. Start with something not too large, for example 1Mb. And read on.
### Background noise removal
Try out the `Preprocessing.py`. CONTINUED takes two dataframes contained abundance of m/z as input, in `.txt` format. One is generated under lockmass and another is generated under unlockmass. In order to get the coordinates of the tissue image, CONTINUED also needs a m/z and a corresponding threshold as input. The m/z and threshold were selected according to the shape of the HE staining slice. In function `find_interest_factor`, you will get a dataframe `df_data` and two `.png` files in the `mz_from_path`, the images can help you adjust the parameters. You can get three `.txt` files finally as input for `filter.mass.R`.
### Filter m/z
See `filter.mass.R`. It filter m/z max intensity is less than 400.
### Parameter selection
The clustering results of samples change with the parameter npc and resolution. In `Parameters.test.R`, you can test a series of npcs and resolutions. Select the final parameters according to the results.
### create seurat object
See `creat_obj.R`. Create three seurat object `.rds` files with parameters selected in `Parameters.test.R`. You can choose to plot the complete clustering image, the single clustering image or the boundary image.


