# CONTINUED
## CONTINUED: Cluster, integration and annotation of *In situ* metabonomics data
CONTINUED is a algorithm to process *In situ* metabonomics data include reconstruction of histological spatial structure based on single sample clustering, integration of multiple samples and annotation with LC-MS data.
### Installation

#### The easy way

Let `devtools` do the installation for you. Type the following in an R shell

```
install.packages("devtools")
devtools::install_github("3DGenomes/binless",subdir="binless")
```

#### Manual installation

You can also install it manually as follows:

- Download the latest release [here](https://github.com/3DGenomes/binless/releases/latest)
- Unpack it in a clean folder on your machine, and hit `R CMD INSTALL binless` in a shell.
- If it complains that some packages are not installed, you must install them in R using `install.packages`

### How does it work?

In the `example/` folder, we provide datasets and scripts to perform a procession, taken from our data. Alternatively you can use your own data. Start with something not too large, for example 1Mb. And read on.
### Background noise removal
Try out the `Preprocessing.py`. CONTINUED takes two dataframes contained abundance of m/z as input, in `.txt` format. One is generated under lockmass and another is generated under unlockmass. In order to get the coordinates of the tissue image, CONTINUED also needs a m/z and a corresponding threshold as input. The m/z and threshold were selected according to the shape of the HE staining slice. In function `find_interest_factor`, you will get a dataframe `df_data` and two `.png` files in the `mz_from_path`, the images can help you adjust the parameters. You can get three `.txt` files finally as input for `filter.mass.R`.
### Filter m/z
See `filter.mass.R`. It filter m/z max intensity is less than 400.
### Parameter selection
The clustering results of samples change with the parameter npc and resolution. In `Parameters.test.R`, you can test a series of npcs and resolutions. Select the final parameters according to the results.
### Reconstruction of spatial distribution pattern of organization
See `creat_obj.R`. Create three seurat object `.rds` files with parameters selected in `Parameters.test.R`. You can choose to plot the complete clustering image, the single clustering image or the boundary image. If you have a lot of samples, `create.obj.with.final.parameters.sh` can help you run more conveniently.
### Plot 3D intensity
See `Total.intensity.3D_plot.R`. Here, we create a way to assess the quality of the data. Data quality can be preliminarily judged by observing contour shape of 3D map. In our data, if the contour shape presents regular horizontal lines, it is likely to be bad. We will eliminate such data to ensure the accuracy of subsequent analysis.
### Integration and annotation of multiple samples of DESI-MS
We achieve the clustering of m/z by fitting the kernel density function for all m/z values. We also annotate m/zs as specific metabolites. See `merge.mass.pl`. This requires LC-MS databases, either from your own lab, or from an online database. You will get a dataframe containing clustering and annotation information with prefix `colon_cancer_desi`. This is the most important result.
### Plot target m/z
See `Plot.target.mass.R`. You can visualize the distribution of the same metabolite in different samples.
### Questions? Problems?
Check out the [FAQ](https://github.com/Tang-Lab-super/CONTINUED/labels/FAQ) in the GitHub issues tracker. If you don't find anything, just file an issue by clicking on the "New issue" button. You can use it to report bugs, ask questions or suggest new features. We are looking forward to your feedback!



