#!/bin/bash


Rscript="/home/yuchen/miniconda3/envs/R4.0/bin/Rscript --vanilla"
create_obj="/data/bingling/rPackageTutorial/Continued/example/7.1 share_metabolic.tumor.R"

# Cluster labels corresponding to each sample tumor region.
${Rscript} ${create_obj} ST88_20210331 1 15 17 NULL NULL NULL NULL NULL NULL
${Rscript} ${create_obj} ST87_20210331 5 8 18 NULL NULL NULL NULL NULL NULL
${Rscript} ${create_obj} ST06_20210716 2 8 9 11 NULL NULL NULL NULL NULL
${Rscript} ${create_obj} ST103_20210718 7 11 16 10 18 NULL NULL NULL NULL
${Rscript} ${create_obj} ST109_20210330 21 17 10 19 1 16 18 5 9
${Rscript} ${create_obj} ST121_20210806 3 13 14 8 12 NULL NULL NULL NULL
${Rscript} ${create_obj} ST124_20211223 8 5 6 NULL NULL NULL NULL NULL NULL
${Rscript} ${create_obj} ST129_20201214 1 NULL NULL NULL NULL NULL NULL NULL NULL
${Rscript} ${create_obj} ST129_20210428 1 NULL NULL NULL NULL NULL NULL NULL NULL
${Rscript} ${create_obj} ST133_20210429 4 NULL NULL NULL NULL NULL NULL NULL NULL
${Rscript} ${create_obj} ST32_20210807 2 10 6 NULL NULL NULL NULL NULL NULL
${Rscript} ${create_obj} ST35_20210401 3 4 9 12 11 8 2 7 10
${Rscript} ${create_obj} ST49_20211220 6 NULL NULL NULL NULL NULL NULL NULL NULL
${Rscript} ${create_obj} ST69_20211222 3 4 6 NULL NULL NULL NULL NULL NULL
${Rscript} ${create_obj} ST84_20211223 2 3 13 NULL NULL NULL NULL NULL NULL

# If you want data on other type of tissue, replace the label for the corresponding tissue region.
