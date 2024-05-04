# Spatial source contribution and interannual variation in deposition of dust aerosols over the Chinese Loess Plateau

This is the workflow used to generate the figures published in ... 
The workflow is customized for the NIRD toolkit enivorment, and uses absolute path for the data. 
Therefore these need to be updated if the workflow is to be run on other machine where paths to the data files are different
The paths used in the workflow is configure in the `Master_thesis_UiO_workflow/config` directory. 

The computing enviroment can be installed using conda there is also a docker enviroment which can be used https://hub.docker.com/repository/docker/ovewh/thesisdocker/general.
The easiest would be to use the docker container, data files are available in this zenodo repository. The workflow can configured to work on any computing system you the same thesisdocker container if the path are updated. 

## Setting up computing evniroment. 

To pull and isntall the docker image:

```
docker pull thesisdocker:latest
```

## Running the workflow

To run the workflow simply within the thesis docker container enviroment or on nird toolkit, simply activite the enviroment by `source activate dust`. Then run the workflow you need to be in the `AGU_JGR_CLP_SOURCE_WORKFLOW/Master_thesis_UiO_workflow` directory, the first step of runnning the workflow is done by first checking if all the files are there:

```shell
snakemake -n plot_all_agu 

```
This will excute a dry run, if thtis work without any errors then the paths have been set correctly. Running the workflow is then just a matter of typing

```shell
snakemake -j2 plot_all_agu
```

The `-j2` flag request 2 processors. The `--rerun-triggers mtime` can be added avoid generating all the intermidiate files for making the plots. 
