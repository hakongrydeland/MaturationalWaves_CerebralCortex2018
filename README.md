# MaturationalWaves_CerebralCortex2018

This code underpins analyses in *Waves of Maturation and Senescence in Micro-Structural MRI Markers of Human Cortical Myelination over the Lifespan*, by Grydeland et al., published in Cerebral Cortex, 2018.

### What's in this repo

The main analysis file is [`spline.derivative.R`](spline.derivative.R).
It is provided in this GitHub repository if you'd like to see the exact `R` code that was run, or to conduct similar analyses on your own data.

#### Data

We are not able to share data from this study as participants did not provide informed consent to make data - including derived data - publicly available.

To help you read the code: the `spline.derivative.R` file expects a matlab `.mat` file that contains the T1/T2 ratio for each of your participants across a parcellation of the brain (note that this file can contain any type of brain imaging data from any parcellation). In this paper we used the [Glasser parcellation](https://doi.org/10.1038/nature18933).

The `myelin_fraction70.mat` file used for this study also contains a few demographic variables such as `age`, `sex`, `IQ`, `motion` etc. If you would like to explore relationships with cortical thickness, a `corrthickness.mat` file is expected.

Data was processed using custom matlab scripts.
These scripts are available on request from Håkon Grydeland: [hakon.grydeland@psykologi.uio.no](mailto:hakon.grydeland@psykologi.uio.no).

Data may be available under a managed access agreeement to researchers looking to extend the work published in the Cerebral Cortex paper.
Please contact Håkon Grydeland at [hakon.grydeland@psykologi.uio.no](mailto:hakon.grydeland@psykologi.uio.no) for more information.

#### Code 

The main analysis file is [`spline.derivative.R`](spline.derivative.R).

The other `.R` files in this repository are called by `spline.derivative.R`.

### Support

This project is not being actively developed as a software package.

However, if you have questions about running similar analyses or better understanding the results published in Cerebral Cortex, please contact Håkon Grydeland at [hakon.grydeland@psykologi.uio.no](mailto:hakon.grydeland@psykologi.uio.no).

### Who wrote the code

Most of the code in this repository was written by Håkon Grydeland but it built on work generously shared by Aaron Alexander-Bloch, Ameera Patel, Rafael Romero-Garcia, František Váša, Petra Vértes, Kirstie Whitaker.

Thank you!
