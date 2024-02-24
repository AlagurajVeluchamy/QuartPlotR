# Welcome to QuartPlotR <img src='logo.png' align="right" height="130"/>

## 1. Introduction


## QuartPlotR: A quarternary phase diagram tool; application on stroke, cancer genomics and metagenomics data

## 2. Installation

```
install.packages('QuartPlotR')
```

## 3. Usage

1. **To run QuartPlotR locally**

This option is particularly **recommended** when you have a large number of QC samples to test.

Run the following script, the Web App will pop up and you can then enjoy using QuartPlotR

```
QuartPlotR::runGui()
```

2. **Use the website**

Alternatively, you can also access QuartPlotR vis [https://bcdd.shinyapps.io/QuartplotR/](https://bcdd.shinyapps.io/quartplotR/). 

> **Note** 
> 
> RawHummus is deployed for free at https://www.shinyapps.io. It allows to use 1024 MB of memory. Therefore, a large numbers of files will not be uploaded and/or analyzed successfully. In this case, please run RawHummus locally.

3. **Demo Data**

A list of demo files have been provided, including log files and raw data files. Please use this [link to download](https://github.com/....alaguraj/quartplotr) the demo data.

## 4. Main Functions

**Figure 1.** Overview of the main functions in...

<img src="https://github.com/Alagu/Fig1.png" width = "75%"/>

                                                          

**How QuartPlotR Works?**

In addition 



## 5. Docker 

Dockerized QuartPlotR App is available on [Docker Hub](https://hub.docker.com/r/alagu/quartplotr)

To run this Shiny App on your computer:

```
docker pull alagu/quartplotr
docker run --rm -p 3838:3838 alagu/quartplotr
```
...
and it will be available at http://0.0.0.0:3838

## 6. How to cite

If you find this shiny App usful, please consider citing it:

- [Alaguraj Veluchamy](https://doi.org/....)
