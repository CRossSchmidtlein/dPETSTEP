ex_dPETSTEP_dynamicImage

An example data file and a script to run a simulation and subsequent model parameter fitting of the simulated data.

1. Unpack the zip file. It contains three files, one with data, one with a main script, and one with loadable settings for GUI usage.

2. Load the data into Matlab:

```
   >> load data_ex_dPETSTEP_dynamicImage.mat
```

3. Run the simulation either:

  3.1. Via GUI in Matlab:

```
   >> dPETSTEPgui_sim
```

  load the settings .xls file by pressing "Load" button in GUI

  3.2. Via the main script from Matlab prompt:

  Open up and adjust the settings in file "Dynamic_setSimParameters.m" prior to running main file.

```
   >> main_ex_dPETSTEP_dynamicImage
```
