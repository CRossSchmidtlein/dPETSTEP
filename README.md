dPETSTEP

This application is written in Matlab and designed as an extension of PETSTEP. Currently, Matlab functions from PETSTEP are used by dPETSTEP, but there is no requirement of CERR (see below).

dPETSTEP is a fast dynamic PET simulator designed for high throughput simulation of dynamic PET images. It allows full simulation of a user defined parametric image, kinetic model, input function, time sampling and more into realistic (noisy) dynamic PET-like images. Furthermore, the dynamic PET data can then be model fitted to produce parameter/parametric image estimates.

This tool is described in Häggström, Ida, et al. "Dynamic PET simulator of tracers via emission projection for kinetic modeling and parametric image studies." Medical Physics 43.6 (2016): 3104-3116.


PETSTEP

Positron Emission Tomography Simulator of Tracers via Emission Projection

This application is written in Matlab and designed as a plugin for CERR (Computational Environment for Radiological Research) also in Github under an GNU GPL license. https://github.com/adityaapte/CERR

PETSTEP is a fast PET simulator designed for high throughput simulation of PET images. It allows a full simulation of a user defined activity distributions (objects) or the insertion of realistic user defined tumors into exisiting patients. Both patient and objects must be in dicom format loaded in CERR. Codes are found in https://github.com/CRossSchmidtlein/PETSTEP.

This tool is described in Berthon, Beatrice, et al. "PETSTEP: Generation of synthetic PET lesions for fast evaluation of segmentation methods." Physica Medica 31.8 (2015): 969-980.