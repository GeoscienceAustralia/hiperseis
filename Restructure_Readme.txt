===============
Software Package for Passive-Seismic Data processing, analysis and modelling
===============

https://docs.google.com/document/d/1vOX8F_FRYE_1v16vODJj-BwEyb88LYHTRKiBqckI0-E/edit?usp=sharing


Restructure, cleanup and document the passive-seismic Github repository
New Target Structure:

Created a folder “scripts”
Move all the miscellaneous, smallish, one-off,  or unknown scripts/folders into scripts
Keep the original folder structure as they are.
All legacy scripts (unsure where to put) can be dumped into the “script” folder.


The folder seismic as a top level main package name
Move major packages of the seismic pipeline into here
In the future, Refactor-rename this seismic folder to  seismicpy package, so that we can release passive-seismic package as seismicpy.*


Notebooks: store jupyter notebooks.


docs folder to be used for more comprehensive docs of this system (sphinx)

tests folder: All unit-tests should be in this folder.


There are three folders, which are still left as where they are in the top-level dir.
Two of them are actually linked to separate github repos.
They are like third-party software.
We need to think what is the best way to handle these folders:
ellip-corr @ 463e898
PhasePApy @ 87bac18
Iloc_rstt
TODO:

Top level Readme, as well as package-module level docs
Sphinx auto build docs and publish with readthedocs.
More unit tests
Automatic build with Travis  CI
Installation scripts,
Workflow pipeline docs and example/demo scripts



..  This text will not be shown 
    https://circleci.com/gh/GeoscienceAustralia/passive-seismic.svg?style=shield
    https://circleci.com/gh/GeoscienceAustralia/passive-seismic



