# UEL_DLO_Article
## Description
This code contains UEL subroutine for simulation of DLO couplded mechanical field. The Abaqus input file is also provided to be used with user element subroutine for modeling coupled multi-physics.

## Usage
Please save the Abaqus input file (.inp) and user element subroutine (.for) in a folder. Open the folder and run \cmd.exe. Change the directory to the folder using 
cd (file path). So, the log file and odb file will be saved in the folder. In the command prompt write the following conde:

abaqus job=(jobe name) user=(user subroutine name)

Then, you can run the Abaqus an open the .odb file to see the results.

## Citation
If you use this code in your research or publications, please cite the following article, which is going to be published in the journal of Mechanics of Materials:

"Constitutive modeling of diffusion-limited oxidation coupled with a large deformation theory for polymer degradation", Hossein Naderi, Roozbeh Dargazany.
