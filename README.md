# UEL_DLO_Article
## Description
This code contains the UEL subroutine for the simulation of the DLO coupled mechanical field. The Abaqus input file is also provided to be used with the user element subroutine for modeling coupled multi-physics.

## Usage
Please save the Abaqus input file (.inp) and user element subroutine (.for) in a folder. Open the folder and run \cmd.exe. Change the directory to the folder using 
cd (file path). So, the log file and odb file will be saved in the folder. In the command prompt, write the following code:

abaqus job=(jobe name) user=(user subroutine name)

Then, you can run Abaqus and open the .odb file to see the results.

## Citation
If you use this code in your research or publications, please cite the following article, which was published in the Journal of Mechanics of Materials:

Naderi, H. and Dargazany, R., 2025. Constitutive modeling of diffusion-limited oxidation coupled with a large deformation theory for polymer degradation. Mechanics of Materials, p.105270.
