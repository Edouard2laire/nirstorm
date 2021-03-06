# NIRSTORM

![NIRSTORM - a Brainstorm plugin for fNIRS data analysis.](https://user-images.githubusercontent.com/24530402/98286213-a4b4a600-1f71-11eb-8a76-8eff3f820e3a.png)

Current features include:
- classical within-subject analysis comprising motion-correction, MBLL and window-averaging.
- statistical analysis using GLM
- optimal montages that optimize the sensitivity to a given cortical
region of interest
- source reconstruction (cMEM, MNE)
- precomputed fluence template based on Colin27

## Authors

 * Thomas Vincent, EPIC center, Montreal Heart Institute, Montreal, Canada
 * Zhengchen Cai, PERFORM Centre and physics dpt., Concordia University, Montreal, Canada
 * Alexis Machado, Multimodal Functional Imaging Lab., Biomedical Engineering Dpt, McGill University, Montreal, Canada
 * Edouard Delaire, PERFORM Centre and physics dpt., Concordia University, Montreal, Canada
 * Robert Stojan, Sportpsychology, Chemnitz University of Technology, Germany
 * Louis Bherer, Centre de recherche EPIC, Institut de Cardiologie de Montréal, Montréal, Canada
 * Jean-Marc Lina, Electrical Engineering Dpt, Ecole de Technologie Supérieure, Montréal, Canada
 * Christophe Grova, PERFORM Centre and physics dpt., Concordia University, Montreal, Canada

## Usage

The main documentation is in the form tutorials available on the [nirstorm github project wiki](https://github.com/Nirstorm/nirstorm/wiki#tutorials).

## Installation

[Brainstorm](http://neuroimage.usc.edu/brainstorm/) must be installed prior to installing nirstorm. 

### From Brainstorm (since Bainstorm v. 3.201110 (10-Nov-2020)) ** Recommended **
In Brainstorm, go to the update menu at the top of the brainstorm windows, then 'Update NIRSTORM' and finally 'Download NIRSTORM'. The latest version of NIRSTORM will then be downloaded and installed.  Updates of NIRSTORM can then be downloaded using the same menu. 

### From Github 

The script `nst_install.m` takes care of copying or linking processes and functions into the brainstorm user folder.

Parts of the nirstorm plugin may already be shipped with the lastest brainstorm version and are available in the process selection menu in the "NIRS" submenu (modified bear lambert law and bad channel tagging).
The current installation will override them.

Run brainstorm before installing nirstorm.
All commands indicated below must be run in the nirstorm folder where the archive was uncrompressed.

### Copy installation (windows, linux)

To copy all processes and functions into the brainstorm user folder, run under matlab:
```matlab
>> nst_install('copy');
```
When updates are downloaded, this installation command must be run again for changes to take effect.

### Linked installation (linux only)

To create symbolic links of all processes and functions into the brainstorm user folder, run under matlab:
```matlab
>> nst_install('link');
```
When updates are downloaded, this installation command has to be run again only if there are new files.


