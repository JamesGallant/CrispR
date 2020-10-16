# Description
Windows application for generating and quering a guide RNA database in bulk and includes primer design. Currently only automatic guide RNA detection and primer design is availible for Streptococcus thermophilus cas9. 

# Requirements
R version 3 or 4 (we tested on version 4), internet connection for certain functions, Windows 10 (older versions may work). 

# How to install
You need to have admin rights to install this software. Ignore warnings, I did not sign this application. 
## Creating the R environment
Install [R](https://cran.r-project.org/bin/windows/base/R-4.0.2-win.exe) on your system and make sure its callable from your command line. To check if this is already availible open command prompt in windows and run the following command:

```
Rscript
```
You should see a menu and no errors if R is callable. If not, R has to be added to your path variable you can add it by following these instructions:

1. Find the location of R's installation on your machine, usually its at `C:\Program Files\R\R-4.0.2\bin\x86`.
2. Open the start menu and type **View advanced system settings**.
3. Click on **Environment variables**.
4. Find **system variables**, **path** and click on **edit**. 
5. Click **new** and add the folder path for R, in this example its `C:\Program Files\R\R-4.0.2\bin\x86`.
6. Verify that R is callable by opening command prompt and running the `Rscript` command.

## Installing stable Crispinator
Download the zip file and extract, double click on the Crispinator-setup.exe and follow the prompts. The program will be installed in the program files 86 directory and the zip file can be deleted

## Uninstall
Crispinator can be uninstalled using the add or remove programs native to windows.
_note_: Windows does not automatically install database files this way. To completely uninstall delete the folder in program files.

# Contribution of code to this package
Crispinator is built using python and R. R is used to run the guide RNA predictions and is based on the CrisprSeek library. For the most part you will be working with the PyQt5 wrapper for Qt. Crispinator.py is the mainWindow. Each button calls a dialog function which can be found in src/dialogs.py. The information from dialog py can be sent back to Crispinator main window (I use dictionaries) and processed. For time consuming functions use threading, the functions used for threading can be found in workers.py. So the worklflow I use is mainwindow -> Dialog -> mainwindow -> worker -> mainwindow

## Prepare the environment
```
git clone -recursive https://github.com/vumc-mmi/CrispR.git
git checkout -b yourversion
pip install -r requirements.txt # installs the python dependencies
```
Some empty directories are required, namely a databases folder in the main folder structure and a folder called local_packages in src

```
Crispinator
     |
     |-databases
     |-gui
         |-Icons
         |-text
     |-src
         |-local_packages
         |-py
         |-R
     |-tests

```
## writing and testing python modules
Some but not all functions have unit tests, using the discover . function may fail some tests. The unit tests are located in the tests directory and can be run from there. Make changes by creating a new unit test class or add to the existing classes and testing these changes using the test runner. Run all unit tests and make sure everything is working before commiting. 
```
py -m unittest discover .
```
## Adding R scripts
The R scripts are located in src/R. Remember to add bioconductor dependencies under install_bioconductor.R and Cran repositories using the following method in the script its self:

```
if(!require(optparse)){
  install.packages("mynewlibrary")
}
```
Call the R scripts using the subprocesses module in python

## Compiling the module
Once done compile the updated application. Any method should be fine but we use pyinstaller

```
pyinstaller --onefile --path /venv/Lib/site-packages --icon icon.ico -w Crispinator.py
```
Copy the Crispinator.exe from the dist folder to the new main directory and delete the dist and build folder. 

I used [inno installer](https://jrsoftware.org/isinfo.php) to make the setup.exe, the script used to compile is inno.iss. Run this using inno installer to create a standalone setup.exe to distribute

# Credits
CrisprSeek, GenomicRanges, BsGenome, PyQt5, qDarkStyle

