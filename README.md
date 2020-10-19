# Description
Windows application for generating and quering a guide RNA database in bulk and includes primer design. Currently only automatic guide RNA detection and primer design is available for Streptococcus thermophilus cas9. 

# Requirements
R version 4 (we tested on version 4), internet connection for certain functions, Windows 10 (older versions may work). 

# How to install
I tried to make it easy but was thwarted by windows trusted developer signatures.  All the components to build from source are here though.

I'll add a zip file with a compiled version and explain the process from there. Under the development section I'll add additional examples of explaining how far I got.
## Creating the R environment
Install [R](https://cran.r-project.org/bin/windows/base/R-4.0.2-win.exe) on your system and make sure its callable from your command line. To check if this is already available open command prompt in windows and run the following command:

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
Download the Crispinator.zip file from the releases page. Extract the contents of this folder somewhere on your PC where you don't need admin rights for convenience. I usually leave it in C:\. There is a lot of files so extract it into a folder.

Find the Crispinator.exe file, this is the application and you can use it from here. However, to make accessing the application convenient, right click on Crispinator.exe and select create shortcut. Move this shortcut to the desktop and rename to Crispinator.exe. You can use this to access the application.
## Uninstall
For now, uninstalling is as simple as deleting the folder containing the Crispinator application as well as the shortcut on your desktop if you made one.

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
pyinstaller --path /venv/Lib/site-packages --icon gui/icon.ico -w Crispinator.py
```
The compiled application is in the dist/Crispinator folder. Copy the src, gui and databases folder into this new directory. Make sure you removed any sqlite3 files from databases or bsgenome packages from src/local_packages. Zip this file and draft a new release

### Notes on more advanced installation attempts. 
Do not use the --onefile parameter of pyinstaller. This slows down the startup but also gets flagged as virus because the application is not signed. In the src directory there are two folders called batch and C# here are the left over scripts to assist in creating launchers for Cripsinator if we ever get this flagging issue fixed. 

The batch creates a pointer to the Crispinator.exe and essentially simulates creating the shortcut. the C# script creates a exe from the batch file and should be used as the launcher. 

To create the launch.exe you need `csc.exe` on your path and the `launch.cmd` as well as the `launch.cs` on in the same folder

This is the coomands to create the file:
```batch
csc.exe /target:winexe launch.cs /win32icon:icon.ico
```

This can be further compiled into a setup.exe using inno installer. However, all of this causes flagging on windows without certs.
# Credits
CrisprSeek, GenomicRanges, BsGenome, PyQt5, qDarkStyle

