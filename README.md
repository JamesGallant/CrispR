# CrispR
 Engine for finding guide RNA's globally. This package uses CRISPRseek to generate the guide RNA's which will be handled in python to select candidates. 
 
# Requirements
R, python

# How to install
## Windows
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

# Contribution of code to this package
Clone the repository and checkout a new branch. 

```
git clone -recursive https://github.com/vumc-mmi/CrispR.git
git checkout -b yourversion
```
Make your changes and run the unit tests to verify that everything still works. Also add unit tests for any new critical functions. 
```
py -m unittest . discover
```
If everything passes then submit a pull request.

# To do
* Create python bindings for the TxDB and guide RNA R scripts. Also unit tests to verify that the files are being made, this should be enough. 
* Create python utility scripts to handle the output from R and choose optimal gRNA's
