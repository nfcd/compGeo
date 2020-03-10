# compGeo
Computational Geosciences resource at UiS. An educational project funded by the Faculty of Science and Technology, Departments of Energy Resources (IER) and Mechanical and Structural Engineering (IMBM), University of Stavanger. Authors: Faculty: Nestor Cardozo (IER), Lisa Watson (IER), Wiktor Weibull (IER), and Knut Giljarhus (IMBM). Students: Angela Hoch (IER), Adham Amer (IMBM), and Vania Mansoor (IER).

Please feel free to use this material for teaching purposes. If you have any comments or want to contribute to the resource, please contact Nestor Cardozo at [mailto](nestor.cardozo@uis.no).

## Manual: Working with this repository
### Getting started - Clone the repository
This step saves the repository to your local machine. It behaves almost like a copy.
1. Open a terminal.
2. Navigate to the folder where you would like to store the local copy of the repository. (`cd <foldername>`)
3. Press the green button 'Clone or download' on the right hand side and copy the path.
4. Execute the terminal command `git clone <fill in path here>`.
5. Now you can start working on the files. They are saved in the folder you chose in the 2. step.

### Continue working - Pull the repository
Whenever you start working on the files after someone else made changes on the repository, you should update your local files with these changes.
1. Open a terminal.
2. Navigate inside the folder of the repository on your machine. (`cd <foldername>`)
3. Execute the terminal command `git pull`

### Save local changes to the common repository - Commit changes
1. Open a terminal.
2. Navigate inside the folder of the repository on your machine. (`cd foldername`)
3. (optional) Check what files you changed by executing `git status`.
4. Add all the modified files to the common repository by executing `git add -A`or add only one specific modified file to the common repository by executing `git add <filename>`.
5. (optional) Check what files will be updated at the common repository by executing `git status`. The green files will be updated, the red ones will not be updated, but stay modified locally.
6. Update the files on the common repository by executing `git commit -m "<type in a message about what you've done>"`.
7. (optional) Check the branch you want to update your files to by executing `git remote -v`. In the list there should currently just occur a path with the name `origin`.
8. Save your changes to the common repository by executing `git push origin master`.
9. (optional) Check if your files are uploaded to the common repository by executing `git status` again. If the files are not listed anymore, your actions were successful.




## Standard libraries for this project
In order to keep our functions consistent, please use the libraries below and import them with the same name in all the functions. Whenever you feel the need to add another library to the list, don't hesitate to do so.

`import numpy as np`

`import matplotlib as plt`

