# Computational Geosciences
Welcome to the Computational Geosciences resource at UiS. This is an educational project funded by the Faculty of Science and Technology, at the University of Stavanger, Norway (UiS). The resource is made by students and faculty from the departments of Energy Resources (IER) and Mechanical and Structural Engineering (IMBM). 

Students: Angela Hoch (IER), Adham Amer (IMBM), Vania Mansoor and Linda Olsen (IER). 

Faculty: Nestor Cardozo (IER), Lisa Watson (IER), Wiktor Weibull (IER), and Knut Giljarhus (IMBM). 

Please feel free to use this material for teaching and research. If you have any comments or want to contribute to the resource, please contact us at nestor.cardozo@uis.no

## Manual: Working with the resource
The resource consists of a book in pdf format (compGeo.pdf), and a source folder where the data, functions and notebooks related to the resource are included.

The programming language of choice is Python, and our approach is as follows: We introduce briefly the theory and applications, implement them in Python functions, and illustrate them using Jupyter notebooks. The best way to work with the resource is to clone this repository.

### Clone the repository
This saves the repository to your local machine. It behaves almost like a copy.
1. Open a terminal.
2. Navigate to the folder where you would like to store the local copy of the repository. (`cd <foldername>`)
3. Press the green button 'Clone or download' on the right hand side and copy the path.
4. Execute the terminal command `git clone <fill in path here>`.
5. Now you can start working with the resource files. They are saved in the folder you chose in the 2. step.

### Updating your local files
Once in a while you should update your local files to the latest changes in the repository. This is important since we will be including new chapters along.
1. Open a terminal.
2. Navigate inside the folder of the repository on your machine. (`cd <foldername>`)
3. Execute the terminal command `git pull`

## Background information
The notebooks follow the directory structure of the resource, which is based on separate data, functions and notebooks folders. We recommend that you follow the same directory structure when running the notebooks. When necessary, we have imported the Numpy and Matplotlib libraries in our functions and notebooks as:

`import numpy as np`

`import matplotlib.pyplot as plt`

We also use other libraries such as Cartopy, pyproj, and uncertainties by Eric O. Lebigot. 

To work with the jupyter notebooks of the resource, open jupyter and navigate to the corresponding folder of your local copy of the resource.

## Current state
Currently, eight chapters are published. Upcoming chapters are: elasticity (Ch. 9), and the inversion problem (Ch. 10).
