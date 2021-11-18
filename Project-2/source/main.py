from sys import path
import yaml
import glob
import os
import numpy as np
# from utils.functions import stiffnessmatrix as sm
# from utils.functions import custom_plot 
import matplotlib.pyplot as plt
import csv

# Make sure your working directory is same as the scripts locations
path = os.path.dirname(__file__)

# Changing the CWD to the current working directory
os.chdir(path)

# loading the YAML files in the the data directory
files = glob.glob(pathname="../data/*.yml")
temp = open(files[0],"r")

# YAML loads the files a dictionary
config = yaml.full_load(temp)

print("This is the end")