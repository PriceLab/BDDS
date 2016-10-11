# bamToBedPreTest.py
#--------------------------------------------------------------------------------
#from subprocess import call
import shutil
#import unittest
#import sys
#import os
#--------------------------------------------------------------------------------
def executableProgramFound(programName):

   print("checking for valid path to executable %s" % programName)
   path = shutil.which(programName)
   return(isinstance(path, str) and len(path) >= len(programName))

#--------------------------------------------------------------------------------
assert(executableProgramFound("wellington_footprints3.py") or executableProgramFound("wellington_footprints.py")) 
