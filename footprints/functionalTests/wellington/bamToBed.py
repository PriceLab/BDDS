# bamToBed.py
#--------------------------------------------------------------------------------
from subprocess import call
import shutil
import unittest
import sys
import os
#--------------------------------------------------------------------------------
executableName = "bamToBed"
#--------------------------------------------------------------------------------
def getProgramPath():

   print("checking for valid path to executable %s" % executableName)
   path = shutil.which(executableName)
   assert(isinstance(path, str))
   assert(len(path) >= len(executableName))

   return(path)

#--------------------------------------------------------------------------------
def runTask(bamFile):
   
   pathToExecutable = getProgramPath()
   call([pathToExecutable, "-i",  bamFile])

#--------------------------------------------------------------------------------
if(len(sys.argv) != 2):
  print("usage: bamToBed.py someFile.bam")
  sys.exit(1)
   
bamInputFile = sys.argv[1]
assert(os.path.isfile(bamInputFile))
runTask(bamInputFile)
 
