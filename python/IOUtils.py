from ROOT import *
import os,sys,json,time,re
import logging
from shutil import copy
from pprint import pformat

def createDir(myDir, log=None, over=True):
  if log!=None:
    log.info('Creating a new directory: %s', myDir)
  if os.path.exists(myDir):
    if log!=None:
      log.warning("\t This directory already exists: %s", myDir)
    if over:
      # Overwrite it...
      if log!=None:
        log.warning("But we will continue anyway (I will --overwrite it!)")
    else:
      if log!=None:
        log.error(' And so I exit this place...')
      print 'The directory exist and we exit. Dir = ', myDir
      sys.exit(1)
  else:
    try: os.makedirs(myDir)
    except OSError:
      if os.path.isdir(myDir): pass
      else: raise
