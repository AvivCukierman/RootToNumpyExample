from numpy import load,log,linspace,digitize,array,mean,std,exp,all,average,sqrt,asarray,sign,zeros,histogram2d,arange,save
import os
import numpy as np
from operator import sub
from optparse import OptionParser
from copy import copy

parser = OptionParser()

# job configuration
parser.add_option("--inputDir", help="Directory containing input files",type=str, default="../data")
parser.add_option("--submitDir", help="Directory containing output files",type=str, default="../output")
parser.add_option("--numEvents", help="How many events to include (set to -1 for all events)",type=int, default=100000)
parser.add_option("-i","--identifier", help="sample identifier",type=str, default="myjets")
parser.add_option("--config", help="Config file containing list of branches to read in",type=str,default="simpleConfig.json")

# Root configuration
## Reconstructed jets and matched truth jets
parser.add_option("--jetpt", help="reco jet pT branch name",type=str, default="j0pt")
parser.add_option("--tjetpt", help="matched truth jet pT branch name",type=str, default="tj0pt")

# analysis configuration
parser.add_option("--minpt", help="min truth pt", type=int, default=0)
parser.add_option("--maxpt", help="max truth pt", type=int, default=float('inf'))

(options, args) = parser.parse_args()

if not os.path.exists(options.inputDir): raise OSError(options.inputDir +' does not exist. This is where the input Root files go.')
if not os.path.exists(options.submitDir):
  print '== Making folder '+options.submitDir+' =='
  os.makedirs(options.submitDir)
if not os.path.exists(options.config): raise OSError(options.config +' does not exist. This is the config file for reading in the Root file.')

identifier = options.identifier

import pdb

import ROOT as r

import json
with open(options.config) as f:
  branchdata = json.load(f)

def findBranch(blist,bname):
  if bname not in blist:
    print '== \''+bname+'\' branch does not exist =='  
    return False
  else:
    print '== \''+bname+'\' branch is being read =='
    return True

def readRoot():
  from sys import stdout,argv
  from math import fabs
  finalmu = options.identifier 

  import glob
  filenamebase = options.inputDir 
  filenames = glob.glob(filenamebase+'/*.root')
  if len(filenames) == 0: raise OSError('Can\'t find any Root files to read in '+options.inputDir) 
  tree = r.TChain('oTree')
  for filename in filenames:
    statinfo = os.stat(filename)
    #if statinfo.st_size < 10000: continue #sometimes batch jobs fail
    print '== Reading in '+filename+' =='
    tree.Add(filename) 

  # make sure the branches are compatible between the two
  branches = set(i.GetName() for i in tree.GetListOfBranches())
  # required:
  if not findBranch(branches,options.jetpt): raise RuntimeError('Need jet pTs.')
  if not findBranch(branches,options.tjetpt): raise RuntimeError('Need truth jet pTs.')
  # optional:
  has_branches = []
  for branch in branchdata:
    branch['exists'] = findBranch(branches,branch['branchname'])
  nentries = tree.GetEntries()

  truepts = np.array([])
  recopts = np.array([])
  isLeadings = np.array([],dtype=bool)

  for branch in branchdata: branch['data'] = np.array([])

  for jentry in xrange(nentries):
    if jentry>options.numEvents and options.numEvents>0: continue
    tree.GetEntry(jentry)
    
    if not jentry%1000:
      stdout.write('== \r%d events read ==\n'%jentry)
      stdout.flush()

    jpts = np.array(getattr(tree,options.jetpt))
    tjpts = np.array(getattr(tree,options.tjetpt))
    if not len(jpts)==len(tjpts):
      raise RuntimeError('There should be the same number of reco jets as truth jets')

    goodinds = np.ones_like(tjpts) #which ones to keep
    for branch in branchdata:
      branch['tempdata'] = 0 #holds jetvals before appending to data
      if not branch['exists']: continue 
      if branch['type']=='event':
        val = getattr(tree,branch['branchname'])
        jetvals = np.ones_like(tjpts)*val
      if branch['type']=='jet':
        jetvals = np.array(getattr(tree,branch['branchname']))
        if not len(jetvals)==len(tjpts): raise RuntimeError(branch['branchname']+' should be the same length as truth jets')
      if 'min' in branch: goodinds = np.all([goodinds,jetvals>branch['min']],axis=0)
      if 'max' in branch: goodinds = np.all([goodinds,jetvals<branch['max']],axis=0)
      branch['tempdata'] = jetvals

    if len(tjpts)>2:
      jet3pt = -np.sort(-tjpts)[2]
      isLeading = tjpts>jet3pt
    else:
      isLeading = np.ones_like(tjpts,dtype=bool)

    truepts = np.append(truepts,tjpts[goodinds])
    recopts = np.append(recopts,jpts[goodinds])
    for branch in branchdata:
      if not branch['exists']: continue 
      branch['data'] = np.append(branch['data'],branch['tempdata'][goodinds])
    isLeadings = np.append(isLeadings,isLeading[goodinds])

  #end loop over entries
  save(options.submitDir+'/truepts_'+finalmu,truepts)
  save(options.submitDir+'/recopts_'+finalmu,recopts)
  save(options.submitDir+'/isLeading_'+finalmu,isLeadings)
  for branch in branchdata:
    if branch['exists']: save(options.submitDir+'/'+branch['name']+'_'+finalmu,branch['data'])
  return

readRoot()
#comment
