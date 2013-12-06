#!/usr/bin/env python


import os
import sys
import fnmatch

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-d","--datfile")
(options,args)=parser.parse_args()

import ROOT as r
r.gROOT.SetBatch(1)
r.gSystem.Load('lib/libBackgroundProfileFitting.so')
r.gROOT.ProcessLine('.L scripts/makeHistsSimpleBias.C+g')

from ROOT import makeHistsSimpleBias

def makeHists(cat=(0,0),options=dict()):
    chain=r.TChain("muTree")
    input = open(options['input'])
    for line in input.readlines():
        if line.startswith('#') or line=='' or line =='\n': continue
        print line.strip()
        if line.find('eos'):
            chain.Add("root://eoscms//"+line.strip())
        else:
            chain.Add(line.strip())
    histDraw=makeHistsSimpleBias(chain)
    histDraw.outFile=options['outfile']
    histDraw.mycat=int(cat[0])
    for val in options['plot'].split(','):
        infos=val.split('(')
        histDraw.plot_vec.push_back(infos[0])
        histProp=infos[1].strip(')').split(':')
        histDraw.nbins_vec.push_back(int(histProp[0]))
        histDraw.low_vec.push_back(float(histProp[1]))
        histDraw.high_vec.push_back(float(histProp[2]))
    histDraw.Loop()

# --main-- here
if not options.datfile:
  print "Please provide a -d datfile"
  exit(-1)
else:
  f = open(options.datfile)
  sw = r.TStopwatch()
  configDict={}
  for line in f.readlines():
    if line.startswith('#') or line=='' or line =='\n': continue
    info = line.strip().split('=')
    if (info[0]=='cat'):
      myOptions={}
      (cat,mu)=info[1].strip('(').strip(')').split(',')
      configDict[(cat,mu)]=myOptions
    else:
      myOptions[info[0]]=info[1]
#  print configDict
  sw.Reset()
  sw.Start()
  for cat in configDict.keys():
      makeHists(cat,configDict[cat])
  sw.Stop()
  print 'Took:', sw.Print()
  f.close()
    
