#!/usr/bin/env python

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-D","--readFromDat",dest="readFromDat",type="str",help="Read these run options from datfile")
parser.add_option("","--executableName",dest="executableName",type="str",default="StandardBiasStudy",help="Executable name")
parser.add_option("-s","--sigfilename",dest="sigfilename",type="str",help="Input signal workspace file")
parser.add_option("-b","--bkgfilename",dest="bkgfilename",type="str",help="Input background/data workspace file")
parser.add_option("-w","--sigwsname",dest="sigwsname",type="str",help="Input signal workspace name",default="cms_hgg_workspace")
parser.add_option("-B","--bkgwsname",dest="bkgwsname",type="str",help="Input background/data workspace name",default="cms_hgg_workspace")
parser.add_option("-d","--datfile",dest="datfile",type="str",help="Config datfile")
parser.add_option("-o","--outerDir",dest="outerDir",help="Dir name")
parser.add_option("-c","--cats",dest="cats",default=[],action="append")
parser.add_option("-m","--expectSignals",default=[],action="append")
parser.add_option("-M","--expectSignalMasses",default=[],action="append")
parser.add_option("-L","--mulow",dest="mulow",type="float",help="Value of mu to start scan")
parser.add_option("-H","--muhigh",dest="muhigh",type="float",help="Value of mu to end scan")
parser.add_option("-S","--mustep",dest="mustep",type="float",help="Value of mu step size")
parser.add_option("-C","--constraintMuWidth",dest="constrainMuWidth",type="float",help="sigma of mu constraint")
parser.add_option("-t","--toysperjob",dest="toysperjob",type="int",help="Number of toys to run per job")
parser.add_option("-n","--njobs",dest="njobs",type="int",help="Number of jobs to run")
parser.add_option("-j","--jstart",dest="jstart",type="int",default=0,help="Start job number here")
parser.add_option("-q","--queue",dest="queue",type="str",default="8nh",help="Which queue to run job in")
parser.add_option("-e","--eosPath",dest="eosPath",help="Write output files to eos")
parser.add_option("--takeOtherFiles",type="str",help="Copy these files over to batch sandbox")
parser.add_option("","--skipPlots",default=False,action="store_true",help="Don\'t plot all the envelopes")
parser.add_option("","--bkgOnly",default=False,action="store_true",help="Background only toys (mu forced to 0)")
parser.add_option("","--useCondor",dest="useCondor",default=False,action="store_true",help="Submit locally using condor")
parser.add_option("","--copyWorkspace",dest="copyWorkspace",default=False,action="store_true",help="Copy the workspaces into the running sandbox")
parser.add_option("","--dryRun",dest="dryRun",default=False,action="store_true",help="Don't submit")
(options,args)=parser.parse_args()

import os
if options.eosPath:
  if '/eos/cms' in options.eosPath:
    options.eosPath = options.eosPath.split('/eos/cms')[1]
  os.system('cmsMkdir %s'%options.eosPath)

def writeSubScript(cat,mlow,mhigh,outdir,muInject,massInject,constraintMu=0,constraintMuWidth=0):
  cat = int(cat)
  muInject=float(muInject)
  massInject=float(massInject)
  constrainMu=int(constraintMu)
  constrainMuWidth=float(constraintMuWidth)

  subline = './%s -s %s -b %s --sigwsname %s --bkgwsname %s -d %s -c %d -L %3.1f -H %3.1f -t %d -D %s --expectSignal=%3.1f --expectSignalMass=%4.1f'%(options.executableName, os.path.basename(options.sigfilename) if options.copyWorkspace else options.sigfilename,os.path.basename(options.bkgfilename) if options.copyWorkspace else options.bkgfilename,options.sigwsname,options.bkgwsname,os.path.basename(options.datfile),cat,mlow,mhigh,options.toysperjob,os.path.abspath(outdir),muInject,massInject)

  if options.skipPlots: subline += ' --skipPlots'
  if options.bkgOnly: subline += ' --bkgOnly'
  if (constrainMu>0 or constrainMuWidth>0):
    subline+= ' --constraintMu --constraintMuWidth '+str(constrainMuWidth)
    
  for j in range(options.jstart,options.njobs):
    f = open('%s/%s/sub_cat%d_job%d.sh'%(os.getcwd(),outdir,cat,j),'w')
    f.write('#!/bin/bash\n')
    if (options.useCondor):
      f.write('cd `mktemp -d`\n')
    f.write('rm -f %s.done\n'%(f.name))
    f.write('rm -f %s.fail\n'%(f.name))
    f.write('rm -f %s.run\n'%(f.name))
    f.write('rm -f %s.log\n'%(f.name))
    f.write('cd %s\n'%os.getcwd())
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('cd -\n')
    if (options.copyWorkspace):
      f.write('cp %s .\n'%os.path.abspath(options.sigfilename))
      f.write('cp %s .\n'%os.path.abspath(options.bkgfilename))
    f.write('cp %s .\n'%os.path.abspath(options.datfile))
    f.write('cp %s/bin/%s .\n'%(os.getcwd(),options.executableName))
    if options.takeOtherFiles:
     for sandbox_file in options.takeOtherFiles.split(','):
       f.write('cp %s . \n'%os.path.abspath(sandbox_file))
    f.write('touch %s.run\n'%(f.name))
    execline = subline + ' -j %d -o StandardBiasStudyOut_cat%d_job%d.root'%(j,cat,j)
    f.write('if ( %s ) then \n'%execline)
    if options.eosPath:
      f.write('\tif ( cmsStage StandardBiasStudyOut_cat%d_job%d.root %s/%s/StandardBiasStudyOut_cat%d_job%d.root ) then \n'%(cat,j,options.eosPath,outdir,cat,j))
    else:
      f.write('\tif ( cp StandardBiasStudyOut_cat%d_job%d.root %s/%s ) then\n'%(cat,j,os.getcwd(),outdir))
    f.write('\t\ttouch %s.done\n'%(f.name))
    f.write('\t\trm -f %s.run\n'%(f.name))
    f.write('\telse\n')
    f.write('\t\ttouch %s.fail\n'%(f.name))
    f.write('\tfi\n')
    f.write('else\n')
    f.write('\ttouch %s.fail\n'%(f.name))
    f.write('fi\n')
    if (options.useCondor):
      f.write('rm -rf *\n')
    f.close()
    os.system('chmod +x %s'%f.name)
    if not options.dryRun:
      if (not options.useCondor):

        os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))
      else:
        fc = open('%s/%s/sub_cat%d_job%d_condor'%(os.getcwd(),outdir,cat,j),'w')
        fc.write('universe = vanilla\n')
        fc.write('executable = %s\n'%(f.name))
        fc.write('output = %s.log\n'%(f.name))
        fc.write('error = %s.log\n'%(f.name))
        fc.write('log = %s.log\n'%(f.name))
        fc.write('queue 1\n')
        os.system('condor_submit %s'%(fc.name))
    else:
      print 'bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name)

if not options.readFromDat:
  for cat in options.cats:
    cat = int(cat)
    for i, mu in enumerate(options.expectSignals):
      mu=float(mu)
      storage_dir = '%s/cat%d_mu%3.1f'%(options.outerDir,cat,mu)
      if options.eosPath:
        os.system('cmsMkdir %s/%s'%(options.eosPath,storage_dir))
      os.system('mkdir -p %s'%storage_dir)
      writeSubScript(cat,options.mulow,options.muhigh,options.mustep,storage_dir,mu,options.expectSignalMasses[i])

else:
  f = open(options.readFromDat)
  configDict={}
  for line in f.readlines():
    if line.startswith('#') or line.startswith(' ') or line.startswith('\n'): continue
    if not line.startswith('cat'):
      #configDict[line.split('=')[0]] = line.split('=')[1].strip('\n')
      setattr(options,line.split('=')[0],line.split('=')[1].strip('\n'))
    else:
      lineConfig = line.split()
      cat = int(lineConfig[0].split('=')[1])
      injectmu = float(lineConfig[1].split('=')[1])
      mlow = float(lineConfig[2].split('=')[1])
      mhigh = float(lineConfig[3].split('=')[1])
#      mstep = float(lineConfig[4].split('=')[1])
      injectmass = float(lineConfig[4].split('=')[1])
      if (len(lineConfig)>5):
        constrainMu = int(lineConfig[5].split('=')[1])
        constrainMuWidth = float(lineConfig[6].split('=')[1])
        storage_dir = '%s/cat%d_mu%3.1f_constrainMu%3.2f_mass%d'%(options.outerDir,cat,injectmu,constrainMuWidth,injectmass)
      elif (len(lineConfig)>3):
        constrainMu = -1
        constrainMuWidth = -1
        storage_dir = '%s/cat%d_mu%3.1f_mass%d'%(options.outerDir,cat,injectmu,injectmass)
      else:
        storage_dir = '%s/cat%d_mu%3.1f_mass%d'%(options.outerDir,cat,injectmu,injectmass)
      if options.eosPath:
        os.system('cmsMkdir %s/%s'%(options.eosPath,storage_dir))
      os.system('mkdir -p %s'%storage_dir)
      writeSubScript(cat,mlow,mhigh,storage_dir,injectmu,injectmass,constrainMu,constrainMuWidth)

