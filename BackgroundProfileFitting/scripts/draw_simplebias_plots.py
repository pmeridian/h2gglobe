#!/usr/bin/env python

import sys
import os
from array import *

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-o","--outdir")
parser.add_option("-i","--infile")
parser.add_option("-t","--additionalText")
parser.add_option("-L","--legLeft",default=False,action="store_true")

(options,args)=parser.parse_args()

import ROOT as r
r.gROOT.SetBatch(1)
r.gStyle.SetOptStat(0)
r.gStyle.SetOptTitle(0)

ranges={}
ranges["mu"]=(-5,5)
ranges["pull"]=(-0.1,0.1)
ranges["err"]=(0,4)

drawRanges={}
drawRanges["mu"]=(-2,2)
drawRanges["pull"]=(-2,2)
drawRanges["err"]=(0,4)

colorMap={}
colorMap['exp1']=r.kRed
colorMap['exp2']=r.kRed
colorMap['exp3']=r.kRed
colorMap['pow1']=r.kGreen+1
colorMap['pow2']=r.kGreen+1
colorMap['pow3']=r.kGreen+1
colorMap['lau1']=r.kBlue+1
colorMap['lau2']=r.kBlue+1
colorMap['lau3']=r.kBlue+1
colorMap['pol1']=r.kBlack
colorMap['pol2']=r.kMagenta
colorMap['pol3']=r.kCyan
colorMap['pol4']=r.kOrange
colorMap['pol5']=r.kMagenta-10
colorMap['pol6']=r.kRed-9
colorMap['pol7']=r.kGreen-7

labelMap={}
labelMap['pull']='#mu pull'
labelMap['mu']='#mu'
labelMap['err']='#sigma_{#mu} from fit'

def drawHists():
  file = r.TFile.Open(options.infile)

  os.system('mkdir -p %s'%options.outdir)
  cats=[]
  truth_mods=[]
  test_mods=[]
  muInj=[]
  plot_types=[]
  masses=[]


  # ---- Decode file content
  for key in file.GetListOfKeys():
    name = key.GetName()
    cat = name.split('_')[0]
    truth = name.split('_')[2]
    test = name.split('_')[4]
    mu = name.split('_')[6]
    mass = name.split('_')[8]
    plot = name.split('_')[9]
    if cat not in cats: cats.append(cat)
    if truth not in truth_mods: truth_mods.append(truth)
    if test not in test_mods: test_mods.append(test)
    if mu not in muInj: muInj.append(mu)
    if mass not in masses: masses.append(mass)
    if plot not in plot_types: plot_types.append(plot)

  print "CATEGORIES: "+str(cats)
  print "TRUTH MODELS:"+str(truth_mods)
  print "TEST MODELS: "+str(test_mods)
  print "MU INJ: "+str(muInj)
  print "MASS TEST POINTS: "+str(masses)
  print "PLOT TYPES: "+str(plot_types)

  # ---- Prepare dictionary
  graphs={}
  for cat in cats:
    for truth in truth_mods:
      for test in test_mods:
        for mu in muInj:
          for plot in plot_types:
            graphs[(cat,truth,test,plot)]=drawGraph(cat,truth,test,mu,plot,masses)



  # ---- Assemble the final plots
  c=r.TCanvas("c","c",800,600)
  for cat in cats:
    for truth in truth_mods:
      for mu in muInj:
        for plt in plot_types:
          axis=r.TH2F("a","a",10,110,150,10,drawRanges[plt][0],drawRanges[plt][1])
          axis.GetXaxis().SetTitle("#gamma#gamma invariant mass")
          axis.GetYaxis().SetTitle(labelMap[plt])
          axis.Draw()
          text=r.TLatex()
          text.SetTextAlign(12)
          text.SetTextSize(0.04)
          text.DrawTextNDC(0.73,0.86,options.additionalText)
          text.DrawTextNDC(0.73,0.82,str(cat))
          text.DrawTextNDC(0.73,0.78,'truth '+str(truth))
          leg=r.TLegend(0.12,0.68,0.4,0.88)
          leg.SetBorderSize(0)
          leg.SetFillColor(0)
          if (plt=="pull"):
            box=r.TBox(110,-0.15,150,0.15)
            box.SetFillColor(r.kGray)
            box.SetFillStyle(3001)
            box.SetLineColor(r.kGray+2)
            box.SetLineStyle(2)
            box.SetLineWidth(2)
            box.Draw()
          for test in test_mods:
            graphs[(cat,truth,test,plt)]['median'].SetMarkerStyle(20)
            graphs[(cat,truth,test,plt)]['median'].SetMarkerSize(0.9)
            graphs[(cat,truth,test,plt)]['median'].SetMarkerColor(colorMap[test])
            graphs[(cat,truth,test,plt)]['median'].SetLineColor(colorMap[test])
            graphs[(cat,truth,test,plt)]['median'].SetLineWidth(2)
            graphs[(cat,truth,test,plt)]['median'].Draw("SAMEPL")
          leg.AddEntry(graphs[(cat,truth,test,plt)]['median'],str(test),'pl')
          leg.Draw()
          axis.Draw('AXISSAME')
          for format in [".png",".pdf"]:
            c.SaveAs("%s/cat%d_truth_%s_muInj_%s_%s_median%s"%(options.outdir,int(cat.lstrip('cat')),truth,srt(mu),plt,format))
          

def drawGraph(cat=0,truth="exp1",test="exp1",mu="0.0",plot="pull",masses=[110,120,130,140,150]):
  file = r.TFile.Open(options.infile)
  histos={}
  plotType=plot
  if (plot=="pull"):
    plotType="pull_mu"
  elif (plot=="err"):
    plotType="err_mu"

  graphDict={}
  graphDict['median']=r.TGraph()
  graphDict['mean']=r.TGraph()
  graphDict['rmsband']=r.TGraphErrors()
  graphDict['band68']=r.TGraphAsymmErrors()
  graphDict['band95']=r.TGraphAsymmErrors()

  graphDict['median'].SetName("cat%d_truth_%s_test_%s_muInj_%s_median"%(int(cat.lstrip('cat')),truth,str(mu),test))
  graphDict['mean'].SetName("cat%d_truth_%s_test_%s_muInj_%s_mean"%(int(cat.lstrip('cat')),truth,str(mu),test))
  graphDict['rmsband'].SetName("cat%d_truth_%s_test_%s_muInj_%s_rmsband"%(int(cat.lstrip('cat')),truth,str(mu).test))
  graphDict['band68'].SetName("cat%d_truth_%s_test_%s_muInj_%s_band68"%(int(cat.lstrip('cat')),truth,str(mu),test))
  graphDict['band95'].SetName("cat%d_truth_%s_test_%s_muInj_%s_band95"%(int(cat.lstrip('cat')),truth,str(mu),test))

  i=0
  print "+++++++ PLOT_TYPE "+str(plot)
  for mass in masses:
    h=file.Get("cat%d_truth_%s_test_%s_muInj_%s_mass_%d_%s"%(int(cat.lstrip('cat')),truth,test,int(mass),str(mu),plotType))
    print h.GetName()
    h.GetXaxis().SetRangeUser(ranges[plot][0],ranges[plot][1])
    quantiles=array('d',[0,0,0,0,0])
    probs=array('d',[0.025,(1.-0.683)/2.,0.5,1-(1.-0.683)/2.,0.975])
    h.GetQuantiles(5,quantiles,probs)
    graphDict['mean'].SetPoint(i,float(mass),float(h.GetMean()))
    graphDict['median'].SetPoint(i,float(mass),float(quantiles[2]))
    graphDict['rmsband'].SetPoint(i,float(mass),float(quantiles[2]))
    graphDict['rmsband'].SetPointError(i,0,float(h.GetRMS()))
    graphDict['band68'].SetPoint(i,float(mass),float(quantiles[2]))
    graphDict['band68'].SetPointEYlow(i,float(quantiles[2]-quantiles[1]))
    graphDict['band68'].SetPointEYhigh(i,float(quantiles[3]-quantiles[2]))
    graphDict['band95'].SetPoint(i,float(mass),float(quantiles[2]))
    graphDict['band95'].SetPointEYlow(i,float(quantiles[2]-quantiles[0]))
    graphDict['band95'].SetPointEYhigh(i,float(quantiles[4]-quantiles[2]))
    i=i+1

  return graphDict

drawHists()
