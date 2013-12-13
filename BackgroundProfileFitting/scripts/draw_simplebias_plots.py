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
ranges["pullMu"]=(-0.1,0.1)
ranges["pullBkg"]=(-0.1,0.1)
ranges["errMu"]=(0,4)
ranges["errBkg"]=(0,1)

drawRanges={}
drawRanges["mu"]=(-2,2)
drawRanges["pullMu"]=(-2,2)
drawRanges["pullBkg"]=(-2,2)
drawRanges["errBkg"]=(0,0.5)
drawRanges["errMu"]=(0,4)

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
labelMap['pullMu']='#mu pull'
labelMap['pullBkg']='bkg in 1FWHM pull'
labelMap['mu']='#mu'
labelMap['errMu']='#sigma_{#mu} from fit'
labelMap['errBkg']='#sigma_{bkg} from fit'

def getMaximumBias(graph):
  maxBias=-999
  for i in range(0,graph.GetN()):
    x=r.Double(0)
    y=r.Double(0)
    graph.GetPoint(i,x,y)
    if (x>=115 and x<=135):
      if (abs(y)>maxBias):
        maxBias=abs(y)
  return maxBias
        
def drawHists():
  file = r.TFile.Open(options.infile)

  os.system('mkdir -p %s'%options.outdir)
  cats=[]
  truth_mods=[]
  test_mods=[]
  muInj=[]
  muConstr=[]
  plot_types=[]
  masses=[]


  # ---- Decode file content
  for key in file.GetListOfKeys():
    name = key.GetName()
    cat = name.split('_')[0]
    truth = name.split('_')[2]
    test = name.split('_')[4]
    mu = name.split('_')[6]
    muC = name.split('_')[8]
    mass = name.split('_')[10]
    plot = name.split('_')[11]
    if cat not in cats: cats.append(cat)
    if truth not in truth_mods: truth_mods.append(truth)
    if test not in test_mods: test_mods.append(test)
    if mu not in muInj: muInj.append(mu)
    if muC not in muConstr: muConstr.append(muC)
    if mass not in masses: masses.append(mass)
    if plot not in plot_types: plot_types.append(plot)

  print "CATEGORIES: "+str(cats)
  print "TRUTH MODELS:"+str(truth_mods)
  print "TEST MODELS: "+str(test_mods)
  print "MU INJ: "+str(muInj)
  print "MU CONSTRAINT: "+str(muConstr)
  print "MASS TEST POINTS: "+str(masses)
  print "PLOT TYPES: "+str(plot_types)

  # ---- Prepare dictionary
  graphs={}
  for cat in cats:
    for truth in truth_mods:
      for test in test_mods:
        for mu in muInj:
          for muC in muConstr:
            for plot in plot_types:
              graphs[(cat,truth,test,plot)]=drawGraph(cat,truth,test,mu,muC,plot,masses)



  # ---- Assemble the final plots
  c=r.TCanvas("c","c",800,600)
  for cat in cats:
    for truth in truth_mods:
      for mu in muInj:
        for muC in muConstr:
          for plt in plot_types:
            axis=r.TH2F("a","a",10,110,150,10,drawRanges[plt][0],drawRanges[plt][1])
            axis.GetXaxis().SetTitle("#gamma#gamma invariant mass")
            axis.GetYaxis().SetTitle(labelMap[plt])
            axis.Draw()
            text=r.TLatex()
            text.SetTextAlign(12)
            text.SetTextSize(0.04)
            if (not options.additionalText == ""):
              text.DrawTextNDC(0.71,0.86,options.additionalText)
            text.DrawTextNDC(0.71,0.82,str(cat)+' muInj '+str(mu))
            text.DrawTextNDC(0.71,0.78,'truth '+str(truth))
            leg=r.TLegend(0.12,0.68,0.4,0.88)
            leg.SetBorderSize(0)
            leg.SetFillColor(0)
            maxBias=0  
            if (plt.find("pull")>-1):
              maxBias=0.15
            if (plt.find("pullBkg")>-1):
              maxBias=0.2
            if (maxBias>0):
              box=r.TBox(110,-maxBias,150,maxBias)
              box.SetFillColor(r.kGray)
              box.SetFillStyle(3001)
              box.SetLineColor(r.kBlack)
              box.SetLineStyle(3)
              box.SetLineWidth(2)
              box.Draw()
              line=r.TLine(110,0,150,0)
              line.SetLineColor(r.kBlack)
              line.SetLineStyle(2)
              line.SetLineWidth(1)
              line.Draw()
              line1=r.TLine(110,maxBias,150,maxBias)
              line1.SetLineColor(r.kBlack)
              line1.SetLineStyle(2)
              line1.SetLineWidth(1)
              line1.Draw()
              line2=r.TLine(110,-maxBias,150,-maxBias)
              line2.SetLineColor(r.kBlack)
              line2.SetLineStyle(2)
              line2.SetLineWidth(1)
              line2.Draw()

            for test in test_mods:
              graphs[(cat,truth,test,plt)]['median'].SetMarkerStyle(20)
              graphs[(cat,truth,test,plt)]['median'].SetMarkerSize(0.9)
              graphs[(cat,truth,test,plt)]['median'].SetMarkerColor(colorMap[test])
              graphs[(cat,truth,test,plt)]['median'].SetLineColor(colorMap[test])
              graphs[(cat,truth,test,plt)]['median'].SetLineWidth(2)
              graphs[(cat,truth,test,plt)]['median'].Draw("SAMEPL")
              maxAbsV=getMaximumBias(graphs[(cat,truth,test,plt)]['median'])
              if (plt.find("pull")>-1):
                leg.AddEntry(graphs[(cat,truth,test,plt)]['median'],str(test)+": maxBias(%3.2f)"%(maxAbsV),'pl')
              else:
                leg.AddEntry(graphs[(cat,truth,test,plt)]['median'],str(test),'pl')
            leg.Draw()
            axis.Draw('AXISSAME')
            for format in [".png",".pdf"]:
              c.SaveAs("%s/cat%d_truth_%s_muInj_%s_muConstr_%s_%s_median%s"%(options.outdir,int(cat.lstrip('cat')),truth,str(mu),str(muC),plt,format))
          

def drawGraph(cat=0,truth="exp1",test="exp1",mu="0.0",muC="0.0",plot="pull",masses=[110,120,130,140,150]):
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

  graphDict['median'].SetName("cat%d_truth_%s_test_%s_muInj_%s_muConstr_%s_median"%(int(cat.lstrip('cat')),truth,str(mu),str(muC),test))
  graphDict['mean'].SetName("cat%d_truth_%s_test_%s_muInj_%s_muConstr_%s_mean"%(int(cat.lstrip('cat')),truth,str(mu),str(muC),test))
  graphDict['rmsband'].SetName("cat%d_truth_%s_test_%s_muInj_%s_muConstr_%s_rmsband"%(int(cat.lstrip('cat')),truth,str(mu),str(muC),test))
  graphDict['band68'].SetName("cat%d_truth_%s_test_%s_muInj_%s_muConstr_%s_band68"%(int(cat.lstrip('cat')),truth,str(mu),str(muC),test))
  graphDict['band95'].SetName("cat%d_truth_%s_test_%s_muInj_%s_muConstr_%s_band95"%(int(cat.lstrip('cat')),truth,str(mu),str(muC),test))

  i=0
  print "+++++++ PLOT_TYPE "+str(plot)
  for mass in masses:
    h=file.Get("cat%d_truth_%s_test_%s_muInj_%s_muConstr_%s_mass_%d_%s"%(int(cat.lstrip('cat')),truth,test,str(mu),str(muC),int(mass),plotType))
    print str(h.GetName())+"\t"+str(h.GetMean())+"\t"+str(h.GetRMS())
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
