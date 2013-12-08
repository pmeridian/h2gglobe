#!/bin/tcsh

eos ls $1 | grep biasStudy | grep root | xargs -I {} echo "scripts/draw_simplebias_plots.py -o $2 -i root://eoscms//$1/{} -t '7 TeV'" 
