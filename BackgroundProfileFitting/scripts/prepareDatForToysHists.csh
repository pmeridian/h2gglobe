#!/bin/tcsh

#set outdir=root://eoscms//eos/cms/store/group/phys_higgs/meridian/BiasStudy/legacy_freeze_v2_massFac_d
set outdir=root://eoscms//eos/cms/store/group/phys_higgs/meridian/BiasStudy/legacy_freeze_v3_massFac_7TeV_test
set outname=biasStudy_legacy_freeze_v3_massFac_7TeV_test

foreach file ( $1/*txt )
    set base=`basename $file`
    set cat=`echo $base | awk -F '_mu'  '{print $1}' | sed -e "s%cat%%g"`
    set mu=`echo $base | awk -F '_mu'  '{print $2}' | awk -F '_' '{print $1}'`
    set muConstrain=`echo $base | awk -F '_constrainMu'  '{print $2}' | awk -F '_' '{print $1}'`
    echo "#######  cat ${cat} with muInj ${mu} and muConstrain ${muConstrain} ######"
    echo $file | awk '{printf "cat=('$cat','${mu}','${muConstrain}')\nplot=mu(200:-10:10),errMu(100:0:5),errBkg(100:0:2),pullMu(200:-2.5:2.5),pullBkg(200:-2.5:2.5)\ninput='$file'\noutfile='${outdir}'/'${outname}'_'$base:r'.root\n"}'
end
