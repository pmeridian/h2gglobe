#!/bin/tcsh

mkdir -p $2

foreach dir (`/afs/cern.ch/project/eos/installation/0.3.1-22/bin/eos.select find -d $1 | sed -e "s%$1%%g" | awk '{print $1}' | awk -F '_mass' '{print $1}' | grep cat | sort | uniq`)
    echo "+++++++ Preparing list for toys $dir"
    /afs/cern.ch/project/eos/installation/0.3.1-22/bin/eos.select find -f $1 | grep root | grep $dir >! $2/${dir}_list.txt
    echo "+++++++ Found `wc -l  $2/${dir}_list.txt` files for toys $dir. List is @  $2/${dir}_list.txt"
end
echo "=========== Summary =========="
wc -l $2/*_list.txt
