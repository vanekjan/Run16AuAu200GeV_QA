#!/bin/bash

#copylist - Run16 Au+Au @ 200  GeV: startLine 1, maxLine 583172
#for VPDMB triggers start at line 4783
#good HFT runs start at line 166378 (RunId > 17062047)
#choose any prediefined number of lines from picoList_all_new.list
sed -n '1,100001 p' ./picoLists/picoList_all_new.list > picoList_test.list



path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )


echo executing submitPicoHFMaker.csh f0r picoList_test.list inside $path

csh starSubmit/submitPicoHFMaker.csh $path picoList_test.list
