#########################################################################
# File Name: getRunnumber.sh
# Created Time: Thu 26 Mar 2015 05:38:55 PM EDT
#########################################################################
#!/bin/csh
rm Runnumber.list
get_file_list.pl -keys 'orda ( runnumber )' -cond production=P16ij,trgsetupname=AuAu200_production2_2016,filetype=daq_reco_picoDst,filename~sst_,storage=local -limit 0 -delim -distinct > Runnumber.list
