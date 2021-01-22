#!/bin/csh
set nRunnumbers=`cat Runnumber.list`
	set i=1
while ( $i <= $#nRunnumbers )
	set iRunnumber = `echo $nRunnumbers[$i]`
	
	set iDay = `echo $iRunnumber | awk '{print substr($0,3,3)}'`
	echo $iDay $iRunnumber
	
	if ( $i == 1 ) then
		get_file_list.pl -keys path,filename -cond production=P16ij,trgsetupname=AuAu_200_production_2016,filetype=daq_reco_picoDst,filename~st_physics,storage=local,runnumber=$iRunnumber -limit 0 -delim / -distinct > ./picoLists/pico.list
	else
		get_file_list.pl -keys path,filename -cond production=P16ij,trgsetupname=AuAu_200_production_2016,filetype=daq_reco_picoDst,filename~st_physics,storage=local,runnumber=$iRunnumber -limit 0 -delim / -distinct >> ./picoLists/pico.list
	endif
	
	@ i = $i + 1
	
end
