declare -a LAT=(50S_ 40S_ 30S_ 20S_ 10S_ 00N_ 10N_ 20N_ 30N_ 40N_)
declare -a LON=(130W 120W 110W 100W 090W 080W 070W 060W 050W 040W 030W 020W 010W 000E 010E 020E 030E 040E 050E 060E 070E 080E 090E 100E 110E 120E 130E 140E 150E 160E 170E)

for lat in "${LAT[@]}"
do
    for lon in "${LON[@]}"
    do
	#echo https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_treecover2000_$lat$lon.tif
	wget https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_treecover2000_$lat$lon.tif
	wget https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_loss_$lat$lon.tif
	wget https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_gain_$lat$lon.tif
	wget https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_lossyear_$lat$lon.tif
	wget https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_datamask_$lat$lon.tif
	wget https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_first_$lat$lon.tif
	wget https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_last_$lat$lon.tif
    done
done
