#!/bin/bash
##usage sh remove_overlapgene.sh input output

while read line; do
	start=`echo "$line" | cut -f 4`;
	end=`echo "$line" | cut -f 5`;

	if [ $start -lt $end ]; then
		echo "$line";
		#sed -i 's/ /\t/g' $line;	
	fi

done < $1 #>> $2
#mv o file
