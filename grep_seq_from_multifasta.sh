#!/bin/bash
# usage: this_script.sh <file_list> <input_filename> <output_filename>

filename=$1

while read lines 
do
#sed -n -e "/$lines/,/>/p" $2 | sed -e "$ s/^>.*//" | sed -e "/^$/ d" >> $3 

sed -n -e "/$lines/,/>/p" $2 | sed -e "$ s/^>.*//" | sed -e "/^$/ d" >> $3

done<$filename