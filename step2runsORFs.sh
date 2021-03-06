#!/bin/bash

#########################################################
# Author: shuhailams					#
# Date created: 28 feb 2014				#
#							#
# create dirlist by using cmd: ls -d */ > dirlist 	#
#							#
#########################################################

#copy directory name to list
FILELIST="dirlist.txt"

mkdir sORF_cluster

while read lines
do	

#enter directory	
cd $lines

#combine scaffolds/chromosomes
	cat *.fna > "$lines"_genomic.fa
	cat *.faa > "$lines"_protein.faa
	cat *_intron.fasta > "$lines"_intron.fasta
	cat *.gff > "$lines"_genomic.gff
	cat *.ffn > "$lines"_coding.fasta
	

#sORF from genome annotation
	grep ">" *_protein.faa >"$lines"_orfp
	infoseq -name -length -only *_protein.faa > "$lines"_lengthp
	awk '$2<=80 {print $1}' "$lines"_lengthp > "$lines"_listp
	grep -Ff "$lines"_listp "$lines"_orfp > "$lines"_sorflistp
	perl ../fasta_getseq.pl "$lines"_sorflistp *_protein.faa > "$lines"_sorfp_genomeannotation


##sORF from ab initio
	##extract igr
		##index genome sequences
		samtools faidx "$lines"_genomic.fa
		#cek genome index success or not
		cat "$lines"_genomic.fa.fai
				
		#python ../get_intergenic.py *_genomic.gff
		##grep gene type from genomic.gff3
		perl ../grep_type_gene.pl "$lines"_genomic.gff > "$lines"_gene.gff
	
		perl ../myigrgff.pl "$lines"_gene.gff "$lines"_genomic.fa.fai > "$lines"_igr.gff
		
		##remove igr with start higher than stop.. 
		sh ../remove_overlapgene.sh "$lines"_igr.gff >"$lines"_igr_removerlapgene.gff
	
		##extract igr sequences
		bedtools getfasta -fi "$lines"_genomic.fa -bed "$lines"_igr_removerlapgene.gff -fo "$lines"_igr.fa
		
	##predict sORF by getorf
		getorf -sequence "$lines"_igr.fa -outseq "$lines"_igr_getorf.fasta -find 3 -maxsize 240
		getorf -sequence "$lines"_genomic.fa -outseq "$lines"_genomic_getorf.fasta -find 3 -maxsize 240
		cat "$lines"_igr_getorf.fasta "$lines"_genomic_getorf.fasta >"$lines"_getorf.fasta
				
	##predict sORF by sorffinder
		##makemodel	
		perl /pathtosORFfinder/src/make_model.pl -c "$lines"_coding.fasta -n "$lines"_intron.fasta -o "$lines"_matrix
		##simulate
		perl /pathtosORFfinder/src/simulate.pl -m "$lines"_matrix -p 0.5 -o "$lines"_simulate0.5
		##sorffinder
		perl /pathtosORFfinder/src/search_sORF.pl -m "$lines"_matrix -p 0.5 -s "$lines"_simulate0.5 -i "$lines"_igr.fa -o "$lines"_sORF0.5.fasta -d b
			
	## filter length sORFfinder
		##remove header result sORFfinder
		sed '1,3d' "$lines"_sORF0.5.fasta >> "$lines"_sORF0.5_rename.fasta
	
		##multifasta to single fasta
		perl -pe '/^>/ ? print "\n" : chomp' "$lines"_sORF0.5_rename.fasta > "$lines"_sORF0.5_oneline.fasta
	
		##filter fasta by length <240nt
		awk '!/^>/ { next } { getline seq } length(seq) <= 240 { print $0 "\n" seq }' "$lines"_sORF0.5_oneline.fasta > "$lines"_sORF0.5_filter.fasta
	
	##rename fasta header according to software used
		##rename fasta header for sorffinder
		awk '/>/{print ">sf"(++i)}!/>/' < "$lines"_sORF0.5_filter.fasta > "$lines"_sORF0.5_filter_rename.fasta
	
		##rename fasta header for getorf-genome
		awk '/>/{print ">getge"(++i)}!/>/' < "$lines"_genome_getorf.fasta >"$lines"_genome_getorf_rename.fasta
	
		##rename fasta header for getorf-igr
		awk '/>/{print ">getigr"(++i)}!/>/' < "$lines"_igr_getorf.fasta >"$lines"_igr_getorf_rename.fasta
	
		echo "rename completed"
	
		##combine all sorf prediction
		cat "$lines"_sORF0.5_filter_rename.fasta "$lines"_genome_getorf_rename.fasta "$lines"_igr_getorf_rename.fasta >"$lines"_sORF_combined.fasta
	
		##cluster based on 100% using cd-hit-est 
		##remove redundant
		cd-hit-est -i "$lines"_sORF_combined.fasta -o "$lines"_sORF_combined_cdhit100 -c 1.0 -n 8
	
		##grab sORF that agree by genome_getorf, intergenic_getorf and sorffinder
		python ../grab2.py "$lines"_sORF_combined_cdhit100.clstr >grab_list
	
		##rename * to sorf
		sed 's/>*... \*/ sorf/g' grab_list >grab_list2
	
		## grep header with sorf
		grep "sorf" grab_list2 >grab_list_sorf
	
		#edit grab_list_sorf
		sed 's/ /\t/' grab_list_sorf >grab_list_sorf_edit
	
		cut -f3 grab_list_sorf_edit > list_sorf
		
		##grep sequences sORFs
		perl ../fasta_getseq.pl list_sorf "$lines"_sORF_combined_cdhit100 >"$lines"_sORF_predicted.fasta

		##blast sORFs against coding sequences
		makeblastdb -in *_protein.faa -dbtype prot
		blastx -db *_protein.faa -query "$lines"_sORF_predicted.fasta -out "$lines"_sORF_abinitiopredictedVSprotein -outfmt 0 -num_threads 20
		perl ../bpSearchIO.pl "$lines"_sORF_abinitiopredictedVSprotein "$lines"_sORF_abinitiopredictedVSprotein.parsed

		##grep sorf seq hit protein
		cut -f1 "$lines"_sORF_abinitiopredictedVSprotein.parsed | sort -u | sed 's/\_1//g'> abinitio_sorf_hits_protein
		
		#grep list id from sORF predicted.fasta (sorf-cluster.sh)
		grep ">" "$lines"_sORF_predicted.fasta > "$lines"_sORF_predicted_list 
		sed 's/>//g' "$lines"_sORF_predicted_list > "$lines"_sORF_predicted_list2
		
		#grep sORF do no have mapped to ORF
		grep -v -w -f abinitio_sorf_hits_protein "$lines"_sORF_predicted_list2 > abinitio_sorf_final
	
		##grep sORFs predicted final
		perl ../fasta_getseq.pl abinitio_sorf_final "$lines"_sORF_combined_cdhit100 > "$lines"_sORF_abinitio_final.fasta

		##translate sORFs ab initio
		transeq "$lines"_sORF_abinitio_final.fasta > "$lines"_sORF_abinitio_final.faa

##combine sORFs genome annotation and ab initio
		##combine sORF abinitio + genome annotation
	cat "$lines"_sORF_abinitio_final.faa "$lines"_sorfp_genomeannotation > "$lines"_sorf_final
	
	awk '/>/{sub(">","&"FILENAME"-");sub(/\.fasta/,x)}1' "$lines"_sorf_final > "$lines"_sorf_final2
	##rename with species name
	for f in "$lines"_sorf_final2
	do 
    	    		
	#remove _sorf_final
	sed -i "s/\_sorf_final2//g" "$f"
		
	done

	cp "$lines"_sorf_final2 ../sORF_cluster

	echo "$lines"

#cd into upper directory
cd ..

done < $FILELIST

##clustering to get conserved sORFs
#enter folder sORF_cluster
cd sORF_cluster

	##combined sORFs from all fungal genomes
		cat *_sorf_final2 > sORF_combined.fasta
	##cluster sORF combined
		cd-hit -i sORF_combined.fasta -o cluster_sORF_combined_cdhit70 -c 0.7 -n 5 -g 1 -G 0 -aS 0.8 -d 500 -p 1 > db_70.log
		echo "cdhit"
	##sh count_organism_uniq.sh and need to change input.clstr name first
	
		csplit -z cluster_sORF_combined_cdhit70.clstr '/Cluster /' {*}
		mkdir cluster_sORF_combined_cdhit70_ori
		cp xx* cluster_sORF_combined_cdhit70_ori
	
	##count how many species in each clustered
		for i in xx*
		do
			j=`more $i | wc -l` 
			j2=$(( j - 1 ))
		
			 echo $i": "$j2 >> count.xx
			 sed -i 's/\-.*//g' $i
			##rename Aspergillus_fumigatus-sf488 nk tukar jd A.fumagitus-sf488
			 sed -i 's/\-.*//g' $i
			organism=`awk '{ print $3 }' $i`	
		
			 echo $organism > $i.organism
			 sed -i 's/ >/\n>/g' $i.organism
		done
		
		for i in *.organism
		do
			sort $i | uniq -c | sort -n > $i.uniq
		done
		
		for i in *.uniq
		do
			k=`more $i | wc -l` 
			echo $i $k >> uniq.xx
		done
		echo "species count complete"

	##output count.xx uniq.xx
	##Select atleast two species in one cluster
	##cut -f2 "2" uniq.xx > listuniq_70
	##cut -f2 | awk '>2 {print}'
	awk '{ if ($2 >= 2) { print } }' uniq.xx > listuniq_70
	echo "select conserved sORFs in clusters"

	##edit cluster name by remove " .organism.uniq "
	sed -i "s/\..*//" listuniq_70

##create file for selected cluster	
	mkdir cluster_listuniq_70
	cd cluster_sORF_combined_cdhit70_ori
	
		#grep cluster based on conserved cluster
		xargs -a ../listuniq_70 cp -t ../cluster_listuniq_70
		
		#count files in folder
		find . -type f | wc -l
	cd ..

#get into cluster_listuniq_70 to grep fasta seq
	cd cluster_listuniq_70
	#cat all sorf cluster
	cat xx* >cluster_sorf
	
	#edit to get gene_id only from list of cluster
	sed 's\>Cluster.*\\g' cluster_sorf >test
	grep ">" test > test1
	sed 's/\...*//g' test1 >test2
	sed 's/, >/\t\>/g' test2 >test3
	cut -f3 test3 >test4
	
	#grep seq
	sh ../../grep_seq_from_multifasta.sh test4 ../sORF_combined.fasta conserved_sorf.fasta
	echo "sorf conserved completed"

	cd ..

#out folder sORF_cluster
cd ..
