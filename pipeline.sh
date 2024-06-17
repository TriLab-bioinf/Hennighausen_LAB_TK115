## retrieve CDS from gbk file
python retrieve_sequences.py Wikim_164_afdsb_consensus.pgap.gbk Wikim_164_afdsb_consensus.pgap.CDS.fasta
python retrieve_sequences.py Wikim_166_afdsc_consensus.pgap.gbk Wikim_166_afdsc_consensus.pgap.CDS.fasta

## download protein sequences from NCBI
module load edirect
cat housekeeping.genes.filter |while read lines
do
echo -e "esearch -db protein -query "$lines" | efetch -format fasta > housekeeping/${lines}.AA.fasta" >>esearch.swarm 
done

swarm -f esearch.swarm -t 2 --module edirect --time 24:00:00

## Calculate the sequences similarity between known genomes and two assemblies
cat housekeeping.genes.filter|while read lines
do
	## retrieve specific protein sequences 
	grep -A 1 $lines Wikim_164_afdsb_consensus.pgap.CDS.fasta >blast/Wikim_164_${lines}.faa
	grep -A 1 $lines Wikim_166_afdsc_consensus.pgap.CDS.fasta >blast/Wikim_166_${lines}.faa
	## build database for blastp
	makeblastdb -in blast/Wikim_164_${lines}.faa -dbtype prot -title blast/Wikim_164_${lines}
	makeblastdb -in blast/Wikim_166_${lines}.faa -dbtype prot -title blast/Wikim_166_${lines}
	## alignment
	blastp -db blast/Wikim_164_${lines}.faa -query housekeeping/${lines}.AA.fasta -outfmt 6 -out blast/Wikim_164_${lines}.res.txt
	blastp -db blast/Wikim_166_${lines}.faa -query housekeeping/${lines}.AA.fasta -outfmt 6 -out blast/Wikim_166_${lines}.res.txt
done

## filter sequences that are 100% matched with two assemblies
for i in `ls blast/*res.txt`
do
	file=`basename $i .res.txt`
	awk '{if($3==100 && $11==0) print}' $i |sort -k3,3gr
done>blast.res.txt

## identify the filtered sequences species
cut -f1 blast.res.txt|sort -u|grep -f - -w housekeeping/*fasta >species.txt

## unique the species
sed 's/^.*\[//g;s/\]//g' species.txt|sort -u >species.filter.txt

## genome comparison
for i in `ls genomes/Weissella/*fna`
do
	file=`basename $i .fna`
	echo -e "nucmer --mum Wikim_164_afdsb_consensus.fasta -p delta/Wikim_164_to_${file} $i" >>nucmer.swarm
	echo -e "nucmer --mum Wikim_166_afdsc_consensus.fasta -p delta/Wikim_166_to_${file} $i" >>nucmer.swarm
done

swarm -f nucmer.swarm -t 36 --module mummer --sbatch "--mem=80g"

for i in `ls delta/*delta`
do
	file=`basename $i .delta`
	## A dotplot between two sequences to reveal their macroscopic similarity
	mummerplot -l $i -f --large -prefix plots/${file} -t png --color
	## identify all the SNPs and indels between the two sequence sets
	show-snps -C $i -T -r > SNP/${file}.snps
done

## filter SNPs
for i in `ls SNP/*snps`
do
	n=`wc -l $i|cut -d ' ' -f1`
	if [ $n -gt 1000 -a $n -lt 30000 ]
	then
		cp $i SNP.filter
	fi
done

## combine the filter SNP files to a unique SNP list
for i in `ls SNP.filter/Wikim_166*snps`
do 
	grep "AFDS" $i|grep -v "TK115"|cut -f1,2,9 
done>Wikim_166.tmp

sort -u Wikim_166.tmp|sort -k3,3 -k1,1n |awk '{if($2!=".")print}'>Wikim_166.sites

for i in `ls SNP.filter/Wikim_164*snps`
do 
	grep "AFDS" $i|grep -v "TK115"|cut -f1,2,9 
done>Wikim_164.tmp

sort -u Wikim_164.tmp|sort -k3,3 -k1,1n |awk '{if($2!=".")print}'>Wikim_164.sites

## get the sequences based on the SNP list
for i in `ls SNP.filter/Wikim_164*snps`
do 	
	file=`basename $i .snps`
	awk '{print $3":"$1"\t"$2}' Wikim_164.sites|sort -k1,1|join -t$'\t' -1 1 - -2 1 <(awk '{if($2!=".")print $9":"$1"\t"$3}' $i|sort -k1,1) -a 1 -o '1.1,1.2,2.1,2.2' -e 'N' |awk 'BEGIN{OFS="\t"}{if($4=="N")$4=$2;print $1,$4}'|sed 's/:/\t/g'|sort -k1,1 -k2,2n >sequence/${file}.tmp
done

for i in `ls SNP.filter/Wikim_166*snps`
do 	
	file=`basename $i .snps`
	awk '{print $3":"$1"\t"$2}' Wikim_166.sites|sort -k1,1|join -t$'\t' -1 1 - -2 1 <(awk '{if($2!=".")print $9":"$1"\t"$3}' $i|sort -k1,1) -a 1 -o '1.1,1.2,2.1,2.2' -e 'N' |awk 'BEGIN{OFS="\t"}{if($4=="N")$4=$2;print $1,$4}'|sed 's/:/\t/g'|sort -k1,1 -k2,2n>sequence/${file}.tmp
done

## create pseudo sequences for each species and combine them into a fasta file
cut -f2 Wikim_166.sites|tr '\n' '\t'|sed 's/\t//g;s/\./N/g'|sed '1i >Wikim_166'>Wikim_166.fasta
echo "" >>Wikim_166.fasta
for i in `ls sequence/Wikim_166*tmp`
do
	file=`basename $i .tmp`
	echo ">$file" >>Wikim_166.fasta
	sequence=`cut -f3 $i|tr '\n' '\t'|sed 's/\t//g;s/\./N/g'`
	echo $sequence >>Wikim_166.fasta
done

cut -f2 Wikim_164.sites|tr '\n' '\t'|sed 's/\t//g;s/\./N/g'|sed '1i >Wikim_164'>Wikim_164.fasta
echo "" >>Wikim_164.fasta
for i in `ls sequence/Wikim_164*tmp`
do
	file=`basename $i .tmp`
	echo ">$file" >>Wikim_164.fasta
	sequence=`cut -f3 $i|tr '\n' '\t'|sed 's/\t//g;s/\./N/g'`
	echo $sequence >>Wikim_164.fasta
done

## multiple sequence alignment
swarm -f mafft.swarm -t 36 --module mafft --sbatch "--mem=100g --time=72:00:00"

## approximately-maximum-likelihood phylogenetic trees generation
swarm -f fasttree.swarm -t 36 --module FastTree --sbatch "--mem=100g --time=72:00:00"
