## this is popmiss_filter_nomac_for_euler.sh
# also removes any RAD-tags that contain SNPs that are excessively heterozygous (Hardy-Weinberg) in any of the populations 

if [[ -z "$2" ]]; then
echo "Usage is pop_missing_filter vcffile popmap percent_missing_per_pop name_for_output"
exit 1
fi

vcftb=/cluster/project/gdc/people/schamath/tools/vcftools/bin


POPS=( `cut -f2 $2 | sort | uniq `)
scratchdirname=scratchdir_${RANDOM}

echo ${scratchdirname}
mkdir ${scratchdirname}
cd ${scratchdirname}

for i in "${POPS[@]}"
do (
grep -w $i ../$2 | cut -f1 > keep.$i

# first check missingness
$vcftb/vcftools --vcf ../$1 --keep keep.$i --missing-site --out $i 
awk '!/CHROM/' $i.lmiss | awk -v x=$3 '$6 < x' | cut -f1,2 > goodloci.${i} # has two columns: contig	pos

# now check Hardy-Weinberg heterozygote excess
$vcftb/vcftools --vcf ../$1 --keep keep.$i --hardy --out $i
cat ${i}.hwe | awk '{ if ($8 <= 0.05 ) print $1 }' | sort | uniq > excessively_heterozygous_refcontigs.${i}.txt # has one column: contig

# now remove the hardy-violators from the list of goodloci; because the awk utput is empty if the excessively_heterozygous_refcontigs.${i}.txt is empty, we have to use this if statemnt to let all goodloci pass:
awk 'NR==FNR{a[$0];next}!($1 in a)' excessively_heterozygous_refcontigs.${i}.txt goodloci.${i} > tmp.${i}

if (( $( cat tmp.${i} | wc -l ) == 0)) ; then
  echo "no excessively heterozygous contigs"
else
	mv tmp.${i} goodloci.${i} 
fi

) &
done

wait

mylist=$(ls goodloci.* | tr "\n" "\t")

#comm -12  <(sort $(echo $mylist | awk '{print $1}') ) <(sort $(echo $mylist | awk '{print $2}')) > common_loci

# iteratively finds common loci of arbitrary number of files:
mylist_array=( $mylist )
sort ${mylist_array[1]} > common_loci

for i in $mylist ;
do
echo $i
comm -12 <(sort $i ) <(sort common_loci) > tmp
mv tmp common_loci
done

wc -l common_loci

cut -f1 ../$2 > keep.overall

$vcftb/vcftools --vcf ../$1 --keep keep.overall --positions common_loci --recode --recode-INFO-all --out $4

mv $4.recode.vcf ../

# cleanup
sleep 20
cd ..
rm -r ${scratchdirname}

#####

