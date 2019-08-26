#!/usr/bin/awk

#BEGINNING OF SCRIPT
{
if ($4~/76M/){
	OFS="\t"
	# The read is a perfect match, label as WT
	print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14,"WT"
}

else {
#CIGAR indicates indel occurring:
# Split the CIGAR strings by spans of matches
split($4,a,/[[:digit:]].?M/,seps)
OFS="\t"
# Begin at the start position of the read:
strpos=$2

# EXAMPLE
# If the string is 50M5D26M:
# The characters stored in "a" are :
# BLANK; 5D; BLANK
# The characters stored in "seps" are :
# 50M; 26M
# Start looking at the split results at second index

for(i=2;i<=length(a);++i){
	if (length(a[i])==0 && i==length(a)){
		# No record in a[i], and at the end of the splits
		# If the code didn't break by now, there is no match:
		print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14,"WT"
		break
	}	
	strpos=strpos+substr(seps[i-1],1,(length(seps[i-1])-1))
	msize=substr(a[i],1,(length(a)-1))
	
	#Assumes the indel size is greater than the 3bp window
	if (($9-strpos < msize) && ($9-strpos >=0)){
		# Indel position contains CRISPR target
		# Break the for loop, print the result, move to next line
		print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14,"INDEL"
		break
	}
	#Assumes the indel size is smaller than the 3bp window
	else if (($10-(strpos+msize)>=0) && (strpos-$9>=0)) {
		# Indel is within the CRISPR target region
		# Break the for loop, print the result, move to next line
		print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14,"INDEL"
                break
	}
	#If there's no hit, add the indel size to string position to keep looking for mutations
	strpos=strpos+msize
	
	#Check to see if looping is done, and print no match
	if (i==length(a)){
		print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14,"WT"
	}	
	
}

}
#END OF ELSE STATEMENT
}
#END OF SCRIPT
