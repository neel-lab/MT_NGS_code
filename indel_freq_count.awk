#!/usr/bin/awk

#BEGINNING OF AWK SCRIPT
{
if (!($5 in gwt) && !($5 in gid)){
	gwt[$5]=0
	gid[$5]=0
	rtype[$5]=$(NF-1)
}
if ($NF~/WT/){
	gwt[$5]=(gwt[$5]+1)
}
else if ($NF~/INDEL/){
	gid[$5]=(gid[$5]+1)
}

}

END {
for (g in gwt){
	OFS="\t"
	print g,rtype[g],gid[g],(gid[g]+gwt[g]),((gid[g])/(gid[g]+gwt[g]))
}
}
#END OF AWK SCRIPT
