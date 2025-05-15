#!/usr/bin/env -S awk -f

(ARGC < 2){
    print "Provide a vcf file and a rpb threshold" > "/dev/stderr"
}

BEGIN {
    th=ARGV[2]
    ARGC--
    OFS="\t"
}

/^##INFO/ {
    split(substr($0,12),a,",");
    split(a[2],b,"=");
    InfoNumType[a[1]]=b[2];
    
}

/^##FORMAT/ {
    split(substr($0,14),a,",");
    split(a[2],b,"=");
    FormatNumType[a[1]]=b[2];
}

/^#/{print; next}

{
    #Determine which alleles pass the threshold
    n=split($8,info,";")
    split("",allelePass,"")
    #tracks the number of alleles to the left of this allele which failed
    split("",allelePassOffset,"") 
    nPass=0;
    nFail=0;
    for(i=1;i<=n;i++){
        split(info[i],a,"=");
        key=a[1]
        val=a[2]
        if(key=="RPB"){
            m=split(val,b,",");
            for(j=1;j<=m;j++){
                allelePass[j]=0;
                allelePassOffset[j]=nFail;
                if(b[j]<=th){
                    nPass++
                    allelePass[j]=1
                } else {
                    nFail++
                }
            }
        }
    }
    #Done if there aren't two passing alleles
    if(nPass < 2){next}
    nAlt = split($5,a,",")
    altStr=""
    for(i=1;i<=nAlt;i++){
        if(allelePass[i+1]){
            altStr=altStr","a[i]
        }
    }
    #Process any multi-allelic INFO tags
    $5=substr(altStr,2)
    infoStr=""
    for(i=1;i<=n;i++){
        split(info[i],a,"=");
        key=a[1]
        val=a[2]
        numType=InfoNumType[key]
        if(numType == "A" || numType == "R" || numType == "G"){
            m=split(val,b,",")
            off = 0
            s=","b[1]
            if(numType == "A") {s="";off = 1}
            for(j=2-off;j<=m;j++){
                if(allelePass[j+off]){
                    s = s "," b[j]
                }
            }
            info[i]=key"="substr(s,2)
        }
        infoStr=infoStr";"info[i]
    }
    #Process any multi-allelic FORMAT tags
    $8=substr(infoStr,2)
    n=split($9,fmt,":")
    split($10,format,":")
    formatStr=""
    for(i=1;i<=n;i++){
        key=fmt[i]
        val=format[i]
        numType=FormatNumType[key]
        if(numType == "A" || numType == "R" || numType == "G"){
            m=split(val,b,",")
            off = 0
            s=","b[1]
            if(numType == "A") {s="";off = 1}
            for(j=2-off;j<=m;j++){
                if(allelePass[j+off]){
                    s = s "," b[j]
                }
            }
            format[i]=substr(s,2)
        }
        #Special processing for Genotype fields
        if(key == "GT"){
            ploidy=split(val,gtArr,"\\||\\/");
            sep="|"
            if(val ~ "/"){sep = "/"}
            gtStr=""
            for(j=1;j<=ploidy;j++){
                gt = gtArr[j]
                if(gt != "."){
                    gt -= allelePassOffset[gt+1]
                }
                gtStr=gtStr sep gt
            }
            format[i]=substr(gtStr,2)
        }
        formatStr=formatStr":"format[i]
    }
    $10 = substr(formatStr,2)
    print
}
