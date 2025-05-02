#!/bin/bash
set -o pipefail

ExecDir="$(dirname "$(readlink -f "$0")")"
source "$ExecDir/BashFunctionLibrary/variables/ExitStates.sh"
source "$ExecDir/BashFunctionLibrary/functions/CheckDependency.sh"
source "$ExecDir/BashFunctionLibrary/functions/CheckFile.sh"
source "$ExecDir/BashFunctionLibrary/functions/IsNumeric.sh"
source "$ExecDir/BashFunctionLibrary/functions/JoinBy.sh"
source "$ExecDir/BashFunctionLibrary/functions/RandomString.sh"
RPBCmd="$ExecDir/util/AddRPBInfoTag.sh"

##Default Values
MaxThreads=$(grep -c '^proc' /proc/cpuinfo)
if [ $? -ne 0 ]; then
    >&2 echo "[WARNING] Could not determine SysMax Threads - Using 1"
    MaxThreads=1;
fi

NThread=1;
WorkDir="MVPipe"
Force=0;
ProcDir="$WorkDir/processedReads"
AlignDir="$WorkDir/alignments"
CallDir="$WorkDir/calls"
declare -A ProcessedFileMap
declare -A AlignedFileMap
declare -A RawVCFMap
declare -A FilteredVCFMap
declare -a IDList
ExclusionBedFile=""
MinMapQual=30
MaxHRUN=5
MinMAC=6
MaxRPB=13


##HANDLE INPUTS

function usage {
	>&2 echo -e "Usage: $(basename "$0") Meta.txt ref.fna[.gz] refIndex\n" \
		"===This script is intended to operate on a set of samples\n" \
		"\twith a common origin. Run multiple instances otherwise\n" \
        "===INPUTS\n" \
        "\tMeta.tab\tsemicolon delim file with columns: ID, Path2Read1, Path2Read2\n" \
        "\t\tNote: will be accessed multiple times, cannot be temporary\n" \
        "\tref.fna\tFasta formatted (optionally gzipped) reference sequence\n" \
        "\trefIndex\tPrefix for the indexed reference to use\n" \
        "===OUTPUT\n" \
        "\tCreates the Working directory with the following subfolder and files:\n" \
        "\t processedReads\tThe quality and adapter trimmed reads\n" \
        "\t alignments\tThe bam files from aligning to the reference\n" \
        "\t calls\tBoth raw and filtered calls for each sample with each caller\n" \
        "\t commonCalls.tsv\tThe mv sites which were agreed upon by all callers\n" \
        "\t expanded.vcf.gz\tThe information at all common sites, even if the site was not\n" \
        "\t\t called in a particular sample\n" \
        "\t\thaplotyped.vcf.gz\tSame info as expanded, but with nearby variants combined\n" \
		"===OPTIONS\n" \
        "\t-m [1,∞)εZ=$MinMAC\tThe minimum number of reads supporting a minor allele\n" \
        "\t-p [0,∞)εR=$MaxRPB\tThe maximum Read Position Bias Value\n" \
        "\t-q [0,∞)εR=$MinMapQual\tMinimum mapping quality for reads\n" \
        "\t-r [0,∞)εZ=$MaxHRUN\tThe maximum allowable homopolymer length near an indel\n" \
        "\t-t [1,$MaxThreads]=$NThread\tNumber of threads to use\n" \
        "\t-w PATH=$WorkDir\tA working directory for intermediate files;\n" \
        "\t\tWill be created if necessary\n" \
        "\t-x PATH\tA bed formatted file defining genomic regions to exclude from analysis\n" \
		"===FLAGS\n" \
        "\t-f\tForce execution of all pipeline steps\n" \
		"\t-h\tDisplay this message and exit" \
        ;
}



### MAIN

function main {
    #Process Options
    while getopts "m:p:q:r:t:w:x:fh" opts; do
    	case $opts in
            m)
                MinMAC="$OPTARG"
                IsNumeric "$MinMAC" Natural || exit "$EXIT_FAILURE"
                ;;
            p)
                MaxRPB="$OPTARG"
                IsNumeric "$MaxRPB" NonNegReal || exit "$EXIT_FAILURE"
                ;;
            q)
                MinMapQual="$OPTARG"
                IsNumeric "$MinMapQual" NonNegReal || exit "$EXIT_FAILURE"
                ;;
            r)
                MaxHRUN="$OPTARG"
                IsNumeric "$MaxHRUN" Whole || exit "$EXIT_FAILURE"
                ;;
            t)
                NThread="$OPTARG";
                IsNumeric "$NThread" PosInt || exit "$EXIT_FAILURE"
                ;;
            w)
                WorkDir="$OPTARG"
                ;;
            x)
                ExclusionBedFile="$OPTARG"
                CheckFile "$ExclusionBedFile" "Exclusion Bed File" || exit "$EXIT_FAILURE"
                ;;
            f)
                Force=1
                ;;
    		h)
    			usage;
                exit "$EXIT_FAILURE";
    			;;
            *)
                >&2 echo "Unrecognized Option ($opts)"
                usage;
                exit "$EXIT_FAILURE";
                ;;
    	esac
    done
    shift $((OPTIND-1));
    #Check Args 
    if [ "$#" -lt 2 ]; then
    	usage;
    	exit "$EXIT_FAILURE";
    fi
    #Check dependencies 
    CheckDependency lofreq      || exit "$EXIT_FAILURE";
    CheckDependency bcftools    || exit "$EXIT_FAILURE";
    CheckDependency samtools    || exit "$EXIT_FAILURE";
    CheckDependency freebayes   || exit "$EXIT_FAILURE";
    CheckDependency bowtie2     || exit "$EXIT_FAILURE";
    CheckDependency fastp       || exit "$EXIT_FAILURE";
    CheckDependency bgzip       || exit "$EXIT_FAILURE";
    #Start
    local metaFile="$1"; shift
    local refFile="$1"; shift
    local refIndex="$1"; shift
    CheckFile "$metaFile" || exit "$EXIT_FAILURE"
    readarray -t IDList < <(cut -f1 -d: "$metaFile");
    #Construct the working dir and vars for various paths to use
    mkdir -p "$WorkDir"
    ProcDir="$WorkDir/processedReads"
    AlignDir="$WorkDir/alignments"
    CallDir="$WorkDir/calls"
    #PreProcess the fastq Files
    PreProcess "$metaFile" || exit "$EXIT_FAILURE"
	#Align the reads
    Align "$metaFile" "$refIndex" || exit "$EXIT_FAILURE"
	#Call MV Sites with different callers
    for vCaller in lofreq freebayes; do
        Call "$metaFile" "$refFile" $vCaller || exit "$EXIT_FAILURE" 
    done
    #Filter each individual set of calls 
    Filter "$metaFile" || exit "$EXIT_FAILURE"
    #TODO:
    #Then we normalize and atomize variants
    #then we cross check
	#Retain Calls from both methods
	#	re-expand sites which were retained in at least one sample
    #Attempt Haplotype Calling
	#Report the final set of calls	
}

### FUNCTION DEFINITIONS

#Run FastP on all samples
#Inputs - MetaFile
#Output - None, creates fastpOut directory and fills it with results
#ExitCode - The Number of samples for which fastpFailed
#TODO: Increase Modularity by splitting into a function that iterates and one that does the work
function PreProcess {
    local  metaFile="$1"; shift
    mkdir -p "$ProcDir"
    #Preprocess the raw data with fastp
	#	Do not merge reads or deduplicate
    local failCount=0;
    while IFS=":" read -r id path1 path2; do
        #Define output file names and store them in the global map
        for segment in 1P 2P 1U 2U; do
            ProcessedFileMap["$id:$segment"]="$ProcDir/${id}_$segment.fq.gz"
        done
        #If not forcing and all output files exists
        if [ "$Force" -eq 0 ]; then
            bPass=1;
            for segment in 1P 2P 1U 2U; do
                f="${ProcessedFileMap["$id:$segment"]}"
                if [ ! -f "$f" ]; then
                    bPass=0;
                    break;
                fi
            done
            [ "$bPass" -eq 1 ] && continue;
        fi
        #Run Fastp without dedup or merging
        fastp -i "$path1" -I "$path2" -o "${ProcessedFileMap["$id:1P"]}" -O "${ProcessedFileMap["$id:2P"]}" \
            --unpaired1 "${ProcessedFileMap["$id:1U"]}" --unpaired2 "${ProcessedFileMap["$id:2U"]}" \
            --detect_adapter_for_pe --overlap_diff_percent_limit 20 --overlap_diff_limit 6 \
            --length_required 30 --average_qual 30 --correction --n_base_limit 0\
            --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 \
            --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
            --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
            --report_title "$id" --thread "$NThread" --html /dev/null --json /dev/null \
            2>| "$ProcDir/${id}.log"
        #Check for successful fastp run, delete any partial output files generated
        if [ "$?" -ne 0 ]; then
            >&2 echo "[ERROR] FastP failure for $id"
            for segment in 1P 2P 1U 2U; do
                rm -f "${ProcessedFileMap["$id:$segment"]}"
            done
            ((failCount++))
        fi
    done < "$metaFile"
    return "$failCount"
}

##Run Bowtie2 on all samples
#Inputs - metaFile
#TODO: Increase Modularity by splitting into a function that iterates and one that does the work
function Align {
    local metaFile="$1"; shift
    local refIndex="$2"; shift
    mkdir -p "$AlignDir"
    local failCount=0;
    for id in "${IDList[@]}"; do
        #Construct and globally store the aligned file name
        AlignedFileMap["$id"]="$AlignDir/$id.bam"
        #If not forcing and the alignment file exists skip this id
        [ "$Force" -eq 0 ] && [ -f "${AlignedFileMap["$id"]}" ] && continue;
        unpairedStr=$(JoinBy , "${ProcessedFileMap["$id:1U"]}" "${ProcessedFileMap["$id:2U"]}")
        #Run Bowtie 2, then filter out unmapped and low mapQ reads and sort
        bowtie2 --threads "$NThread" --very-sensitive -x "$refIndex" \
                -1 "${ProcessedFileMap["$id:1P"]}" -2 "${ProcessedFileMap["$id:2P"]}" \
                -U "$unpairedStr" 2>| "$AlignDir/$id.log" |
            samtools view -h -F0x4 -q "$MinMapQual" - |
            samtools sort >| "${AlignedFileMap["$id"]}"
        #Check for failure
        if [ "$?" -ne 0 ]; then
            >&2 echo "[ERROR] Bowtie2/samtools failure for $id"
            rm -f "${AlignedFileMap["$id"]}"
            ((failCount++))
        fi
    done
    return "$failCount"
}

#Dispatcher for variant calling
#Inputs - the meta file to process
#       - the variant caller to use
#Output - None
#ExitCode, the number of ids for which calling failed
#TODO: Increase Modularity by splitting into a function that iterates and one that does the work
function Call {
    local metaFile="$1"; shift
    local refFile="$1"; shift
    local caller="$1"; shift
    mkdir -p "$CallDir"
    callFunc="";
    case "$caller" in
        lofreq)
            callFunc="Call_lofreq"
            ;;
        freebayes)
            callFunc="Call_freebayes"
            ;;
        *)
            >&2 echo "[ERROR] Unrecognized variant caller ($caller)"
            return "$(wc -l "$metaFile")"
            ;;
    esac
    local failCount=0;
    #TODO: This could potentially be parallelized
    for id in "${IDList[@]}"; do
        RawVCFMap["$caller:$id"]="$CallDir/$id-$caller-raw.vcf.gz"
        #If the call file already exists, skip
        [ "$Force" -eq 0 ] && [ -f "${RawVCFMap["$caller:$id"]}" ] && continue;
        #Make the variant calls with the caller
        "$callFunc" "$id" "$refFile" "${AlignedFileMap["$id"]}" \
            2>| "$CallDir/$id-$caller-raw.log" |
            bgzip >| "${RawVCFMap[$caller:$id]}"
        #Check for failure
        if [ "$?" -ne 0 ]; then
            >&2 echo "[ERROR] lofreq calling failure for $id"
            rm -f "${RawVCFMap["$caller:$id"]}"
            ((failCount++))
        fi
    done
    return "$failCount"
}

#Run the variant caller lofreq on the samples
function Call_lofreq {
    local id=$1; shift;
    local refFile=$1; shift
    local alnFile=$1; shift
    #Index the alignment
    if ! lofreq index "$alnFile"; then
        >&2 "[ERROR] Failure to index alignment for $id";
        return "$EXIT_FAILURE"
    fi
    lofreq indelqual --dindel --ref "$refFile" "$alnFile" |
        lofreq call --ref "$refFile" --call-indels - |
        awk ' #Add FORMAT/AD to match 
            BEGIN {OFS="\t"}
            /^##/{print; next}
            /^#/ {
                print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Number of observations for each allele\">"
                print $0,"FORMAT","unknown"
                next
            }
            {
                match($8,/DP4=([0-9]+,){3}[0-9]+/);
                split(substr($8,RSTART+4,RLENGTH-4),dp4,",");
                ad=dp4[1]+dp4[2]","dp4[3]+dp4[4]
                print $0,"AD",ad
            }

        '
}

#Run the variant caller freebayes on the samples
function Call_freebayes {
    local id=$1; shift;
    local refFile=$1; shift
    local alnFile=$1; shift
    local tmpFile;
    tmpFile="$CallDir/fb_${id}_$(RandomString 8).tmp.vcf"
    #Call freebayes, then filter non-minor variant sites (mac < 1)
    freebayes -f "$refFile" --max-complex-gap 75 -p 1 --pooled-continuous "$alnFile" |
        FilterVCFByMac /dev/stdin 1 >| "$tmpFile"
    #Check if freeebays worked before continuing
    if [ "$?" -ne 0 ]; then
        >&2 echo "[ERROR] freebayes failure for $id"
        rm -f "$tmpFile"
        return "$EXIT_FAILURE"
    fi
    #Calculate the value of HRUN for the INDELS
    AddHRUN2freebayes "$tmpFile" "$refFile" || return "$EXIT_FAILURE"
    rm -f "$tmpFile"
}

#Filters out sites where the MAC is lower than some thereshold
#   specifically it ensures there are at least 2 allelic depths greater than the threshold
#   this allows a sample with two non-reference alleles to pass
#Inputs - an uncompressed vcf file with a FORMAT AD field
#       - the MAC threshold to use
#Output - print to stdout the filtered vcf file
function FilterVCFByMAC {
    local vcf="$1"; shift
    local thresh="$1"; shift
    awk -v mac="$thresh" '
        /^#/{print;next}
        {   
            n=split($9,fmt,":");
            ADidx=-1; for(i=1;i<=n;i++){if(fmt[i]=="AD"){ADidx=i;break}}
            if(ADidx==-1){exit 1}
            split($10,val,":");
            n=split(val[ADidx],AD,",");
            nPass=0; for(i=1;i<=n;i++){if(AD[i]>=mac){nPass++}}
        }
        (nPass>1)
    '  "$vcf"
}

#Determines the HRUN value for freebayes variants
#Inputs - a freebayes generated vcf, uncompressed
#       - a reference sequence
#Output - a vcffile with added INFO fields for HRUN
function AddHRUN2freebayes {
    local vcfFile="$1"; shift
    local refFile="$1"; shift
    local tmpFile;
    tmpFile="$CallDir/fb_${id}_$(RandomString 8)_HRUN.tmp.tab"
    bcftools norm -m- --force -a "$vcfFile" |
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" |
        awk -v slop=$((MaxHRUN * 2 + 2)) ' #Convert to bed region of interest
            BEGIN {OFS="\t"}
            {rl=length($3); al=length($4);}
            (rl != al){print $1,$2,$2+slop}
        ' |
        bedtools getfasta -fi "$refFile" -bed - |
        awk ' #Get the HRUN Length
            BEGIN{OFS="\t"}
            /^>/{split(substr($1,2),a,":|-"); chrom=a[1];pos=a[2];next}
            {l=length($0); sub(substr($0,1,1)"+","",$0); print chrom,pos,l-length($0)}
        ' > "$tmpFile"
    if [ "$?" -ne 0 ]; then
        >&2 echo "[ERROR] Failure to Calculate Freebayes HRUN for $id"
        rm -f "$tmpFile"
        return "$EXIT_FAILURE"
    fi
    awk '
        BEGIN{OFS="\t"}
        (ARGIND == 1){HRUN[$1,$2]=$3;next}
        /^##/{print;next}
        /^#CHROM/{
            Desc = "Homopolymer length to the right of report indel position"
            print "##INFO=<ID=HRUN,Number=1,Type=Integer,Description=\""Desc"\">";
            print; next
        }
        (HRUN[$1,$2]){$8=$8";HRUN="HRUN[$1,$2]}1
    ' "$tmpFile" "$vcfFile" || return "$EXIT_FAILURE"
    #Todo add the annotations to the vcf file
    rm -f "$tmpFile" 
}

#Runs a series of filters on a raw vcf files
# The vcf files must have AD in the format field, and HRUN Info for Indels
# Read Position Bias Information will be Added at this step
# Exclusion regions will be removed at this step
#Inputs - a metadata file 
#TODO: Increase Modularity by splitting into a function that iterates and one that does the work
function Filter {
    metafile="$1" shift;
    local label=""
    local id=""
    local vCaller="";
    #First need to filter out any variants in the excluded regions
    #Then We add RPB values to the results
    #Then we filter on HRUN, RPB, and MAC: in that order
    #MAC must come last as the RPB filter can reveal non-minor variant sites which the MAC filter can clean up
    local failCount=0
    for label in "${!RawVCFMap[@]}"; do 
        IFS=":" read -r id vCaller <<< "$label";
        FilteredVCFMap["$label"]="$CallDir/$id-$caller-filt.vcf.gz"
        #If the call file already exists, skip
        [ "$Force" -eq 0 ] && [ -f "${FilteredVCFMap["$caller:$id"]}" ] && continue;
        local rawFile="${RawVCFMap["$label"]}" 
        local filtFile="${FilteredVCFMap["$label"]}" 
        local tmpFile;
        tmpFile="$CallDir/filter_${label}_$(RandomString 8).tmp.vcf.gz"
        #Check if an exclusion BED file was provided, if so filter out variants at those sites
        if [ -n "$ExclusionBedFile" ]; then 
            if ! cp "$rawFile" "$tmpFile"; then
                >&2 echo "[ERROR] Failure to copy rawFile for filtering of $label"
                rm -f "$tmpFile"
                ((failCount++))
                continue;
            fi
        else
            if ! bcftools view -T "^$ExclusionBedFile" -o "$tmpFile" -O z "$rawFile"; then
                >&2 echo "[ERROR] Failure to filter excluded regions for $label"
                rm -f "$tmpFile"
                ((failCount++))
                rm -f "$tmpFile"
                continue;
            fi
        fi
        #Add RPB tags
        if ! "$RPBCmd" "${AlignedFileMap["$id"]}" "$tmpFile"; then
            >&2 echo "[ERROR] Failure to add RPB tag for $label"
            ((failCount++))
            continue;
        fi
        #Filter by HRUN, RPB, then MAC tags, then Bgzip to final file
        bcftools filter -e "HRUN > $MaxHRUN" "$tmpFile" | 
            FilterVCFByRPB /dev/stdin "$MaxRPB" |
            FilterVCFByMAC /dev/stdin "$MinMAC" |
            bgzip >| "$filtFile"
        if [ "$?" -ne 0 ]; then
            >&2 echo "[ERROR] Failure to filter by RPB/MAC/HRUN for $label"
            ((failCount++))
            rm -f "$filtFile"
            continue;
        fi
        rm -f "$tmpFile"
    done
    return "$failCount"
}

#Filters out sites where the RPB is higher than some thereshold
#   specificially it determines which allelese have RPB less than threshold
#   sites with less than 2 passing allelese are filtered completely
#   any failing alleles, and their info tags are removed
#Inputs - an uncompressed vcf file with an INFO RPB field
#       - the RPB threshold to use
#Output - print to stdout the filtered vcf file
function FilterVCFByRPB {
    local vcf="$1"; shift;
    local rpb="$1"; shift;
    awk -v th="$rpb" '
        BEGIN {OFS="\t"}
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
            n=split($8,info,";")
            split("",allelePass,"")
            nPass=0;
            for(i=1;i<=n;i++){
                split(info[i],a,"=");
                key=a[1]
                val=a[2]
                if(key=="RPB"){
                    m=split(val,b,",");
                    for(j=1;j<=m;j++){
                        allelePass[j]=0;
                        if(b[j]<=th){
                            nPass++
                            allelePass[j]=1
                        }
                    }
                }
            }
            if(nPass < 2){next}
            nAlt = split($5,a,",")
            altStr=""
            for(i=1;i<=nAlt;i++){
                if(allelePass[i+1]){
                    altStr=altStr","a[i]
                }
            }
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
                formatStr=formatStr":"format[i]
            }
            $10 = substr(formatStr,2)
            print
        }
    ' "$vcf"
}

### CALL MAIN TO SIMULATE FORWARD DECLARATIONS
if [ "${BASH_SOURCE[0]}" == "${0}" ]; then 
    main "$@"
fi
