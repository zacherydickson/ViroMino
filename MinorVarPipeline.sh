#!/bin/bash
set -o pipefail

FullCall="$0 $*"

ExecDir="$(dirname "$(readlink -f "$0")")"
source "$ExecDir/BashFunctionLibrary/variables/ExitStates.sh"
source "$ExecDir/BashFunctionLibrary/functions/CheckDependency.sh"
source "$ExecDir/BashFunctionLibrary/functions/CheckFile.sh"
source "$ExecDir/BashFunctionLibrary/functions/IsNumeric.sh"
source "$ExecDir/BashFunctionLibrary/functions/JoinBy.sh"
source "$ExecDir/BashFunctionLibrary/functions/RandomString.sh"
RPBCmd="$ExecDir/utils/AddRPBInfoTag.sh"

##Default Values

if ! MaxThreads=$(grep -c '^proc' /proc/cpuinfo); then
    >&2 echo "$(date "+%H:%M:%S") [WARNING] Could not determine SysMax Threads - Using 1"
    MaxThreads=1;
fi

NThread=1;
WorkDir="MVPipe"
Force=0;
ProcDir="$WorkDir/processedReads"
AlignDir="$WorkDir/alignments"
CallDir="$WorkDir/calls"
CommonCallsFile="$WorkDir/CommonCalls.tab"
ExpandedVCFFile="$WorkDir/expanded.vcf.gz"
HaplotypedVCFFile="$WorkDir/haplotyped.vcf.gz"
declare -A ProcessedFileMap
declare -A AlignedFileMap
declare -A RawVCFMap
declare -A FilteredVCFMap
declare -a IDList
declare -a VCallerList=(lofreq freebayes);
ExclusionBedFile=""
MinMapQual=30
MaxHRUN=5
MinMAC=6
MaxRPB=13
MinNCallers="${#VCallerList[@]}"
MaxPileupDepth=1000
Verbose=0

##HANDLE INPUTS

function usage {
	>&2 echo -e "Usage: $(basename "$0") Meta.txt ref.fna[.gz] refIndex\n" \
		"===This script is intended to operate on a set of samples\n" \
		"\twith a common origin. Run multiple instances otherwise\n" \
        "===INPUTS\n" \
        "\tMeta.txt\tsemicolon delim file with columns: ID, Path2Read1, Path2Read2\n" \
        "\t\tNote: will be accessed multiple times, cannot be temporary\n" \
        "\tref.fna\tFasta formatted (optionally gzipped) reference sequence\n" \
        "\trefIndex\tPrefix for the indexed reference to use\n" \
        "===OUTPUT\n" \
        "\tCreates the Working directory with the following subfolder and files:\n" \
        "\t $(basename "$ProcDir")\tThe quality and adapter trimmed reads\n" \
        "\t $(basename "$AlignDir")\tThe bam files from aligning to the reference\n" \
        "\t $(basename "$CallDir")\tBoth raw and filtered calls for each sample with each caller\n" \
        "\t $(basename "$CommonCallsFile")\tThe mv sites which were agreed upon by all callers\n" \
        "\t $(basename "$ExpandedVCFFile")\tThe information at all common sites, even if the site was not\n" \
        "\t\t called in a particular sample\n" \
        "\t $(basename "$HaplotypedVCFFile")\tSame info as expanded, but with nearby variants combined\n" \
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
    while getopts "m:p:q:r:t:w:x:fvh" opts; do
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
            v)
                Verbose=1
                ;;
    		h)
    			usage;
                exit "$EXIT_FAILURE";
    			;;
            *)
                Log ERROR "Unrecognized Option ($opts)"
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
    [ "$Verbose" -eq 1 ] && Log INFO "Initialized"
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
    CommonCallsFile="$WorkDir/CommonCalls.tab"
    ExpandedVCFFile="$WorkDir/expanded.vcf.gz"
    HaplotypedVCFFile="$WorkDir/haplotyped.vcf.gz"
    #PreProcess the fastq Files - Fills in ProcessedFileMap
    PreProcess "$metaFile" || exit "$EXIT_FAILURE"
	#Align the reads - Fills in AlignedFileMap
    Align "$metaFile" "$refIndex" || exit "$EXIT_FAILURE"
	#Call MV Sites with different callers - Fills in Raw VCF Map
    for vCaller in "${VCallerList[@]}"; do
        Call "$metaFile" "$refFile" "$vCaller" || exit "$EXIT_FAILURE" 
    done
    #Filter each individual set of calls  - Fills in FilteredVCFMap
    FilterCalls "$metaFile" || exit "$EXIT_FAILURE"
    #Get the set of sites which are common to all callers - creates CommonCalls
    ReconcileCalls "$metaFile" || exit "$EXIT_FAILURE"
    #Re-expand sites which were retained in at least one sample - creates ExpandedVCF
    ReExpand "$metaFile" "$refFile" || exit "$EXIT_FAILURE"
    ##Attempt Haplotype Calling - creates HaplotypedVCF
    ##TODO:
    [ "$Verbose" -eq 1 ] && Log INFO "Done"
}

### FUNCTION DEFINITIONS

function Log {
    level="$1"; shift
    msg="$1"; shift
    >&2 echo -e "$(date "+%H:%M:%S") [$level] $msg"
}

#Run FastP on all samples
#Inputs - MetaFile
#Output - None, creates fastpOut directory and fills it with results
#ExitCode - The Number of samples for which fastpFailed
#TODO: Increase Modularity by splitting into a function that iterates and one that does the work
function PreProcess {
    [ "$Verbose" -eq 1 ] && Log INFO "Beginning Preprocessing"
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
            if [ "$bPass" -eq 1 ]; then
                [ "$Verbose" -eq 1 ] && Log INFO "\tReads for $id already preprocessed"
                continue;
            fi
        fi
        [ "$Verbose" -eq 1 ] && Log INFO "\tPreprocessing reads for $id ..."
        local code=0;
        #Run Fastp without dedup or merging
        fastp -i "$path1" -I "$path2" -o "${ProcessedFileMap["$id:1P"]}" -O "${ProcessedFileMap["$id:2P"]}" \
            --unpaired1 "${ProcessedFileMap["$id:1U"]}" --unpaired2 "${ProcessedFileMap["$id:2U"]}" \
            --detect_adapter_for_pe --overlap_diff_percent_limit 20 --overlap_diff_limit 6 \
            --length_required 30 --average_qual 30 --correction --n_base_limit 0\
            --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 \
            --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
            --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
            --report_title "$id" --thread "$NThread" --html /dev/null --json /dev/null \
            2>| "$ProcDir/${id}.log" || code="$?"
        #Check for successful fastp run, delete any partial output files generated
        if [ "$code" -ne 0 ]; then
            Log ERROR "\tFastP failure for $id"
            for segment in 1P 2P 1U 2U; do
                rm -f "${ProcessedFileMap["$id:$segment"]}"
            done
            ((failCount++))
        fi
    done < "$metaFile"
    [ "$Verbose" -eq 1 ] && Log INFO "Preprocessing complete"
    return "$failCount"
}

##Run Bowtie2 on all samples
#Inputs - metaFile
#TODO: Increase Modularity by splitting into a function that iterates and one that does the work
function Align {
    [ "$Verbose" -eq 1 ] && Log INFO "Beginning Alignment"
    local metaFile="$1"; shift
    local refIndex="$1"; shift
    local failCount=0;
    local code=0;
    mkdir -p "$AlignDir"
    for id in "${IDList[@]}"; do
        #Construct and globally store the aligned file name
        AlignedFileMap["$id"]="$AlignDir/$id.bam"
        #If not forcing and the alignment file exists skip this id
        if [ "$Force" -eq 0 ] && [ -f "${AlignedFileMap["$id"]}" ]; then
            Log INFO "\tReads already Aligned for $id"
            continue;
        fi
        [ "$Verbose" -eq 1 ] && Log INFO "\tAligning reads for $id ..."
        unpairedStr=$(JoinBy , "${ProcessedFileMap["$id:1U"]}" "${ProcessedFileMap["$id:2U"]}")
        #Run Bowtie 2, then filter out unmapped and low qual reads and sort
        bowtie2 --threads "$NThread" --very-sensitive -x "$refIndex" -1 "${ProcessedFileMap["$id:1P"]}" -2 "${ProcessedFileMap["$id:2P"]}" -U "$unpairedStr" 2>| "$AlignDir/$id.log" |
            samtools view -h -F0x4 -q "$MinMapQual" - |
            samtools sort - >| "${AlignedFileMap["$id"]}" ||
            code="$?"
        #Check for failure
        if [ "$code" -ne 0 ]; then
            Log ERROR "\tBowtie2/samtools failure for $id"
            rm -f "${AlignedFileMap["$id"]}"
            ((failCount++))
        fi
    done
    [ "$Verbose" -eq 1 ] && Log INFO "Alignment Complete"
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
    [ "$Verbose" -eq 1 ] && Log INFO "Beginning MV Calling with $caller"
    local code=0;
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
            Log ERROR "Unrecognized variant caller ($caller)"
            return "$(wc -l "$metaFile")"
            ;;
    esac
    local failCount=0;
    #TODO: This could potentially be parallelized
    for id in "${IDList[@]}"; do
        RawVCFMap["$caller:$id"]="$CallDir/$id-$caller-raw.vcf.gz"
        #If the call file already exists, skip
        if [ "$Force" -eq 0 ] && [ -f "${RawVCFMap["$caller:$id"]}" ]; then
            Log INFO "\tMVs already called with $caller for $id"
            continue;
        fi
        [ "$Verbose" -eq 1 ] && Log INFO "\tCalling $caller MVs for $id ..."
        #Make the variant calls with the caller
        "$callFunc" "$id" "$refFile" "${AlignedFileMap["$id"]}" \
            2>| "$CallDir/$id-$caller-raw.log" |
            bgzip >| "${RawVCFMap[$caller:$id]}" || code="$?"
        #Check for failure
        if [ "$code" -ne 0 ]; then
            Log ERROR "\t$caller calling failure for $id"
            rm -f "${RawVCFMap["$caller:$id"]}"
            ((failCount++))
        fi
    done
    [ "$Verbose" -eq 1 ] && Log INFO "MV Calling with $caller Complete"
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
        AddFormat2lofreq /dev/stdin
}

#Adds GT and AD format fields to an 'unknown' sample in lofreq generated output
function AddFormat2lofreq {
    inVCF="$1";shift
    awk ' #Add FORMAT/AD,GT to match 
        BEGIN {OFS="\t"}
        /^##/{print; next}
        /^#/ {
            print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Number of observations for each allele\">"
            print $0,"FORMAT","unknown"
            next
        }
        {
            match($8,/DP4=([0-9]+,){3}[0-9]+/);
            gt=0;
            split(substr($8,RSTART+4,RLENGTH-4),dp4,",");
            ad1 = dp4[1]+dp4[2]
            ad2 = dp4[3]+dp4[4]
            if(ad2 > ad1){
                gt = 1
            }
            ad=ad1","ad2
            print $0,"GT:AD",gt":"ad
        }
    ' "$inVCF"
}

#Run the variant caller freebayes on the samples
function Call_freebayes {
    local id=$1; shift;
    local refFile=$1; shift
    local alnFile=$1; shift
    local tmpFile;
    local code=0;
    tmpFile="$CallDir/fb_${id}_$(RandomString 8).tmp.vcf"
    #Call freebayes, then filter non-minor variant sites (mac < 1)
    freebayes -f "$refFile" --max-complex-gap 75 -p 1 --pooled-continuous "$alnFile" |
        FilterVCFByMac /dev/stdin 1 >| "$tmpFile" || code="$?"
    #Check if freeebays worked before continuing
    if [ "$code" -ne 0 ]; then
        Log ERROR "\tfreebayes failure for $id - Check $tmpFile"
        rm -f "$tmpFile" "$tmpFile.csi"
        return "$EXIT_FAILURE"
    fi
    #Calculate the value of HRUN for the INDELS
    AddHRUN2freebayes "$tmpFile" "$refFile" || code="$?"
    if [ "$code" -ne 0 ]; then
        Log ERROR "\tAddHRUN2freebayes failure - Check $tmpFile"
        return "$EXIT_FAILURE"
    fi
    rm -f "$tmpFile" "$tmpFile.csi"
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
#           Note the VCF file is accessed twice and so cannot be temporary
#       - a reference sequence
#Output - a vcffile with added INFO fields for HRUN
function AddHRUN2freebayes {
    local vcfFile="$1"; shift
    local refFile="$1"; shift
    local tmpFile;
    local code=0;
    tmpFile="$CallDir/fb_${id}_$(RandomString 8)_HRUN.tmp.tab"
    bcftools norm -m- --force -a "$vcfFile" 2> >(grep -v '^Lines' >&2) |
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" - |
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
        ' > "$tmpFile" || code="$?"
    if [ "$code" -ne 0 ]; then
        Log ERROR "\tFailure to Calculate Freebayes HRUN for $id"
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
    ' "$tmpFile" "$vcfFile" || code="$?"
    if [ "$code" -ne 0 ]; then 
        Log ERROR "\tFailure to add HRUN field to $vcfFile - Check $tmpFile"
        return "$EXIT_FAILURE"
    fi
    #Todo add the annotations to the vcf file
    rm -f "$tmpFile" 
}

#Runs a series of filters on a raw vcf files
# The vcf files must have AD in the format field, and HRUN Info for Indels
# Read Position Bias Information will be Added at this step
# Exclusion regions will be removed at this step
#Inputs - a metadata file 
function FilterCalls {
    metafile="$1" shift;
    [ "$Verbose" -eq 1 ] && Log INFO "Beginning MV Filtering"
    local label=""
    local id=""
    local vCaller="";
    local code=0;
    #First need to filter out any variants in the excluded regions
    #Then We add RPB values to the results
    #Then we filter on HRUN, RPB, and MAC: in that order
    #MAC must come last as the RPB filter can reveal non-minor variant sites which the MAC filter can clean up
    local failCount=0
    for label in "${!RawVCFMap[@]}"; do 
        IFS=":" read -r vCaller id <<< "$label";
        FilteredVCFMap["$label"]="$CallDir/$id-$vCaller-filt.vcf.gz"
        #If the call file already exists, skip
        if [ "$Force" -eq 0 ] && [ -f "${FilteredVCFMap["$label"]}" ]; then
            Log INFO "\tMV calls for $label already Filtered"
            continue;
        fi
        [ "$Verbose" -eq 1 ] && Log INFO "\tFiltering MV calls for $label ..."
        local filtFile="${FilteredVCFMap["$label"]}" 
        Filter "$label" "${AlignedFileMap["$id"]}" "${RawVCFMap["$label"]}" |
            bgzip >| "$filtFile" || code="$?"
        if [ "$code" -ne 0 ]; then
            Log ERROR "\tFiltering or compression failure for $label"
            ((failCount++))
            rm -f "$filtFile"
            continue;
        fi
    done
    [ "$Verbose" -eq 1 ] && Log INFO "MV Call Filtering Complete"
    return "$failCount"
}

#Perform all filtering operations for a single set of mv site calls
#Inputs -
#Output - prints filtered, uncompressed vcf file to stdout
function Filter {
    local label="$1"; shift
    local alnFile="$1"; shift
    local rawFile="$1"; shift
    local tmpFile;
    local code=0;
    tmpFile="$CallDir/filter_${label}_$(RandomString 8).tmp.vcf.gz"
    #Check if an exclusion BED file was provided, if so filter out variants at those sites
    if [ -n "$ExclusionBedFile" ]; then 
        if ! cp "$rawFile" "$tmpFile"; then
            Log ERROR "\tFailure to copy rawFile for filtering of $label"
            rm -f "$tmpFile" "$tmpFile.csi"
            return "$EXIT_FAILURE"
        fi
    else
        if ! bcftools view -T "^$ExclusionBedFile" -o "$tmpFile" -O z "$rawFile"; then
            Log ERROR "\tFailure to filter excluded regions for $label"
            rm -f "$tmpFile" "$tmpFile.csi"
            return "$EXIT_FAILURE"
        fi
    fi
    #Add RPB tags
    if ! "$RPBCmd" "$alnFile" "$tmpFile"; then
        Log ERROR "\tFailure to add RPB tag for $label - Check $tmpFile"
        return "$EXIT_FAILURE"
    fi
    #Filter by HRUN, RPB, then MAC tags, then Bgzip to final file
    bcftools filter -e "HRUN > $MaxHRUN" "$tmpFile" | 
        FilterVCFByRPB /dev/stdin "$MaxRPB" |
        FilterVCFByMAC /dev/stdin "$MinMAC" || code="$?"
    if [ "$code" -ne 0 ]; then
        Log ERROR "\nFailure in filtering by RPB/MAC/HRUN for $label - Check $tmpFile"
        return "$EXIT_FAILURE"
    fi
    rm -f "$tmpFile" "$tmpFile.csi"
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

#Take the set of calls for an individual sample and output the sites upon which all callers agree
#If any sample fails, no common calls file is generated, as incomplete data is undesirable
#Inputs - colon-delim metaFile which includes smaple ids in first column
#Output - None, writes to the CommonCallsFile
function ReconcileCalls {
    local metaFile="$1"; shift
    #If the common calls file already exists, we can skip it
    if [ "$Force" -eq 0 ] && [ -f "$CommonCallsFile" ]; then
        Log INFO "MV calls already Reconciled"
        return "$EXIT_SUCCESS";
    fi
    [ "$Verbose" -eq 1 ] && Log INFO "Beginning MV call Reconciliation"
    rm -f "$CommonCallsFile"
    #Iterate over ids
    local id;
    local path1;
    local path2;
    local failCount=0;
    echo -e "CHROM\tPOS\tREF\tALT\tCallers\tSample" >| "$CommonCallsFile"
    for id in "${IDList[@]}"; do
        [ "$Verbose" -eq 1 ] && Log INFO "\tReconciling MV calls for $id ..."
        declare -a fileArr=();
        #Get the list of filtered calls for this id
        for vCaller in "${VCallerList[@]}"; do
            local label="$vCaller:$id"
            fileArr+=("${FilteredVCFMap["$label"]}")
        done
        #Attempt reconcilliation
        if ! Reconcile "$id" "${fileArr[@]}" >> "$CommonCallsFile"; then
            Log ERROR "\tFailure to reconcile calls for $id"
            rm -f "$CommonCallsFile"
            ((failCount++))
            continue;
        fi
    done < "$metaFile"
    [ "$Verbose" -eq 1 ] && Log INFO "Completed Reconciliation"
    return "$failCount";

}

#Performs the reconciliation for a single sample
#Inputs - an id for the sample
#       - a list of filtered vcf files
#Output - a tab delim file with 5 columns:
#           chrom, pos, ref, alt, callerstr, sample id
function Reconcile {
    local id="$1";shift
    local code=0
    local tmpDir;
    tmpDir="$WorkDir/reconcile_${id}_tmp_$(RandomString 8)" || return "$EXIT_FAILURE"
    mkdir -p "$tmpDir" || return "$EXIT_FAILURE"
    failCount=0
    for filtVCF in "$@"; do
        outFile="$tmpDir/$(basename "$filtVCF")"
        bcftools norm -a -m- "$filtVCF" 2> >(grep -v '^Lines' >&2) |
            bgzip >| "$outFile" || code="$?"
        if [ "$code" -ne 0 ]; then
            Log ERROR "\tFailure to normalize $(basename "$filtVCF")";
            rm -f "$outFile"
            ((failCount++))
            continue;
        fi
        if ! bcftools index "$outFile"; then
            Log ERROR "\tFailure to index normalized $(basename "$filtVCF")";
            rm -r "$outFile.csi"
            ((failCount++))
            continue;
        fi
    done
    #Do not continue if normalization failed
    [ "$failCount" -gt 0 ] && return "$failCount"
    bcftools isec -n="$MinNCallers" "$tmpDir"/*.vcf.gz 2> >(grep -v "Note: -w" >&2) |
        awk -v id="$id" '{print $0"\t"id}' 
        #CollapseCommonMultiallelicSites "$id" /dev/stdin || return "$EXIT_FAILURE"
    rm -rf "$tmpDir"
}

##Takes the output from bcftools isec (site List) and combines multiple alt alleles into a single site call
##Inputs - an id to append to records
##       - a site list file from bcftools isec
##Output - a site list file with the id appended as the last column
#function CollapseCommonMultiallelicSites {
#    #TODO: This assumes that only one reference allele form will be present at any site,
#    #   but if both a SNP and an INDEL both exists at one site
#    # the ref forms would be different
#    #   Need to consolidate into one, all encompassing reference and alleles in the same form
#    local id="$1"; shift
#    local file="$1"; shift
#    sort -k6,6 -k1,1 -k2,2n "$file" |
#        awk -v id="$id" '
#            function output(    n,altStr,i){
#                n=asorti(altSet,altList);
#                altStr=altList[1];
#                for(i=2;i<=n;i++){altStr=altStr","altList[i]}
#                print lastChr,lastPos,commonRef,altStr,flag,id
#            }
#            BEGIN{OFS="\t"}
#            (lastChr==$1 && lastPos==$2){
#                ref=$3;
#                refLen = length(ref);
#                alt=$4
#                if(refLen < cRefLen){
#                    suffix = substr(commonRef,refLen+1); 
#                    alt=alt suffix;
#                }
#                if(refLen > cRefLen){
#                    suffix = substr(ref,cRefLen+1)
#                    commonRef=ref;
#                    cRefLen=refLen;
#                    n=asorti(altSet,altList);
#                    split("",altSet,"");
#                    for(i=1;i<=n;i++){
#                        altSet[altList[i] suffix] = 1
#                    }
#                }
#                altSet[alt]=1
#                next;
#            }
#            (lastChr){output()}
#            {
#                split("",altSet,"");
#                lastChr=$1;lastPos=$2;flag=$5
#                commonRef=$3;
#                cRefLen=length(commonRef);
#                altSet[$4]=1;
#                nAlt=1;
#            }
#            END{if(lastChr){output()}}
#        '
#}

function ReExpand {
    local metaFile="$1"; shift
    local refFile="$1"; shift
    local code=0
    #If the common calls file already exists, we can skip it
    if [ "$Force" -eq 0 ] && [ -f "$ExpandedVCFFile" ]; then
        Log INFO "MV calls already ReExpanded"
        return "$EXIT_SUCCESS";
    fi
    [ "$Verbose" -eq 1 ] && Log INFO "ReExpanding MV calls ..."
    rm -f "$ExpandedVCFFile"
    #Pileup the results and clean everything up
    bcftools mpileup --regions-file <(cut "$CommonCallsFile" -f1,2 | tail -n+2 | sort | uniq) --fasta-ref "$refFile" \
        --bam-list <( printf "%s\n" "${AlignedFileMap[@]}" ) --max-depth "$MaxPileupDepth" \
        -q "$MinMapQual" --ignore-RG --annotate "FORMAT/AD" --threads "$NThread" \
        --no-version 2> >(grep -Pv '(samples in [0-9]+ input files)|(maximum number of reads per)' >&2 ) |
        CleanAnnotations /dev/stdin | 
        AddNCallAnnote "$WorkDir" "$CommonCallsFile" /dev/stdin |
        bcftools view --trim-unseen-allele --no-version - |
        bcftools reheader --samples <(printf "%s\n" "${!AlignedFileMap[@]}") - |
        AddCCAnnote "$CommonCallsFile" /dev/stdin |
        awk -v call="$FullCall" '(FNR==1){print; print "##source="call;next}1' |
        bgzip > "$ExpandedVCFFile" || code="$?"
    if [ "$code" -ne 0 ]; then
        Log ERROR "Failure to generate expanded vcf"
        rm -f "$ExpandedVCFFile"
        return "$EXIT_FAILURE"
    fi
    [ "$Verbose" -eq 1 ] && Log INFO "MV call ReExpansion complete"
}

#Remove unwanted annotations from the VCF file
#Inputs - a vcf file to clean
#Output - Prints the modified VCF to stderr
function CleanAnnotations {
    inVcf="$1";shift
    #Get the list of default annotations to remove
    declare -a deannotateList;
    local deannotateStr
    #Get the default list of annotations
    readarray -t deannotateList < \
        <(bcftools mpileup --annotate "?" 2>&1 | grep '^\* INFO' | cut -f2 -d' ')
    #Remove the auxillary calling tags and the genotype likelihoods
    deannotateList+=(INFO/DP INFO/I16 INFO/QS FORMAT/PL)
    deannotateStr=$(JoinBy "," "${deannotateList[@]}");
    bcftools annotate --no-version -x "$deannotateStr" "$inVcf" || return "$EXIT_FAILURE"
}

#Given a file containing called MV sites and a vcf file, adds INFO/NCALLS
#   Which specifies the number of samples where the site was confidently called
#Inputs - a directory in which to place temporary files
#       - A common calls text file
#       - a VCF file, can be a temporary file
#Output - Prints the modified VCF to stdout
function AddNCallAnnote {
    local workDir=$1; shift
    local commonCalls=$1; shift
    local in="$1"; shift
    local code=0;
    local tmpFile;
    tmpFile="$workDir/reexpand-ncall_$(RandomString 8).tmp.tab.gz" 
    #Build annotation file
    awk '
        (FNR==1){next}
        {Count[$1"\t"$2]++}
        END{
            n=asorti(Count,posList);
            for(i=1;i<=n;i++){print posList[i]"\t"Count[posList[i]]
            }
        }
    ' "$commonCalls" | bgzip > "$tmpFile" || code="$?"
    #Check that construction worked and attempt indexing
    if [ "$code" -ne 0 ] || ! tabix -b2 -e2 "$tmpFile"; then
        Log ERROR "Failure to generate and index NCALL annotation"
        rm -f "$tmpFile" "$tmpFile.tbi"
        return "$EXIT_FAILURE"
    fi
    #Other annotation vars
    columnStr="CHROM,POS,INFO/NCALL"
    headerLine='##INFO=<ID=NCALL,Number=1,Type=Integer,Description="Number of samples in which this site was called">'
    #Attempt Annotation
    if ! bcftools annotate --no-version -a "$tmpFile" -c "$columnStr" -h <( echo "$headerLine") "$in"; then 
        >&2  echo "[ERROR] Failure to add INFO/NCALL Annotation - Check $tmpFile"
        return "$EXIT_FAILURE"
    fi
    #Cleanup
    rm -f "$tmpFile" "$tmpFile.tbi"
}

#Using the commonCalls File generates CC annotations to add to the input vcf
#   Annotations are added with awk, so the input must be uncompressed
#Inputs - A common calls text file
#       - a VCF file, can be a temporary file
#Output - Prints the modified VCF to stdout
function AddCCAnnote {
    local commonCalls=$1; shift
    local in="$1"; shift
    local code=0;
    awk -F '\\t' '
        BEGIN {OFS="\t"}
        (ARGIND == 1){
            if(FNR > 1){
                CCSite[$1,$2,$6]=$4;
            }
            next
        }
        /^##/{print; next}
        /^#CHROM/{
            Desc = "The Confidently Called Alt Allele(s), if any, in this sample"
            print "##FORMAT=<ID=CCA,Number=1,Type=String,Description=\""Desc"\">"
            for(i=9;i<=NF;i++){
                SampleName[i]=$i;
            }
            print;next;
        }
        {
            emptyFormat=($9==".")
            CCSite[$1,$2,"FORMAT"]="CCA"
            for(i=9;i<=NF;i++){
                CC="."
                smpl=SampleName[i]
                if(CCSite[$1,$2,smpl]) { CC=CCSite[$1,$2,smpl] }
                if(emptyFormat){ $i=CC } else { $i=$i":"CC }
            }
            print
        }
    ' "$commonCalls" "$in" || code="$?"
    if [ "$code" -ne 0 ]; then
        Log ERROR "Failure to add CC Annotation to VCF file"
        return "$EXIT_FAILURE"
    fi
}

### CALL MAIN TO SIMULATE FORWARD DECLARATIONS
if [ "${BASH_SOURCE[0]}" == "${0}" ]; then 
    main "$@"
fi
