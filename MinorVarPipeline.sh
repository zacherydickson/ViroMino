#!/bin/bash
set -o pipefail

source "$BASHVARLIB/ExitStates.sh"
source "$BASHFUNCLIB/CheckDependency.sh"
source "$BASHFUNCLIB/CheckFile.sh"
source "$BASHFUNCLIB/IsNumeric.sh"
source "$BASHFUNCLIB/JoinBy.sh"

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
declare -a IDList
ExclusionBedFile=""
MinMapQual=30


##HANDLE INPUTS

function usage {
	>&2 echo -e "Usage: $(basename $0) Meta.txt ref.fna[.gz] refIndex\n" \
		"===This script is intended to operate on a set of samples\n" \
		"\twith a common origin. Run multiple instances otherwise\n" \
        "===INPUTS\n" \
        "\tMeta.tab\tsemicolon delim file with columns: ID, Path2Read1, Path2Read2\n" \
        "\t\tNote: will be accessed multiple times, cannot be temporary\n" \
        "\tref.fna\tFasta formatted (optionally gzipped) reference sequence\n" \
        "\trefIndex\tPrefix for the indexed reference to use\n" \
		"===OPTIONS\n" \
        "\t-q [0,∞)=$MinMapQual\tMinimum mapping quality for reads\n" \
        "\t-t [1,$MaxThreads]=$NThread\tNumber of threads to use\n" \
        "\t-w PATH=$WorkdDir\tA working directory for intermediate files;\n" \
        "\t\tWill be created if necessary\n" \
        "\t-x PATH\tA bed formatted file defining genomic regions to exclude from analysis\n" \
		"===FLAGS\n" \
        "\t-f\tForce execution of all pipeline steps\n" \
		"\t-h\tDisplay this message and exit"
        ;
}

while getopts "q:t:w:fh" opts; do
	case $opts in
        q)
            MinMapQual="$OPTARG"
            IsNumeric "$MinMapQual" NonNegReal || exit "$EXIT_FAILURE"
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
            exit $EXIT_FAILURE;
			;;
	esac
done
shift $((OPTIND-1));

if [ "$#" -lt 2 ]; then
	usage;
	exit $EXIT_FAILURE;
fi

CheckDependency lofreq      || exit "$EXIT_FAILURE";
CheckDependency bcftools    || exit "$EXIT_FAILURE";
CheckDependency samtools    || exit "$EXIT_FAILURE";
CheckDependency freebayes   || exit "$EXIT_FAILURE";
CheckDependency bowtie2     || exit "$EXIT_FAILURE";
CheckDependency fastp       || exit "$EXIT_FAILURE";
CheckDependency bgzip       || exit "$EXIT_FAILURE";

### MAIN

function main {
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
    Call "$metaFile" "$refFile" lofreq || exit "$EXIT_FAILURE" 
    Call "$metaFile" "$refFile" freebayes || exit "$EXIT_FAILURE" 
	#Apply Filters to each set of calls, exclusion regions
    #Need to split results to Allelic Primatives and RPB
	#Retain Calls from both methods
	#	re-expand sites which were retained in at least one sample
	#Report the final set of calls	
}

### FUNCTION DEFINITIONS

#Run FastP on all samples
#Inputs - MetaFile
#Output - None, creates fastpOut directory and fills it with results
#ExitCode - The Number of samples for which fastpFailed
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
function Align {
    local metaFile="$1"; shift
    local refIndex="$2"; shift
    mkdir -p "$AlignDir"
    failCount=0;
    for id in "${IDList[@]}"; do
        #Construct and globally store the aligned file name
        AlignedFileMap["$id"]="$AlignDir/$id.bam"
        #If not forcing and the alignment file exists skip this id
        [ "$Force" -eq 0 ] && [ -f "$AlignedFileMap["$id"]" ] && continue;
        unpairedStr=$(JoinBy , ${ProcessedFileMap["$id:1U"]} ${ProcessedFileMap["$id:2U"]})
        #Run Bowtie 2, then filter out unmapped and low mapQ reads and sort
        bowtie2 --threads "$NThread" --very-sensitive -x "$refIndex" \
                -1 "${ProcessedFileMap["$id:1P"]}" -2 "${ProcessedFileMap["$id:2P"]}" \
                -U "$unpairedStr" 2>| "$AlignDir/$id.log" |
            samtools view -h -F0x4 -q "$MinMapQaul" - |
            samtools sort > "${AlignedFileMap["$id"]}"
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
function Call {
    local metaFile="$1"; shift
    local refFile="$1"; shift
    local caller="$1"; shift
    mkdir -p "$CallDir"
    exitCode=0;
    callFunc="";
    case "$caller" in
        lofreq)
            callFunc="Call_lofreq"
            Call_lofreq "$metaFile" "$refFile"
            exitCode="$?"
            ;;
        freebayes)
            callFunc="Call_freebayes"
            Call_freebayes "$metaFile" "$refFile"
            exitCode="$?"
            ;;
        *)
            >&2 echo "[ERROR] Unrecognized variant caller ($caller)"
            return $(wc -l "$metaFile")
            ;;
    esac
    failCount=0;
    #TODO: This could be parallelized
    for id in "${IDList[@]}"; do
        RawVCFMap["$caller:$id"]="$CallDir/$id-$caller-raw.vcf.gz"
        #If the call file already exists, delete it
        [ "$Force" -eq 0 ] && [ -f "${RawVCFMap["$caller:$id"]}" ] && continue;
        #Make the variant calls with teh caller
        "$callFunc" "$id" "$refFile" "${AlignedFileMap["$id"]}" \
            2>| "$CallDir/$id-$caller-raw.log" |
            bgzip > "$RawVCFMap[$caller:$id]"
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
        lofreq call --ref "$refFile" --call-indels -
    #TODO: Ensure consistent output fields
}

#Run the variant caller freebayes on the samples
function Call_freebayes {
    local id=$1; shift;
    local refFile=$1; shift
    local alnFile=$1; shift
    freebayes -f "$refFile" --max-complex-gap 75 -p 1 --pooled-continuous "$alnFile"
    #TODO: Ensure consistent output fields
}



### CALL MAIN TO SIMULATE FORWARD DECLARATIONS
main "$@"
