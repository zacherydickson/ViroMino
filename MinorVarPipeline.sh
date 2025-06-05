#!/bin/bash
set -o pipefail


FullCall="$0 $*"

ExecDir="$(dirname "$(readlink -f "$0")")"
source "$ExecDir/BashFunctionLibrary/variables/ExitStates.sh"
source "$ExecDir/BashFunctionLibrary/functions/CheckDependency.sh"
source "$ExecDir/BashFunctionLibrary/functions/CheckFile.sh"
source "$ExecDir/BashFunctionLibrary/functions/CheckVersion.sh"
source "$ExecDir/BashFunctionLibrary/functions/IsNumeric.sh"
source "$ExecDir/BashFunctionLibrary/functions/JoinBy.sh"
source "$ExecDir/BashFunctionLibrary/functions/RandomString.sh"
RPBCmd="$ExecDir/utils/AddRPBInfoTag.sh"
PileupCmd="$ExecDir/utils/CustomPileup.pl"
FilterVCFByRPBCmd="$ExecDir/utils/FilterVCFByRPB.awk"
VERSION="$(cat "$ExecDir/Version.txt")"

##Default Values

if ! MaxThreads=$(grep -c '^proc' /proc/cpuinfo); then
    >&2 echo "$(date "+%H:%M:%S") [WARNING] Could not determine SysMax Threads - Using 1"
    MaxThreads=1;
fi

NThread=1;
WorkDir="MVPipe"
Force=0;
ForceFrom="";
ProcDir="$WorkDir/processedReads"
AlignDir="$WorkDir/alignments"
CallDir="$WorkDir/calls"
ConfigFile="$WorkDir/Config.txt"
ErrEstFile="$WorkDir/ErrorRateEstimate.tab"
CommonCallsFile="$WorkDir/CommonCalls.tab"
ExpandedVCFFile="$WorkDir/expanded.vcf.gz"
HaplotypedVCFFile="$WorkDir/haplotyped.vcf.gz"
declare -A ProcessedFileMap
declare -A AlignedFileMap
declare -A RawVCFMap
declare -A FilteredVCFMap
declare -a IDList
declare -a VCallerList=(lofreq freebayes);
declare -a PipelineStepList=(init process align errest call filter reconcile expand haplotype)
declare -A PipelineStepIdxMap;
NPipelineStep=1;
for step in "${PipelineStepList[@]}"; do
    PipelineStepIdxMap["$step"]="$NPipelineStep";
    ((NPipelineStep++))
done
ExclusionBedFile=""
MinMapQual=30
MaxHRUN=5
MACAlpha=0.01
MinMAF=0.01
MaxRPB=13
ReadLen=150
MinNCallers="${#VCallerList[@]}"
MaxPileupDepth=1000
Verbose=0

##HANDLE INPUTS

function usage {
	>&2 echo -e "Usage: $(basename "$0") Meta.txt ref.fna[.gz] refIndex\n" \
		"===This script is intended to operate on a set of samples\n" \
		"\twith a common origin. Run multiple instances otherwise\n" \
        "===INPUTS\n" \
        "\tMeta.txt\tcolon delim file with columns: ID, Path2Read1, Path2Read2\n" \
        "\t\tNote: will be accessed multiple times, cannot be temporary\n" \
        "\tref.fna\tFasta formatted (optionally gzipped) reference sequence\n" \
        "\trefIndex\tPrefix for the bowtie2 indexed reference to use\n" \
        "===OUTPUT\n" \
        "\tCreates the Working directory with the following subfolders and files:\n" \
        "\t $(basename "$ProcDir")\tThe quality and adapter trimmed reads\n" \
        "\t $(basename "$AlignDir")\tThe bam files from aligning to the reference\n" \
        "\t $(basename "$CallDir")\tBoth raw and filtered calls for each sample with each caller\n" \
        "\t $(basename "$ConfigFile")\tRecord of settings and version used to generate the results\n" \
        "\t $(basename "$ErrEstFile")\tEstimated Per-Base Error Rates for each sample\n" \
        "\t $(basename "$CommonCallsFile")\tThe mv sites which were agreed upon by all callers\n" \
        "\t $(basename "$ExpandedVCFFile")\tThe information at all common sites, even if the site was not\n" \
        "\t\t called in a particular sample\n" \
        "\t $(basename "$HaplotypedVCFFile")\tSame info as expanded, but with nearby variants combined\n" \
		"===OPTIONS\n" \
        "\t-f STR\tForce execution of all steps including and after this one\n" \
        "\t\tSTR must be one of ${PipelineStepList[*]}\n" \
        "\t-l [1,∞)εZ=$ReadLen\tThe length of input reads\n" \
        "\t-m (0,1)εZ=$MACAlpha\tThe confidence level for Minimum Minor Allele Count\n" \
        "\t-M B|[0,1]εR=$MinMAF\tThe minimum minor allele frequency filter mode; One of\n" \
        "\t\tB - Greater than Per Base Error Rate; PBER\n" \
        "\t\tA value between 0 and 1, exclusive - use a fixed minimum MAF\n" \
        "\t-p [0,∞)εR=$MaxRPB\tThe maximum Read Position Bias Value\n" \
        "\t-q [0,∞)εR=$MinMapQual\tMinimum mapping quality for reads\n" \
        "\t-r [0,∞)εZ=$MaxHRUN\tThe maximum allowable homopolymer length near an indel\n" \
        "\t-t [1,$MaxThreads]=$NThread\tNumber of threads to use\n" \
        "\t-w PATH=$WorkDir\tA working directory for intermediate files;\n" \
        "\t\tWill be created if necessary\n" \
        "\t-x PATH\tA bed formatted file defining genomic regions to exclude from analysis\n" \
		"===FLAGS\n" \
        "\t-F\tForce execution of all pipeline steps; alias for -f init\n" \
        "\t-v\tVerbose Logging information\n" \
        "\t-V\tPrint the Version and exit\n" \
		"\t-h\tDisplay this message and exit" \
        ;
}



### MAIN

function main {
    #Process Options
    while getopts "f:l:m:M:p:q:r:t:w:x:FvVh" opts; do
    	case $opts in
            f)
                ForceFrom="$OPTARG";
                if [ -z "${PipelineStepIdxMap["$ForceFrom"]}" ]; then
                    Log ERROR "Attempt to force after uncrecognized pipeline phase ($ForceFrom)"
                    exit "$EXIT_FAILURE"
                fi
                ;;
            l)
                ReadLen="$OPTARG"
                IsNumeric "$ReadLen" PosInt || exit "$EXIT_FAILURE"
                ;;
            m)
                MACAlpha="$OPTARG"
                IsNumeric "$MACAlpha" OpenUnitIV || exit "$EXIT_FAILURE"
                ;;
            M)
                MinMAF="$OPTARG"
                if [ "$MinMAF" != "B" ]; then
                    IsNumeric "$MinMAF" UnitIV || exit "$EXIT_FAILURE"
                fi
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
                ProcDir="$WorkDir/processedReads"
                AlignDir="$WorkDir/alignments"
                CallDir="$WorkDir/calls"
                ConfigFile="$WorkDir/Config.txt"
                ErrEstFile="$WorkDir/ErrorRateEstimate.tab"
                CommonCallsFile="$WorkDir/CommonCalls.tab"
                ExpandedVCFFile="$WorkDir/expanded.vcf.gz"
                HaplotypedVCFFile="$WorkDir/haplotyped.vcf.gz"
                ;;
            x)
                ExclusionBedFile="$OPTARG"
                CheckFile "$ExclusionBedFile" "Exclusion Bed File" || exit "$EXIT_FAILURE"
                ;;
            F)
                ForceFrom=init
                Force=1
                ;;
            v)
                Verbose=1
                ;;
            V)
                >&2 echo "V$VERSION"
                exit "$EXIT_FAILURE"
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
    CheckAllDependencies || exit "$EXIT_FAILURE";
    #Start
    [ "$Verbose" -eq 1 ] && Log INFO "Initialized"
    [ "$ForceFrom" == "init" ] && Force=1;
    local metaFile="$1"; shift
    local refFile="$1"; shift
    local refIndex="$1"; shift
    CheckFile "$metaFile" || exit "$EXIT_FAILURE"
    mkdir -p "$WorkDir"
    CheckConfig "$metaFile" "$refFile" "$refIndex" || exit "$EXIT_FAILURE"
    #Index the input reference
    samtools faidx "$refFile"
    readarray -t IDList < <(cut -f1 -d: "$metaFile");
    #Construct the working dir
    #PreProcess the fastq Files - Fills in ProcessedFileMap
    [ "$ForceFrom" == "process" ] && Force=1;
    PreProcess "$metaFile" || exit "$EXIT_FAILURE"
	#Align the reads - Fills in AlignedFileMap
    [ "$ForceFrom" == "align" ] && Force=1;
    Align "$metaFile" "$refIndex" || exit "$EXIT_FAILURE"
	#Call MV Sites with different callers - Fills in Raw VCF Map
    [ "$ForceFrom" == "errest" ] && Force=1;
    EstimateError "$metaFile" || exit "$EXIT_FAILURE"
    [ "$ForceFrom" == "call" ] && Force=1;
    for vCaller in "${VCallerList[@]}"; do
        Call "$metaFile" "$refFile" "$vCaller" || exit "$EXIT_FAILURE" 
    done
    #Filter each individual set of calls  - Fills in FilteredVCFMap
    [ "$ForceFrom" == "filter" ] && Force=1;
    FilterCalls "$metaFile" || exit "$EXIT_FAILURE"
    #Get the set of sites which are common to all callers - creates CommonCalls
    [ "$ForceFrom" == "reconcile" ] && Force=1;
    ReconcileCalls "$metaFile" || exit "$EXIT_FAILURE"
    #Re-expand sites which were retained in at least one sample - creates ExpandedVCF
    [ "$ForceFrom" == "expand" ] && Force=1;
    ReExpand "$metaFile" "$refFile" || exit "$EXIT_FAILURE"
    ##Attempt Haplotype Calling - creates HaplotypedVCF
    [ "$ForceFrom" == "haplotype" ] && Force=1;
    ##TODO:
    [ "$Verbose" -eq 1 ] && Log INFO "Done"
}

### FUNCTION DEFINITIONS

function Log {
    level="$1"; shift
    msg="$1"; shift
    >&2 echo -e "$(date "+%H:%M:%S") [$level] $msg"
}

function CheckAllDependencies {
    local code=0;
    declare -a dependList=( lofreq bcftools samtools freebayes \
                            bowtie2 fastp bgzip bedtools Rscript \
                            perl)
    #Check if each dependency exists
    for depend in "${dependList[@]}"; do
        CheckDependency "$depend" || return "$EXIT_FAILURE";
    done
    #Get the current versions of all dependencies
    declare -A versions;
    #lofreq
    versions["cur:lofreq"]="0";
    versions["min:lofreq"]="0";
    versions["max:lofreq"]="";
    #bcftools
    versions["cur:bcftools"]="$(bcftools --version |
                                awk '(NR==1){print $2}')";
    versions["min:bcftools"]="1.21";
    versions["max:bcftools"]="";
    #samtools
    versions["cur:samtools"]="$(samtools --version | 
                                awk '(NR==1){print $2}')";
    versions["min:samtools"]="1.21";
    versions["max:samtools"]="";
    #freebayes
    versions["cur:freebayes"]="$(   freebayes --version | 
                                    awk '(FNR == 1){print $2}')"
    versions["min:freebayes"]="1.3.6";
    versions["max:freebayes"]="";
    #bowtie2
    versions["cur:bowtie2"]="$( bowtie2 --version | 
                                awk '(FNR == 1){print $3}')";
    versions["min:bowtie2"]="2.5.4";
    versions["max:bowtie2"]="";
    #fastp
    versions["cur:fastp"]="$(   fastp --version 2>&1 |
                                awk '(FNR == 1){print $2}')";
    versions["min:fastp"]="0.23.4";
    versions["max:fastp"]="";
    #bgzip
    versions["cur:bgzip"]="$(   bgzip --version |
                                awk '(FNR == 1){print $3}')";
    versions["min:bgzip"]="1.21";
    versions["max:bgzip"]="";
    #bedtools
    versions["cur:bedtools"]="$(bedtools --version |
                                awk '(FNR == 1){print $2}')";
    versions["min:bedtools"]="2.31.1";
    versions["max:bedtools"]="";
    #Rscript
    versions["cur:Rscript"]="$( Rscript --version |
                                awk '(FNR == 1){print $4}')";
    versions["min:Rscript"]="4.3.3";
    versions["max:Rscript"]="";
    #Perl
    versions["cur:perl"]="$(perl --version |
                            grep -P 'v[0-9]+' |
                            awk -F '\\(|)' '{print $2}')"
    versions["min:perl"]="5.32.1";
    versions["max:perl"]="";
    #Check Versions
    for depend in "${dependList[@]}"; do
        code=0;
        CheckVersion    "${versions["cur:$depend"]}" \
                        "${versions["min:$depend"]}" \
                        "${versions["max:$depend"]}"; code="$?";
        if [ "$code" -ne 0 ]; then
            Log ERROR "Installed $depend (${versions["cur:$depend"]}) does not meet requirements ${versions["min:$depend"]} ${versions["max:$depend"]}"
            return "$EXIT_FAILURE";
        fi
    done
    #Special Perl Module Checks
    declare -a moduleList=(Bio::Seq Bio::DB::HTS);
    versions["min:Bio::Seq"]="1.7.8"
    versions["max:Bio::Seq"]=""
    versions["min:Bio::DB::HTS"]="3.01"
    versions["max:Bio::DB::HTS"]=""
    for module in "${moduleList[@]}"; do
        script="print \$$module::VERSION"
        if ! version=$(perl -M"$module" -e "$script" 2> /dev/null); then
            Log ERROR "Could not locate Perl Module $module";
            return "$EXIT_FAILURE"
        fi
        CheckVersion    "$version" \
                        "${versions["min:$module"]}"\
                        "${versions["max:$module"]}"; code="$?"
        if [ "$code" -ne 0 ]; then
            Log ERROR "Installed perl module $module ($version) does not meet requirements ${versions["min:$module"]} ${versions["max:$module"]}";
            return "$EXIT_FAILURE"
        fi
    done
}

function CheckConfig {
    local metaFile;
    local refFile;
    local refIndex;
    local excBedFile;
    metaFile=$(readlink -f "$1"); shift
    refFile=$(readlink -f "$1"); shift
    refIndex=$(readlink -f "$1"); shift
    excBedFile=$(readlink -f "$ExclusionBedFile");
    local code=0;
    #If the Config File Doesn't exist it can be created
    if ! [ -s "$ConfigFile" ]; then
        WriteConfig "$metaFile" "$refFile" "$refIndex" "$excBedFile"; code="$?";
        return $code;
    fi
    bUpdateConfig=0;
    #Iterate over the config file and check for differences
    while IFS="=" read -r key value; do
        target="";
        forceFromMin=init;
        case "$key" in
            Version)
                target="$VERSION";
                forceFromMin=init;
                ;;
            metaFile)
                target="$metaFile"
                forceFromMin=init;
                ;;
            refFile)
                target="$refFile"
                forceFromMin=align;
                ;;
            refIndex)
                target="$refIndex"
                forceFromMin=align;
                ;;
            readLen)
                target="$ReadLen"
                forceFromMin=esterr;
                ;;
            macAlpha)
                target="$MACAlpha"
                forceFromMin=filter;
                ;;
            minMAF)
                target="$MinMAF"
                forceFromMin=filter;
                ;;
            maxRPB)
                target="$MaxRPB"
                forceFromMin=filter;
                ;;
            minMapQ)
                target="$MinMapQual"
                forceFromMin=align;
                ;;
            maxHRUN)
                target="$MaxHRUN"
                forceFromMin=call;
                ;;
            excludeBED)
                target="$excBedFile"
                forceFromMin=filter;
                ;;
        esac
        forceFromIdx="$NPipelineStep"
        [ -n "$ForceFrom" ] && forceFromIdx="${PipelineStepIdxMap[$ForceFrom]}"
        #Check if the call value matches the config value
        # if it doesn't match check if force after is set to the matching step
        if [ "$value" != "$target" ]; then
            bUpdateConfig=1;
            if [ "${PipelineStepIdxMap[$forceFromMin]}" -lt "$forceFromIdx" ]; then
                Log ERROR "Previous Run's $key ($value) does not match current value ($target)";
                Log INFO "Use '-f $forceFromMin' to overwrite or select a different Working directory (-w) to proceed"
                return "$EXIT_FAILURE"
            fi
        fi
    done < "$ConfigFile"
    if [ "$bUpdateConfig" -eq 1 ]; then
        WriteConfig "$metaFile" "$refFile" "$refIndex" "$excBedFile"; code="$?";
    fi
}

function WriteConfig {
    local metaFile="$1"; shift
    local refFile="$1"; shift
    local refIndex="$1"; shift
    local excBedFile="$1"; shift
    printf "%s\n" "Version=$VERSION" "metaFile=$metaFile" "refFile=$refFile" "refIndex=$refIndex" \
                    "readLen=$ReadLen" "macAlpha=$MACAlpha" "minMAF=$MinMAF" "maxRPB=$MaxRPB" \
                    "minMapQ=$MinMapQual" "maxHRUN=$MaxHRUN" "excludeBED=$excBedFile" \
        >| "$ConfigFile"
}

#Run FastP on all samples
#Inputs - MetaFile
#Output - None, creates fastpOut directory and fills it with results
#ExitCode - The Number of samples for which fastpFailed
#TODO: Increase Modularity by splitting into a function that iterates and one that does the work
function PreProcess {
    [ "$Verbose" -eq 1 ] && Log INFO "Beginning Preprocessing - PreProcDir=$ProcDir"
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
        local fastpLog="$ProcDir/${id}.log"
        #Run Fastp without dedup or merging
        fastp -i "$path1" -I "$path2" -o "${ProcessedFileMap["$id:1P"]}" -O "${ProcessedFileMap["$id:2P"]}" \
            --unpaired1 "${ProcessedFileMap["$id:1U"]}" --unpaired2 "${ProcessedFileMap["$id:2U"]}" \
            --detect_adapter_for_pe --overlap_diff_percent_limit 20 --overlap_diff_limit 6 \
            --length_required 30 --average_qual 30 --correction --n_base_limit 0\
            --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 \
            --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
            --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
            --report_title "$id" --thread "$NThread" --html /dev/null --json /dev/null \
            2>| "$fastpLog" ; code="$?"
        #Check for successful fastp run, delete any partial output files generated
        if [ "$code" -ne 0 ]; then
            Log ERROR "\tFastP failure for $id - Check $fastpLog"
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
    [ "$Verbose" -eq 1 ] && Log INFO "Beginning Alignment - AlnDir=$AlignDir"
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
            [ "$Verbose" -eq 1 ] && Log INFO "\tReads already Aligned for $id"
            continue;
        fi
        [ "$Verbose" -eq 1 ] && Log INFO "\tAligning reads for $id ..."
        unpairedStr=$(JoinBy , "${ProcessedFileMap["$id:1U"]}" "${ProcessedFileMap["$id:2U"]}")
        local bowtie2Log="$AlignDir/$id.log"
        #Run Bowtie 2, then filter out unmapped and low qual reads and sort
        bowtie2 --threads "$NThread" --very-sensitive -x "$refIndex" -1 "${ProcessedFileMap["$id:1P"]}" -2 "${ProcessedFileMap["$id:2P"]}" -U "$unpairedStr" 2>| "$bowtie2Log" |
            samtools view -h -F0x704 - |
            samtools sort -n - |
            samtools fixmate -m - - |
            samtools view -hq "$MinMapQual" - |
            samtools sort - |
            samtools markdup -r - -  >| "${AlignedFileMap["$id"]}";  code="$?"
        #Check for failure
        if [ "$code" -ne 0 ]; then
            Log ERROR "\tBowtie2/samtools failure for $id - Check $bowtie2Log"
            rm -f "${AlignedFileMap["$id"]}"
            ((failCount++))
        fi
    done
    [ "$Verbose" -eq 1 ] && Log INFO "Alignment Complete"
    return "$failCount"
}

function EstimateError {
    [ "$Verbose" -eq 1 ] && Log INFO "Beginning Error Rate Estimation - OutFile=$ErrEstFile"
    local metaFile="$1"; shift
    #If the error rate file already exists, we can skip it
    if [ "$Force" -eq 0 ] && [ -f "$ErrEstFile" ]; then
        Log INFO "Error Rates already estimated"
        return "$EXIT_SUCCESS";
    fi
    local code=0;
    #Pileup the viral aligned reads and pull out the allelelic depths
    #Count the number of reads at each site which do not support the consensus allele
    #Add pseudocounts (1 non-consensus and 1 consensus read) to deal with division by zero
    bcftools mpileup -d "$MaxPileupDepth" --no-reference "$AlignDir"/*.bam 2> >(grep -Pv '(samples in [0-9]+ input files)|(maximum number of reads per)' >&2 ) |
        bcftools query -HH -f "[%AD\t]\n" | 
        awk -v rl="$ReadLen" 'BEGIN{OFS="\t"; print "Sample\tPerBase\tPerRead"}
            /^#/{for(i=1;i<=NF;i++){n=split($i,a,"/"); Label[i]=substr(a[n],1,length(a[n])-7)};next}
            {for(i=1;i<=NF;i++){
                n=split($i,a,",");
                t=0;m=-1;
                for(j=1;j<=n;j++){t+=a[j];if(a[j]>m){m=a[j]}}
                T[i]+=t-m;N[i]+=t}
            }
            END{for(k in T){
                pber=(T[k]+1)/(N[k]+2);
                prer=1-(1-pber)^rl;
                print Label[k],pber,prer
            }}' >| "$ErrEstFile"; code="$?"
    if [ "$code" -ne 0 ]; then
        Log ERROR "Error Rate Estimation Failure"
        rm -f "$ErrEstFile"
        return "$EXIT_FAILURE"
    fi

}

#Dispatcher for variant calling, runs the variant caller functions and adds RPB values
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
        local label="$caller:$id"
        RawVCFMap[$label]="$CallDir/$id-$caller-raw.vcf.gz"
        #If the call file already exists, skip
        if [ "$Force" -eq 0 ] && [ -f "${RawVCFMap["$label"]}" ]; then
            [ "$Verbose" -eq 1 ] && Log INFO "\tMVs already called with $caller for $id"
            continue;
        fi
        [ "$Verbose" -eq 1 ] && Log INFO "\tCalling $caller MVs for $id ..."
        local callLog="$CallDir/$id-$caller-raw.log";
        #Make the variant calls with the caller
        "$callFunc" "$id" "$refFile" "${AlignedFileMap["$id"]}" \
            2>| "$callLog" |
            bgzip >| "${RawVCFMap[$label]}"; code="$?"
        #Check for failure
        if [ "$code" -ne 0 ]; then
            Log ERROR "\t$caller calling failure for $id - Check $callLog"
            rm -f "${RawVCFMap["$label"]}"
            ((failCount++))
        fi
        #Add RPB tags
        if ! "$RPBCmd" "${AlignedFileMap["$id"]}" "${RawVCFMap[$label]}"; then
            Log ERROR "\tFailure to add RPB tag for $label"
            rm -f "${RawVCFMap["$label"]}"
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
        AddContigs2lofreq "$refFile" /dev/stdin |
        AddFormat2lofreq /dev/stdin
}

#Uses the faidx index for the reference genome to add contig entries to the lofreq results
function AddContigs2lofreq {
    refFile="$1"; shift
    inVCF="$1"; shift
    awk -F '\\t' '
        (ARGIND == 1){ contigStr[++nContig] = "##contig=<ID="$1",length="$2">";next}   
        (FNR == 1){ print; for(i=1;i<=nContig;i++){print contigStr[i]};next}
        1
    ' "$refFile.fai" "$inVCF"
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
        FilterVCFByMAC /dev/stdin count 1 >| "$tmpFile"; code="$?"
    #Check if freeebays worked before continuing
    if [ "$code" -ne 0 ]; then
        Log ERROR "\tfreebayes failure for $id - Check $tmpFile"
        rm -f "$tmpFile" "$tmpFile.csi"
        return "$EXIT_FAILURE"
    fi
    #Calculate the value of HRUN for the INDELS
    AddHRUN2freebayes "$tmpFile" "$refFile"; code="$?"
    if [ "$code" -ne 0 ]; then
        Log ERROR "\tAddHRUN2freebayes failure - Check $tmpFile"
        return "$EXIT_FAILURE"
    fi
    rm -f "$tmpFile" "$tmpFile.csi"
}

#Filters out sites where the MAC is lower than some thereshold
#   specifically it ensures there are at least 2 allelic depths greater than the threshold
#   this allows a sample with two non-reference alleles to pass
#NOTE: It is required that the input VCF have a DP info field
#Inputs - an uncompressed vcf file with a FORMAT AD field
#       - the threshold mode, count if the threshold is to be interpreted as a read count otherwise interpreted
#           as an alpha value
#       - a threshold to use
#Output - print to stdout the filtered vcf file
function FilterVCFByMAC {
    local vcf="$1"; shift
    local mode="$1"; shift
    local thresh="$1"; shift
    local tmpFile;
    tmpFile="$CallDir/macfilt_$(RandomString 8).tmp.vcf"
    if ! cat "$vcf" >| "$tmpFile"; then
        rm -f "$tmpFile"
        return "$EXIT_FAILURE"
    fi
    local code=0;
    awk '
        BEGIN{MACLine=1}
        (ARGIND == 1){MAC[FNR]=$1;next}
        /^#/{print;next}
        {   
            n=split($9,fmt,":");
            ADidx=-1; for(i=1;i<=n;i++){if(fmt[i]=="AD"){ADidx=i;break}}
            if(ADidx==-1){exit 1}
            split($10,val,":");
            n=split(val[ADidx],AD,",");
            nPass=0; for(i=1;i<=n;i++){if(AD[i]>=MAC[MACLine]){nPass++}}
            MACLine++
        }
        (nPass>1)
    '  <(CalcMAC "$tmpFile" "$mode" "$thresh") "$tmpFile"; code="$?"
    if [ "$code" -ne 0 ]; then
        Log ERROR "\tFailure to Filter by MAC ($mode $th) - Check $tmpFile"
        return "$EXIT_FAILURE"
    fi
    rm -f "$tmpFile"
}

function CalcMAC {
    local vcf="$1"; shift
    local mode="$1"; shift
    local th="$1"; shift
    local code=0;
    #Mode where a specific value is used at all sites
    if [ "$mode" == "count" ]; then
        bcftools view -H "$vcf" |
            awk -v mac="$th" '{print mac}'; code="$?"
        return "$code"
    #Mode where a minor allele count is calculated from a confidence level
    elif [ "$mode" == "alpha" ]; then
        { echo "$th"; bcftools query -f "%DP\n" "$vcf"; } |
            Rscript -e 'n <- file("stdin") |> readLines() |> as.numeric(); zsq <- qnorm(n[1])^2; (n[-1]*zsq/(n[-1]+zsq)) |> ceiling() |> as.character() |> writeLines()'; code="$?"
        return "$code"
    #Mode where a minimum allele frequency is used
    elif [ "$mode" == "af" ]; then
        { echo "$th"; bcftools query -f "%DP\n" "$vcf"; } |
            Rscript -e 'n <- file("stdin") |> readLines() |> as.numeric(); af <- n[1]; n <- n[-1]; (n*af) |> ceiling() |> as.character() |> writeLines()'; code="$?"
        return "$code";
    fi
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
        ' > "$tmpFile"; code="$?"
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
    ' "$tmpFile" "$vcfFile"; code="$?"
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
    #Then we filter on HRUN, RPB, and MAC: in that order
    #MAC must come last as the RPB filter can reveal non-minor variant sites which the MAC filter can clean up
    local failCount=0
    for label in "${!RawVCFMap[@]}"; do 
        IFS=":" read -r vCaller id <<< "$label";
        af="$MinMAF"
        #Allelle frequency threshold based on Per-Base Error Rate
        if [ "$af" == "B" ]; then
            af=$(awk -v s="$id" '($1 == s){print $2}' "$ErrEstFile")
        fi
        FilteredVCFMap["$label"]="$CallDir/$id-$vCaller-filt.vcf.gz"
        #If the call file already exists, skip
        if [ "$Force" -eq 0 ] && [ -f "${FilteredVCFMap["$label"]}" ]; then
            [ "$Verbose" -eq 1 ] && Log INFO "\tMV calls for $label already Filtered"
            continue;
        fi
        [ "$Verbose" -eq 1 ] && Log INFO "\tFiltering MV calls for $label ..."
        local filtFile="${FilteredVCFMap["$label"]}" 
        Filter "$label" "${AlignedFileMap["$id"]}" "${RawVCFMap["$label"]}" "$af" |
            bgzip >| "$filtFile" ; code="$?"
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
    local af="$1"; shift
    local tmpFile;
    local code=0;
    tmpFile="$CallDir/filter_${label}_$(RandomString 8).tmp.vcf.gz"
    #Check if an exclusion BED file was provided, if so filter out variants at those sites
    if [ -z "$ExclusionBedFile" ]; then 
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
    #Filter by HRUN, RPB, then MAC tags, then Bgzip to final file
    bcftools filter -e "HRUN > $MaxHRUN" "$tmpFile" | 
        "$FilterVCFByRPBCmd" /dev/stdin "$MaxRPB" |
        FilterVCFByMAC /dev/stdin "af" "$af" | 
        FilterVCFByMAC /dev/stdin "alpha" "$MACAlpha"; code="$?"
    if [ "$code" -ne 0 ]; then
        Log ERROR "\nFailure in filtering by RPB/MAC/HRUN for $label - Check $tmpFile"
        return "$EXIT_FAILURE"
    fi
    rm -f "$tmpFile" "$tmpFile.csi"
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
            ((failCount++))
            continue;
        fi
    done < "$metaFile"
    if [ "$failCount" -gt 0 ]; then
        Log ERROR "Could not generate Common Calls"
        rm -f "$CommonCallsFile"
        return "$failCount"
    fi
    [ "$Verbose" -eq 1 ] && Log INFO "Completed Reconciliation"
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
            bgzip >| "$outFile"; code="$?"
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
    if [ "$(CountSites "$tmpDir"/*.vcf.gz)" -gt 0 ]; then
        bcftools isec -n="$MinNCallers" "$tmpDir"/*.vcf.gz 2> >(grep -v "Note: -w" >&2) |
            awk -v id="$id" '{print $0"\t"id}' 
            #CollapseCommonMultiallelicSites "$id" /dev/stdin || return "$EXIT_FAILURE"
    fi
    rm -rf "$tmpDir"
}

function CountSites {
    for f in "$@"; do 
        bcftools view -H "$f";
    done | wc -l
}

#Takes the output from bcftools isec (site List) and combines multiple alt alleles into a single site call
#Inputs - an id to append to records
#       - a site list file from bcftools isec
#Output - a site list file with the id appended as the last column
function CollapseCommonMultiallelicSites {
    local id="$1"; shift
    local file="$1"; shift
    sort -k6,6 -k1,1 -k2,2n "$file" |
        awk -v id="$id" '
            function output(    n,altStr,i){
                n=asorti(altSet,altList);
                altStr=altList[1];
                for(i=2;i<=n;i++){altStr=altStr","altList[i]}
                print lastChr,lastPos,commonRef,altStr,flag,id
            }
            BEGIN{OFS="\t"}
            (lastChr==$1 && lastPos==$2){
                ref=$3;
                refLen = length(ref);
                alt=$4
                if(refLen < cRefLen){
                    suffix = substr(commonRef,refLen+1); 
                    alt=alt suffix;
                }
                if(refLen > cRefLen){
                    suffix = substr(ref,cRefLen+1)
                    commonRef=ref;
                    cRefLen=refLen;
                    n=asorti(altSet,altList);
                    split("",altSet,"");
                    for(i=1;i<=n;i++){
                        altSet[altList[i] suffix] = 1
                    }
                }
                altSet[alt]=1
                next;
            }
            (lastChr){output()}
            {
                split("",altSet,"");
                lastChr=$1;lastPos=$2;flag=$5
                commonRef=$3;
                cRefLen=length(commonRef);
                altSet[$4]=1;
                nAlt=1;
            }
            END{if(lastChr){output()}}
        '
}

function ReExpand {
    local metaFile="$1"; shift
    local refFile="$1"; shift
    local code=0
    #If the common calls file already exists, we can skip it
    if [ "$Force" -eq 0 ] && [ -f "$ExpandedVCFFile" ]; then
        [ "$Verbose" -eq 1 ] && Log INFO "MV calls already ReExpanded"
        return "$EXIT_SUCCESS";
    fi
    [ "$Verbose" -eq 1 ] && Log INFO "ReExpanding MV calls ..."

    rm -f "$ExpandedVCFFile"
    #Pileup the results and clean everything up
    "$PileupCmd" "$refFile" "$CommonCallsFile" <(ConstructBamMapFile) |
        awk -v call="$FullCall" '(FNR==1){print; print "##source="call;next}1' |
        bgzip > "$ExpandedVCFFile"; code="$?"
    if [ "$code" -ne 0 ]; then
        Log ERROR "Failure to generate expanded vcf"
        rm -f "$ExpandedVCFFile"
        return "$EXIT_FAILURE"
    fi
    [ "$Verbose" -eq 1 ] && Log INFO "MV call ReExpansion complete"
}

function ConstructBamMapFile {
    for id in "${!AlignedFileMap[@]}"; do
        echo -e "$id\t${AlignedFileMap[$id]}";
    done
}

### CALL MAIN TO SIMULATE FORWARD DECLARATIONS
if [ "${BASH_SOURCE[0]}" == "${0}" ]; then 
    main "$@"
fi
