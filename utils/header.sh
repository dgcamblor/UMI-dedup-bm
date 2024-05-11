BLUE="\033[0;34m"
NC="\033[0m"

function header {
    step=$1
    echo -e "${BLUE}------------------------------------------------------------------------${NC}"
    echo -e "${BLUE}Benchmarking step:${NC} $step"
    echo -e "${BLUE}Initiating analysis with profile:${NC} $PROFILE"  #$ Two different profiles: BRP and UMIVAR, defined in config.conf
    echo -e "${BLUE}Target:${NC} $TARGET"
    echo -e "${BLUE}Reference genome:${NC} $REF_GENOME"
    echo -e "${BLUE}Read structure:${NC} $READ_STRUCTURE"

    if [ $PREPROCESSING -eq 1 ]; then
        PREPROCESSING_WORD="Yes"
    else
        PREPROCESSING_WORD="No"
    fi

    echo -e "${BLUE}Preprocessing:${NC} $PREPROCESSING_WORD"
    echo -e "${BLUE}------------------------------------------------------------------------${NC}"
}