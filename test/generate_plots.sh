#! /bin/bash

JETSCAPE="/home/jetscape-user/${REPO_NAME}"
OUTPUT_DIR="${JETSCAPE}/test/pp/output/new"

#  Parse command line options
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -j)
    JETSCAPE="$2"
    shift # remove argument
    shift # remove value
    ;;
    -o)
    OUTPUT_DIR="$2"
    shift
    shift
    ;;
esac
done

echo ""
echo "Constructing observable histograms in $OUTPUT_DIR ..."
echo ""

if [ ! -d $OUTPUT_DIR ]
then
  mkdir -p $OUTPUT_DIR
fi

# Source environment
source ${JETSCAPE}/test/jetscape_analysis/init.sh
# Install some python packages -- TODO: add to docker image
pip install numba pyarrow awkward python-pptx

#--------------------------------------------------------
# Construct and plot observables 
# We adapt existing machinery in the JETSCAPE-analysis repository from the STAT group,
# since this includes calculation of many observables and comparison to experimental data.
# (In principle one could write a simpler machinery, directly computing observables
# from parquet tables -- which is likely simpler to extend to new observables)
#--------------------------------------------------------

# Construct observables
cd ${JETSCAPE}/test
python jetscape_analysis/analysis/analyze_events_validation.py -c jetscape_analysis/config/PP19_validation.yaml -i ${JETSCAPE}/build/test_out_final_state_hadrons.dat -o ${OUTPUT_DIR}

# Construct histograms
python plot/histogram_results.py -c jetscape_analysis/config/PP19_validation.yaml -i ${OUTPUT_DIR}/test_out_observables.parquet -o ${OUTPUT_DIR}

# Plot histograms and compare to previous histograms
python plot/plot_results.py -c jetscape_analysis/config/PP19_validation.yaml -i ${OUTPUT_DIR}/test_out_histograms.root -r ${OUTPUT_DIR}/../latest/test_out_histograms.root