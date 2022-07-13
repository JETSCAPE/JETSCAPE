#! /bin/bash

JETSCAPE="/home/jetscape-user/JETSCAPE"
ANALYSIS_CONFIG="${JETSCAPE}/test/pp/config/jetscapeTestConfig.yaml"
OUTPUT_DIR="${JETSCAPE}/test/pp/output/new"
REFERENCE_DIR="${JETSCAPE}/test/pp/output/latest"

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
    -c)
    ANALYSIS_CONFIG="$2"
    shift
    shift
    ;;
    -o)
    OUTPUT_DIR="$2"
    shift
    shift
    ;;
    -r)
    REFERENCE_DIR="$2"
    shift
    shift
    ;;
esac
done

echo ""
echo "Running Tests..."
echo ""

cd ${JETSCAPE}/test/jetscape_analysis/generate
python jetscape_events.py -c ${ANALYSIS_CONFIG} -o ${OUTPUT_DIR} -j ${JETSCAPE}

echo ""
echo "Comparing tests in $OUTPUT_DIR to Reference in $REFERENCE_DIR ..."
echo ""

if [ ! -d $OUTPUT_DIR ]
then
  mkdir -p $OUTPUT_DIR
fi

N=0
N_PASSED_HEPMC=0
N_PASSED_ASCII=0
cd $PREFIX/$OUTPUT_DIR
for dir in */ ; do
  N=$((N+1))
  
  DIFF_HEPMC=$(diff $OUTPUT_DIR/$dir/test_out.hepmc $REFERENCE_DIR/${dir}test_out.hepmc)
  if [ $? -ne 0 ]
  then
    echo "Error: Check whether you have used the same YAML config for the Test and the Reference"
    echo "New file: $OUTPUT_DIR/$dir/test_out.hepmc"
    echo "Reference file: $REFERENCE_DIR/${dir}test_out.hepmc"
    exit 1
  fi
  
  DIFF_ASCII=$(diff $OUTPUT_DIR/$dir/test_out.dat $REFERENCE_DIR/${dir}test_out.dat)
  if [ $? -ne 0 ]
  then
    echo "Error: Check whether you have used the same YAML config for the Test and the Reference"
    echo "New file: $OUTPUT_DIR/$dir/test_out.hepmc"
    echo "Reference file: $REFERENCE_DIR/${dir}test_out.hepmc"
    exit 1
  fi

  if [ "${DIFF_HEPMC}" == "" ]
  then
    N_PASSED_HEPMC=$((${N_PASSED_HEPMC}+1))
  else
    echo "Test $dir failed for HepMC"
  fi

  if [ "${DIFF_ASCII}" == "" ]
  then
    N_PASSED_ASCII=$((${N_PASSED_ASCII}+1))
  else
    echo "Test $dir failed for Ascii"
  fi
    
done

N_FAILED_HEPMC=$(($N-$N_PASSED_HEPMC))
N_FAILED_ASCII=$(($N-$N_PASSED_ASCII))
if [[ $N_FAILED_HEPMC -eq 0 && $N_FAILED_ASCII -eq 0 ]]
then
  echo "All $N tests passed! :)"
else
  echo ""
  echo "Tests FAILED :("
  echo "$N_FAILED_HEPMC/$N tests FAILED for HepMC"
  echo "$N_FAILED_ASCII/$N tests FAILED for Ascii"
  exit 1
fi
echo ""
