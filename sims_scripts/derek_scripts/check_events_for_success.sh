#this will search for through each file in *.log 

for f in *.log
do
    if ! grep -q "successfully" $f; then
    echo "event did not finish successfully : $f"
    fi

    if grep -q "failed" $f; then
    echo "event failed : $f"
    fi
done
