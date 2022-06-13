#! /bin/sh

./KQED

echo "diffing KQED created files and regressed files"

## regression test these old files
for exno in 2 4 5 6 7 10; do
    echo "diff-ing example${exno}.txt with REGexample${exno}.txt"
    diff example${exno}.txt REGexample${exno}.txt
done

## diff all the thread files (have to be the same!)
for th in Thread_* ; do
    echo "diff-ing $th with Thread_0"
    diff $th Thread_0
done
