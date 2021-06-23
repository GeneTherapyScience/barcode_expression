#!/bin/bash
# tmpfile="/tmp/tmpfile.$$"
suffix=".uniq"
newsuffix=".hamm2.uniq"
progsuffix=".hamm2.progress.$$"
progdir="progress"
merge_script="../codes/merge_readerror.py"
python="pypy"
mkdir -p ${progdir}
for filename in "$@"
do
    basename=`basename ${filename} ${suffix}`
    output=${basename}${newsuffix}
    progress=${progdir}/${basename}${progsuffix}
    cat ${filename} \
    | awk 'length($2)==30 {print $2 "\t" $1 "\t" 0}' \
    | ${python} ${merge_script} -w /dev/null --noheader --skipN --hamming --max_error 2 \
      2> ${progress} \
    | awk '{print $2 "\t" $1}' > ${output}
done
