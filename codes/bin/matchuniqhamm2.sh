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
    progress=${progdir}/${basename}${progsuffix}
    output=${basename}${newsuffix}
    if [ ! -e ${output} ]; then
        cat ${filename} \
        | awk 'length($2)==30 {print $2 "\t" $1 "\t" 0}' \
        | ${python} ${merge_script} -w /dev/null --noheader --skipN --hamming --max_error 2 \
        2> ${progress} \
        | awk '{print $2 "\t" $1}' \
        | sort -n -r > ${output}
    fi
done
