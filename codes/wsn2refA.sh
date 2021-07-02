#!/bin/bash
# tmpfile="/tmp/tmpfile.$$"
maxerror="1"
suffix=".wsn"
newsuffix=".ref.uniq"
progsuffix=".ref.progress.$$"
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
        | awk 'length($1)==30 {print}' \
        | ${python} ${merge_script} -w /dev/null --skipN --drop --noheader --hamming --max_error ${maxerror} \
        2> ${progress} \
        | awk '{print $2 "\t" $1}' \
        | sort -n -r > ${output}
    fi
done
