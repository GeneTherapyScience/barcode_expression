#!/bin/bash
# tmpfile="/tmp/tmpfile.$$"
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
        tail +2 ${filename} \
        | awk 'length($1)==30 {print $2 "\t" $1}' \
        | ${python} ${merge_script} -w /dev/null --noheader --hamming --max_error 2 \
        2> ${progress} \
        | sort -n -r > ${output}
    fi
done
