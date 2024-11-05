#!/usr/bin/env bash

if [ "${1}" == "-f" ]; then
    gtype="filegraph"
elif [ "${1}" == "-r" ]; then
    gtype="rulegraph"
else
    echo " usage:"
    echo " graph [-r/-f] <snakefile>"
    echo "   -r  produce a cleaned (no all rule) snakemake rule graph"
    echo "   -f  produce a cleaned snakemake file graph"
    echo "   <snakefile> is the name of the snakefile to process "
    exit 0
fi

cfile=${2}

snakemake --snakefile ${cfile} --${gtype} > full_${gtype}.dot

# Remove references to 'all' rule to avoid clutter in diagram
rulenum=$(grep "label = \"all" full_rulegraph.dot | cut -d[ -f1)
rulenum=$(echo $rulenum) # strip whitespace
grep -v "label = \"all\"" full_${gtype}.dot > temp.dot
grep -v "\-> ${rulenum}" temp.dot > ${gtype}.dot

dot -Tpng ${gtype}.dot > ${gtype}.png

