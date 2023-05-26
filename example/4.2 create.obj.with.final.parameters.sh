#!/bin/bash

sample_list='sample.list.selected.txt'
cd output_dir
while read line; do
    mkdir ${line}
done < ${sample_list}

Rscript="R4.0/bin/Rscript --vanilla"
create_obj="4. create_obj.R"
output_dir=output_dir

"""
example:
${Rscript} ${create_obj} ST100 ScaleData 30 0.8 ${output_dir}
"""
${Rscript} ${create_obj} sample1 ScaleData(or ST) npcs res ${output_dir}
