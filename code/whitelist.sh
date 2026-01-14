# umi_tools whitelist --stdin=/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/trimmed_read3_and_1.fastq.gz \
#                     --extract-method=string \
#                     --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNN \
#                     --plot-prefix \
#                     --expect-cells=25000 \
#                     --method=reads \
#                     --knee-method=distance \
#                     --error-correct-threshold=1 \
#                     --ed-above-threshold=correct \
#                     --stdout=../input/whitelist2.txt


umi_tools whitelist --stdin=/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/trimmed_read3_and_1.fastq.gz \
                    --extract-method=string \
                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNN \
                    --plot-prefix=distance \
                    --method=umis \
                    --knee-method=distance \
                    --error-correct-threshold=1 \
                    --ed-above-threshold=discard \
                    --stdout=../input/whitelist_distance.txt
