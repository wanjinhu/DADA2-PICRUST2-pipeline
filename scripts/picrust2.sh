#! /bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate picrust2
picrust2_pipeline.py -s $1 -i $2 -o $3 -p 16
# EC
add_descriptions.py -i $3/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o $3/EC_metagenome_out/Picrust2_EC_descrip.tsv
# KO
add_descriptions.py -i $3/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o $3/KO_metagenome_out/Picrust2_KO_descrip.tsv
# METACYC
add_descriptions.py -i $3/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o $3/pathways_out/Picrust2_METACYC_descrip.tsv
# KEGG_PATHWAY
pathway_pipeline.py -i $3/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
                    -o $3/KEGG_pathways --no_regroup \
                    --map /root/biosoft/picrust2-2.5.0/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv
gunzip -k $3/KEGG_pathways/path_abun_unstrat.tsv.gz
mv $3/KEGG_pathways/path_abun_unstrat.tsv $3/KEGG_pathways/Picrust2_KEGG_PATHWAY.tsv
conda deactivate

## demo
# source activate
# conda activate picrust2
# picrust2_pipeline.py -s ASV.fasta -i ASV_table.xls -o Picrust2 -p 16
# cd Picrust2
# EC
# add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
#                     -o EC_metagenome_out/Picrust2_EC_descrip.tsv
# KO
# add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
#                     -o KO_metagenome_out/Picrust2_KO_descrip.tsv
# METACYC
# add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
#                     -o pathways_out/Picrust2_METACYC_descrip.tsv
# KEGG_PATHWAY
# pathway_pipeline.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
#                     -o KEGG_pathways --no_regroup \
#                     --map /root/biosoft/picrust2-2.5.0/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv
# conda deactivate
