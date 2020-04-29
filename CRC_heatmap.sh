species="Roseburia_intestinalis Roseburia_sp_CAG_303 Streptococcus_oralis Bacteroides_nordii Intestinimonas_butyriciproducens Prevotella_intermedia Fretibacterium_fastidiosum Methanobrevibacter_smithii Bacteroides_fragilis Mogibacterium_diversum Porphyromonas_asaccharolytica Eisenbergiella_tayi Bulleidia_extructa Peptostreptococcus_anaerobius Solobacterium_moorei Peptostreptococcus_stomatis Gemella_morbillorum Fusobacterium_nucleatum Dialister_pneumosintes Parvimonas_micra"

conda activate hclust2
grep -P ${species// /|}"|dataset_name|study_condition|sampleID" /shares/CIBIO-Storage/CM/news/users/francesco.beghini/chocophlan_paper/METAPHLAN3_CRC_DATA.tsv | \
sed 's/^.*s__//g; s/sampleID//g; s/dataset_name/Dataset/g; s/study_condition/Condition/g;' | \
awk 'NR < 4 {print $0 | "sort -k1r" }' | less -S 
> /shares/CIBIO-Storage/CM/news/users/francesco.beghini/chocophlan_paper/MetaPhlan3_CRC_biomarkers.tsv

hclust2.py \
    -i /shares/CIBIO-Storage/CM/news/users/francesco.beghini/chocophlan_paper/MetaPhlan3_CRC_biomarkers.tsv \
    --metadata_rows 2 \
    --no_slabels \
    --cell_aspect_ratio 90 \
    -l \
    --dpi 300 \
    --max_flabel_len 100 \
    -o MetaPhlan3_CRC_biomarkers.svg \
    --legend_file MetaPhlan3_CRC_biomarkers_legend.svg