Running this application
====================

For each input bed file you need to process it first using Misha Bilinky's tool for computing region coverage:

```
java -jar RegionsCoverageFromBEDCalculator.jar \
-b input_bedgraph_file -o ./ \
-s data/mm9.chrom.sizes \
-r data/Ref-seq-exons.bed \
-n NAME (eg. CASTEiJ_C57BL6J-CASTEiJ)
```

Next you can run the python script:

```
python post_alea.py -g data/refFlat.txt -1 data/120307-TE-dels_CAST.bed -2 data/120307-TE_CAST_sorted.bed \
--m1 maternalStrain1.coverage \
--p2 paternalStrain2.coverage \
--p1 paternalStrain1.coverage \
--m2 maternalStrain2.coverage
```
