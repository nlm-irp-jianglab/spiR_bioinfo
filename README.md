# spiR-manuscript-bioinfo
Code to reproduce Figure 6 and Figure S9 in Arp et al. 2026, "SpiR is a gut microbial enzyme that drives cholesterol conversion"

Tables/ contains all of the metagenomic data needed to replicate Figure 6

Scripts/ contains all of the scripts used to determine cholesterol converter status and generate the plot, and generate the ROC curve 

## ROC Curve
The ROC Curves in Figure 6f and Figure S9 can be replicated using Scripts/Plot-ROC.py

```
cd Scripts/

python Plot-ROC.py
```

## Citation
Arp, G., Levy, S., Jiang, A.K. et al. SpiR is a gut microbial enzyme that drives cholesterol conversion. Nat Commun 17, 3495 (2026). https://doi.org/10.1038/s41467-026-70820-6
