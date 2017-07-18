# Output Specification for Simulation Data Workflow

### Output folder structure
In the decompressed data folder, there are folders with name pattern SVRNA${i}_${j}, where ${i} represents the number of SVs simulated, and ${j} represents the replicate. Alignment data and TSV calling results are stored in each SVRNA${i}_${j} folder. The structure of each folder follows:
- SVRNA${i}_${j}
	- WholeGenome          // folder for genome sequence and STAR and BWA indexes
	- reads                // simulated RNA-seq reads
	- Alignments           // folder for alignments
		- Star_rearranged
			- Aligned.sortedByCoord.out.bam
			- Chimeric.out.bam
			- Merged.bam
			- Merged.bam.bai
		- SpeedSeq
			- Aligned.bam
			- Aligned.bam.bai
			- Aligned.discordants.bam
			- Aligned.discordants.bam.bai
			- Aligned.splitters.bam
			- Aligned.splitters.bam.bai
		- TSVcall
			- DELLY_BWA         // result of DELLY2 with SpeedSeq main aligned BAM
				- delly_DEL.bcf // deletions
				- delly_DUP.bcf // duplications
				- delly_INV.bcf // inversions
				- delly_TRA.bcf // translocations
				- delly_SV.vcf  // merged TSVs of DEL, DUP, INV, TRA
				- hit.txt       // accuracy file
			- DELLY_STAR        // result of DELLY2 with merged STAR concordant and chimeric alignment BAM
				- delly_DEL.bcf // deletions
				- delly_DUP.bcf // duplications
				- delly_INV.bcf // inversions
				- delly_TRA.bcf // translocations
				- delly_SV.vcf  // merged TSVs of DEL, DUP, INV, TRA
				- hit.txt       // accuracy file
			- LUMPY                     // LUMPY result with SpeedSeq alignment files 
				- lumpy_thresh${k}.vcf  // predictions under read support threshold to be ${k}
				- hit_thresh${k}.txt    // accuracy of the predictions with read support threshold ${k}
			- SQUID_BWA                 // SQUID result with SpeedSeq main alignment BAM
				- squid_w${k}_sv.txt       // predictions with read support threshold ${k}
				- squid_w${k}_sv_hit.txt   // accuracy of predictions with read support threshold ${k}
			- SQUID_STAR                // SQUID result with STAR concordant and chimeric BAM
				- squid_w${k}_sv.txt       // predictions with read support threshold ${k}
				- squid_w${k}_sv_hit.txt   // accuracy of predictions with read support threshold ${k}
			- TransABySS             // TransABySS result with Gmap and MUMMER3
				- mergedassembly.fa  // assembled sequence by TransABySS
				- gmap.out           // alignment of transcript sequences to reference genome (genome_rearranged.fa) by Gmap
				- gmap_res.bedpe     // processed TSVs from Gmap alignment
				- gmap_res_hit.txt   // accuracy file of Gmap TSVs
				- nucm.delta         // alignment of transcript sequences to reference genome (genome_rearranged.fa) by MUMMER3
				- nucm_text.txt      // alignment of transcript sequences to reference genome by MUMMER3 in human readable format
				- nucm_res.bedpe     // processed TSVs from MUMMER3 alignment
				- nucm_res_hit.txt   // accuracy file of MUMMER3 TSVs
				- other internal files and folder of TransABySS
			- Trinity                // Trinity result with Gmap and MUMMER3
				- Trinity.fasta      // assembled sequence by Trinity
				- gmap.out           // alignment of transcript sequences to reference genome (genome_rearranged.fa) by Gmap
				- gmap_res.bedpe     // processed TSVs from Gmap alignment
				- gmap_res_hit.txt   // accuracy file of Gmap TSVs
				- nucm.delta         // alignment of transcript sequences to reference genome (genome_rearranged.fa) by MUMMER3
				- nucm_text.txt      // alignment of transcript sequences to reference genome by MUMMER3 in human readable format
				- nucm_res.bedpe     // processed TSVs from MUMMER3 alignment
				- nucm_res_hit.txt   // accuracy file of MUMMER3 TSVs
				- other internal files and folder of TransABySS

### Accuracy file specification
- First line: overall accuracy and sensitivity of predictions. Note that when calculating accuracy of DELLY2 and LUMPY, we exclude from the denominator the deletion type of predictions, but still keep the record in file.
- Records contain the following columns:
	1. chromosome of the first TSV breakpoint
	2. start position of junction sequence of the first breakpoint
	3. end position of junction sequence of the first breakpoint
	4. strand of junction sequence of the first breakpoint
	5. chromosome of the second TSV breakpoint
	6. start position of junction sequence of the second breakpoint
	7. end position of junction sequence of the second breakpoint
	8. strand of junction sequence of the second breakpoint
	9. whether prediction is correct (0 for wrong, 1 for correct)
	10. the ID of corresponding correct TSVs or -1