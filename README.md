# FHM-duper

Python script to run fast homozygosity mapping.

--INPUT--<br />
Please edit line 71-74.<br />
line 71: string:working directory<br />
line 72-74: string: VCF filenames for parent1, parent2 and offspring pools.

--OUTPUT--<br />
homoScore_win[:WINSIZE:]_out.txt	 - average homozygosity score for each window<br />
homoScore_win[:WINSIZE:]_merged_out.txt  - homozygosity score for merged segments, high and low end cutoff can be edited at line 63&64.


