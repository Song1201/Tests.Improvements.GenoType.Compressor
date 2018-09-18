# i is the seed, and never should be the same during different simulations. Used seed: . The files contents the names of samples are named as 'samples$i'.
for i in {1..2}
do

	# Simulate
	vtools simulate my --seed $i

	# Clear simulated dump 
	rm -r cache extracted.vcf my_resample* resample.log resample.proj 

	# gzip new simulated Population and delete the .vcf file.
	bgzip -c newSimulatedPopulation.vcf > newSimulatedPopulation.vcf.gz ; tabix -p vcf newSimulatedPopulation.vcf.gz && rm newSimulatedPopulation.vcf

	# Change the names of the samples of new simulated data.
	bcftools reheader -s samples$i -o simulatedPopulation$i.vcf.gz newSimulatedPopulation.vcf.gz 
	tabix -p vcf simulatedPopulation$i.vcf.gz && rm newSimulatedPopulation*

	#vcftools --gzvcf simulatedPopulation$i.vcf.gz

done
