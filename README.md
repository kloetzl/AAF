AAF (Alignment and Assembly Free)
===

This is a package for constructing phylogeny without doing alignment or assembly.

###Bootstrap
The most feedback I received about AAF are around bootsrap. It is very computationally intensive to do the two-step nonparametric bootstrap. In case you have a higher coverage (>8X), we assume that the incomplete coverage problem is minor. To reduce the computational load, you can choose only to carry out the seconde step of the bootstrap (nonparametric\_bootstrap\_s2only.py): sample the kmer table with replacement 1/k of the number of the rows of the table. To further reduce the computaiton, here is a version to sample from the shared kmer table (nonparametric\_bootstrap\_s2only_skt.py). Singletons from each sample (i.e. kmers that only appear in one sample) are calculated from the difference between the total diversity file and the shared kmer table. Then those singletons are added back during the the calculation of pariwise distance, following a poisson distribution with a mean of 1/k of each singletone number.  
####BetaVersion/nonparametric\_bootstrap\_s2only.py

####BetaVersion/nonparametric\_bootstrap\_s2only_skt.py

This only does ONE boostrap. It is designed this way since some users use high throughput facilities. For high performance facility users, increase the ram and threads so each boostrap takes less time. You can wrap this script with a shell script. Be sure not to overwrite the boostrap tree generated each time.

Example:

	python singletonCalculator.py phylokmer.dat.gz kmer_diversity.wc -t 10
	[This would produce a file containing the number of singletons in each sample, in this case phylokmer_singleton.wc]
	for i in {1:100}: #boostrap 100 times
	do
		python nonparametric_bootstrap_s2only_skt.py phylokmer.dat.gz phylokmer_singleton.wc -t 10
		cat phylokmer_bootstrap.tre >> phylokmer_bootstrap
	done
	consense #use phylokmer_bootstrap_trees as infile
		

### FAQ
1. Dear User: If I have paired end (sample.1.fq, sample.2.fq) files for each sample, should I merge
them as input for AAF or should I keep them separately in the ./data/ folder?

	Huan: If you have multiple files for one sample, please put them in the same folder. AAF detects things in one folder as one sample and take the name of the folder as the sample name. Unfortunately AAF does not deal with a mixture of folders and files in the data directory. Therefore if you have one sample that has multiple input files, the rest need to be in folders as well, even if some of them only have one sequence file. Of course you could merge input files from one sample into one so there are only files in the data directory. This way no subdirectories need to be made. Either way it should work. Just no mix of files and folders. I hope I’m not making this sounds more complicated than it needs to be. 

2. Dear User: Should I use the BetaVersion?

	Huan: Like any BetaVersion, it might not work on your machine and most importantly, it might not be consistant with the user manual. But let's be reckless and give it a try! Please email me or report an issue if it does not work. Thanks for your help!