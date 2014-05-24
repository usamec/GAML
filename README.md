Requirements
============
- Reasonable new gcc
- Boost libraries
- Velvet

Parameters forVelvet
======================

In our experiments we used only -cov_cutoff setting. Either set it to auto or
something reasonably small. We recommend you to do the same.

Compiling GAML
==============
cmake .
make

Running GAML
============

./gaml <config file>


Config file outline (check example.cfg):
========================================

Global configuration

[readset1]
readset1 configuration

[readset2]
readset2 configuration

...


Global configuration
====================

- graph=pathtograph     Required. Run Velvet on your reads and set this to
LastGraph in your Velvet output directory.
- threshold=number      Threshold for the length of long contigs. Default 500.

Read set configuration
======================
- type=read type        Required. Can be "single", "paired" or "pacbio".
- filename=filename     Required for single and pacbio reads. Filename where your reads are. Must be in fastq format.
- filename1=filename    Required for paired reads. Filename where one end of your paired
reads are. Must be in fastq format and reads are assumed to be innies.
- filename2=filename    Required for paired reads. Filename where the other end of your paired
reads are. Must be in fastq format and reads are assumed to be innies.
- insert\_mean=number   Required for paired reads. Mean insert length for paired reads.
- insert\_std=number    Required for paired reads. Standard deviation of the insert length
for paired reads.
- cache\_prefix=path    Optional. Where to put cached data from likelihood calculation.
Defaults to read set name.
- weight=number         Optional. Weight of the read set in the likelihood calculation.
Default 1.
- advice=whatever       Optional. If set then we use this read set as advice during walk extending.
- mismatch\_prob=number Optional. Probablity of errors in your reads. Defaults to 0.01.
- min\_prob\_start=number Optional. Constant c in the minimum probablity calculation.
Defaults to -10.
- min\_prob\_per\_base=number Optional. Constant k in the minimum probablity calculation.
Defaults to -0.7.
- penalty_constant=nubmer  Optional. Alpha constant in penalty for assemblies which are not 
connected enough. 
- penalty_step=number Optional. Constant k in penalty for assemblies which are not connected
enough.
- output_prefix=filename Optional. Prefix for the output files. Defaults to "output".
- starting_assembly=filename Optional. Fasta file with starting assembly.
- max_iterations=number Optional. Maximum number of iterations for simulated annealing.
Defaults to 50000.
