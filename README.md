# sequencestatistics

A set of functions to calculate the probability to find a certain sequence by chance in a genome/proteome. Optional parameters are the equiprobability/distribution of the nucleotides/amino acid or what the sequence similarity can be.

Examples:
What is the probability of finding 'ATCTCT' by chance in a genome with uniform distribution of ACTG? 
Answer: (1/4)^6
Function: seqprobDNA('ATCTCT')
Return: 0.000244140625

What is the probability of finding 'ATCTCT' by chance in a genome with distribution {'A': 0.5, 'T': 0.25, 'C': 0.125, 'G': 0.125}? 
Answer: (1/2)(1/4)^3(1/8)^2
Function: seqprobDNA('ATCTCT', equiprob={'A': 0.5, 'T': 0.25, 'C': 0.125, 'G': 0.125})
Return: 0.0001220703125

What is the probability of finding 'ATCTCT' with a similarity of exactly 50% (f.e. ATCGGG, GTCTGA, ...) by chance in a genome with uniform distribution of ACTG?
Answer: comb(6, 3)(1/4)^3(1/4)^3(3/4)^3 (Ways of picking 3 residues out of 6 multiplied with probability of having 3 nucleotide matches and 3 non-matching residues)
Function: seqprobDNA('ATCTCT', similarity = 0.5)

**Under construction**
What is the probability of finding 'ATCTCT' with a similarity of minimum 50% (f.e. ATCGGG, GTCTGA, but also ATCTCG, AGGTCT, ...) by chance in a genome with uniform distribution of ACTG?
Answer: Idem as previous, but sum probabilities of full match, 95% match, 
Function:
