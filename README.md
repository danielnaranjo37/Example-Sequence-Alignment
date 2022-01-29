# Example-Sequence-Alignment

## Protein global sequence alignment of the alpha subunit of the ammonia monooxygenase gene from AOAs, Nitrosopumilus maritimus and Nitrosomonas europeae 

### Retrieve and import FASTA protein sequence from Uniprot
library(seqinr)
library(Biostrings)

require(seqinr)

### Load in Nitrosopumilus and Nitrosomonas AMO-AFASTA files
nitrosopumilus_amoA_file = file.choose()
nitrosopumilus_amoA = read.fasta(file = nitrosopumilus_amoA_file, seqtype = "AA")[[1]]

nitrosomonas_amoA_file = file.choose()
nitrosomonas_amoA = read.fasta(file = nitrosomonas_amoA_file, seqtype = "AA")[[1]]

print(nitrosopumilus_amoA)
print(nitrosomonas_amoA)

pumilus_amoA_aaseq = readAAStringSet(nitrosopumilus_amoA_file)
monas_amoA_aaseq = readAAStringSet(nitrosomonas_amoA_file)

pumilus_string = c2s(pumilus_amoA_aaseq)
monas_string = c2s(monas_amoA_aaseq)
dotPlot(nitrosopumilus_amoA, nitrosomonas_amoA, main = "Dot plot of a protein\nwsize = 1, wstep = 1, nmatch = 1")

### ... as Needleman-Wunsch using the BLOSUM50, 62, and 80 matrices for amino acids
globalAligns1 = Biostrings::pairwiseAlignment(pumilus_string, 
                                  monas_string,
                                  substitutionMatrix = "BLOSUM80",
                                  gapOpening=3, 
                                  gapExtension=1)

globalAligns1
#### Gave a score of 175 using BLOSUM80 w/ p-val = .971
#### Gave a score of 194 using BLOSUM62 w/ p-val = .971
#### Gave a score of 331 using BLOSUM50 and later on a p-val of .001 ... this difference compared to BLOSUM62 is interesting

printPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE)
{
  require(Biostrings)           # This function requires the Biostrings package
  seq1aln <- pattern(alignment) # Get the alignment for the first sequence
  seq2aln <- subject(alignment) # Get the alignment for the second sequence
  alnlen  <- nchar(seq1aln)     # Find the number of columns in the alignment
  starts  <- seq(1, alnlen, by=chunksize)
  n       <- length(starts)
  seq1alnresidues <- 0
  seq2alnresidues <- 0
  for (i in 1:n) {
    chunkseq1aln <- substring(seq1aln, starts[i], starts[i]+chunksize-1)
    chunkseq2aln <- substring(seq2aln, starts[i], starts[i]+chunksize-1)
    # Find out how many gaps there are in chunkseq1aln:
    gaps1 <- countPattern("-",chunkseq1aln) # countPattern() is from Biostrings package
    # Find out how many gaps there are in chunkseq2aln:
    gaps2 <- countPattern("-",chunkseq2aln) # countPattern() is from Biostrings package
    # Calculate how many residues of the first sequence we have printed so far in the alignment:
    seq1alnresidues <- seq1alnresidues + chunksize - gaps1
    # Calculate how many residues of the second sequence we have printed so far in the alignment:
    seq2alnresidues <- seq2alnresidues + chunksize - gaps2
    if (returnlist == 'FALSE')
    {
      print(paste(chunkseq1aln,seq1alnresidues))
      print(paste(chunkseq2aln,seq2alnresidues))
      print(paste(' '))
    }
  }
  if (returnlist == 'TRUE')
  {
    vector1 <- s2c(substring(seq1aln, 1, nchar(seq1aln)))
    vector2 <- s2c(substring(seq2aln, 1, nchar(seq2aln)))
    mylist <- list(vector1, vector2)
    return(mylist)
  }
}


#### Print global alignment 

printPairwiseAlignment(globalAligns1, chunksize = 60, FALSE)


##### Find best local alignment between our amoA subunit genes 
localAlignAmoA = pairwiseAlignment(pumilus_string, monas_string,  gapOpening = 0, gapExtension = 10, scoreOnly = FALSE, type="local")

localAlignAmoA #### Print out optimal local alignment and its score of 15.2793

printPairwiseAlignment(localAlignAmoA, chunksize =  20)

print(localAlignAmoA@pattern) ### Common pattern is IVGATY using a gapOpening of -2 and gapExtension of -8

#### Stays IVGATY with gapOpening = 0 and gapExtension of -10 through +10

### Test statistical significance of our alignment using a multimodal model for protein sequences 
#### -- Probabilities for different amino acids are set equal to their frequencies in an input sequnece

generateSeqsWithMultinomialModel <- function(inputsequence, X)
{
  #### Change the input sequence into a vector of letters
  require("seqinr") # This function requires the SeqinR package.
  inputsequencevector <- s2c(inputsequence)
  ####Find the frequencies of the letters in the input sequence "inputsequencevector":
  mylength <- length(inputsequencevector)
  mytable <- table(inputsequencevector)
  #### Find the names of the letters in the sequence
  letters <- rownames(mytable)
  numletters <- length(letters)
  probabilities <- numeric() # Make a vector to store the probabilities of letters
  for (i in 1:numletters)
  {
    letter <- letters[i]
    count <- mytable[[i]]
    probabilities[i] <- count/mylength
  }
  #### Make X random sequences using the multinomial model with probabilities "probabilities"
  seqs <- numeric(X)
  for (j in 1:X)
  {
    seq <- sample(letters, mylength, rep=TRUE, prob=probabilities) # Sample with replacement
    seq <- c2s(seq)
    seqs[j] <- seq
  }
  #### Return the vector of random sequences
  return(seqs)
}

### Generate random seqs to compare ours to 
pumilus_random_seqs = generateSeqsWithMultinomialModel(pumilus_string, 1000)
monas_random_seqs = generateSeqsWithMultinomialModel(monas_string, 1000)

### Test our protein sequences against these randomly generated sequences
#### -- Nitrosopumilus vs Nitrosomonas' random seqs

pairwiseAlignment(pumilus_string, monas_random_seqs[1], 
                  substitutionMatrix = "BLOSUM80",
                  gapOpening=3, 
                  gapExtension=1, 
                  scoreOnly = TRUE)

randomscores = double(1000)
for (i in 1:1000)
{
  score <- pairwiseAlignment(pumilus_string, monas_random_seqs[i], 
                             substitutionMatrix = "BLOSUM80",
                             gapOpening = 3, gapExtension = 1, scoreOnly = TRUE)
  randomscores[i] <- score
}

hist(randomscores, col = "blue")

sum(randomscores >= globalAligns1@score) ### BLOSUM50: 1/1000, BLOSUM62: 1/1000, BLOSUM80: 1/1000   

##### This gives us a value of 1, which means exactly 1 of the random seqs generated gives an alignment score greater than or equal to the real alignment score between Nitrosopumilus and Nitrosomonas amoA genes
#### -- Repeat with Nitrosomonas vs Nitrosopumilus' random seqs

pairwiseAlignment(monas_string, pumilus_random_seqs[1], 
                  substitutionMatrix = "BLOSUM80",
                  gapOpening=3, 
                  gapExtension=1, 
                  scoreOnly = TRUE)

randomscores = double(1000)
for (i in 1:1000)
{
  score <- pairwiseAlignment(monas_string, pumilus_random_seqs[i], 
                             substitutionMatrix = "BLOSUM80",
                             gapOpening = 3, gapExtension = 1, scoreOnly = TRUE)
  randomscores[i] <- score
}
print(randomscores[1000])

hist(randomscores, col = "red")

sum(randomscores >= globalAligns1@score) ### BLOSUM50: 1/1000, BLOSUM62: 1/1000, BLOSUM80: 0/1000

#### This gives us a value of 1, which means exactly 1 of the random seqs generated gives an alignment score greater than or equal to the real alignment score between Nitrosopumilus and Nitrosomonas amoA genes
#### ... i.e. that our p-value is 1/1000 = 0.001, and this maintains between both amino acid sequences

#### To do next: play with parameters to get the most accurate readings. Current parameters are fairly arbitrary and not based on starting methionine or anything reasonable like that. 

### Final notes: a lot of this code was adapted from an online tutorial that can be found at: 
#### ... https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter4.html#:~:text=For%20example%2C%20if%20amino%20acid,%3D50%20and%20y%20%3D53.

