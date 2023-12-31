% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_seed.R
\name{get_seed}
\alias{get_seed}
\title{Get the target seed sequence given a canonical seed name and input sequence}
\usage{
get_seed(guide.seq, seed.name = "mer7m8", start.pos = 1, stop.pos = 8)
}
\arguments{
\item{guide.seq}{A character string greater than 8 bp and oriented 5'-> 3'.}

\item{seed.name}{The seed name of interest. Options: mer8,  mer7A1, mer7m8,
mer6. If not in the default list, the start.pos and stop.pos arguments
will be used to define the seed.}

\item{start.pos}{The start position for a custom seed definition}

\item{stop.pos}{The stop position for a custom seed definition}
}
\value{
An object with the entries:
\itemize{
\item \code{Guide}: Input guide sequence. Input is expected to be RNA.
\item \code{Seed.Name}: The seed name.
\item \code{Seed.Seq.RNA}: The seed sequence as a \code{RNAString}
\item \code{Seed.Seq.DNA}: The seed sequence as a \code{DNAString}
\item \code{Target.Seq}: The target DNA sequence based on the reverse complement of
the seed as a \code{DNAString}
}
}
\description{
Given a sequence greater than 8 bp oriented 5' -> 3' and a seed
definition, this function will return an object containing seed-specific sequence
information. Users can input a custom seed name, but must provide the start
position (\code{start.pos}) and stop position (\code{stop.pos}) that define the
range of the seed sequence.

Built-in options: \code{mer8},  \code{mer7A1}, \code{mer7m8}, \code{mer6}

Note: The seed definitions \code{mer8} and \code{mer7A1} force a U at position g1.
This results in an A in the target sequence being searched.
}
\examples{
# Example Ttr from Schlegel et al. 2022
guide.seq = "UUAUAGAGCAAGAACACUGUUUU"

# Get seed match
seed.seq = get_seed(guide.seq, "mer7m8")
}
