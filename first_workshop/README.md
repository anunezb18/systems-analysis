# Artificial DNA Sequence Generator

This program generates a database of artificial DNA sequences with specified properties. The sequences are composed of the typical nucleotide bases (A, C, G, T) and are generated based on given probabilities for each base. The program also includes a filter to remove sequences with low entropy, ensuring a more diverse and chaotic set of sequences.

## Features

- Generates a specified number of DNA sequences of a given length.
- Allows setting the probabilities for each nucleotide base (A, C, G, T).
- Filters out sequences with low Shannon entropy to ensure diversity.
- Finds and reports the most frequent motif of a specified size in the generated sequences.
- Saves the generated sequences to a `.txt` file.

## Usage

### Parameters

- `numSequences`: Number of sequences to generate.
- `sequenceLength`: Length of each sequence.
- `probabilities`: List of probabilities for each nucleotide base (A, C, G, T).
- `motifSize`: Size of the motif to find in the sequences.
- `entropyThreshold`: Minimum entropy required for a sequence to be included.