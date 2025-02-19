from Bio import SeqIO
import Levenshtein
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager
import csv
from tqdm import tqdm


def find_fuzzy_matches_in_record_all(record, query_kmers, k, max_mismatches):
    """Find k-mer matches with a maximum mismatch allowance in a single transcript record."""
    transcript_seq = str(record.seq)
    local_matches = {}

    for kmer in query_kmers:
        for i in range(len(transcript_seq) - k + 1):
            transcript_kmer = transcript_seq[i:i + k]
            mismatch_count = Levenshtein.distance(kmer, transcript_kmer)
            if mismatch_count <= max_mismatches:
                if kmer not in local_matches:
                    local_matches[kmer] = []
                local_matches[kmer].append((record.id, i, mismatch_count))

    return local_matches

def substitution_only_distance(seq1, seq2):
    """Calculate the number of substitutions (ignoring insertions and deletions) between two equal-length sequences."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length to calculate substitution-only distance.")
    return sum(1 for a, b in zip(seq1, seq2) if a != b)


def find_fuzzy_matches_in_record_subonly(record, query_kmers, k, max_mismatches):
    """Find k-mer matches with substitution-only mismatches in a single transcript record."""
    transcript_seq = str(record.seq)
    local_matches = {}

    for kmer in query_kmers:
        for i in range(len(transcript_seq) - k + 1):
            transcript_kmer = transcript_seq[i:i + k]
            mismatch_count = substitution_only_distance(kmer, transcript_kmer)
            if mismatch_count <= max_mismatches:
                if kmer not in local_matches:
                    local_matches[kmer] = []
                local_matches[kmer].append((record.id, i, mismatch_count))

    return local_matches


def find_fuzzy_kmer_matches_parallel(query_seq, transcriptome_file, k=21, max_mismatches=1, max_processes=20, subonly=True):
    """Find k-mer matches with multiprocessing, collecting results from each future as they complete."""
    # Generate unique k-mers from the query sequence
    query_kmers = {query_seq[i:i + k] for i in range(len(query_seq) - k + 1)}

    # Shared dictionary for final aggregated results
    manager = Manager()
    matches = manager.dict()  # Use a managed dictionary

    # Using ProcessPoolExecutor to submit tasks
    with ProcessPoolExecutor(max_workers=max_processes) as executor:
        records = list(SeqIO.parse(transcriptome_file, "fasta"))
        total_records = len(records)

        if subonly:
            with tqdm(total=total_records, desc="Processing transcripts") as pbar:
                futures = [
                    executor.submit(find_fuzzy_matches_in_record_subonly, record, query_kmers, k, max_mismatches)
                    for record in records
                ]
                # Collect results as each future completes
                for future in as_completed(futures):
                    local_matches = future.result()  # Get result from each completed future
                    # Aggregate local matches into the shared matches dictionary
                    for kmer, info in local_matches.items():
                        if kmer not in matches:
                            matches[kmer] = manager.list()
                        matches[kmer].extend(info)

                    # Update the progress bar
                    pbar.update(1)
        else:
            with tqdm(total=total_records, desc="Processing transcripts") as pbar:
                futures = [
                    executor.submit(find_fuzzy_matches_in_record_all, record, query_kmers, k, max_mismatches)
                    for record in records
                ]

                for future in as_completed(futures):
                    local_matches = future.result()  # Get result from each completed future
                    # Aggregate local matches into the shared matches dictionary
                    for kmer, info in local_matches.items():
                        if kmer not in matches:
                            matches[kmer] = manager.list()
                        matches[kmer].extend(info)

                    # Update the progress bar
                    pbar.update(1)
    # Convert managed dict to a regular dict for final output
    return {kmer: list(info) for kmer, info in matches.items() if info}


def save_matches_to_csv(kmer_matches, output_file):
    """Save k-mer matches to a CSV file with mismatch information."""
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(["K-mer", "Transcript ID", "Position", "Mismatches"])

        # Write each k-mer match entry
        for kmer, matches in kmer_matches.items():
            for transcript_id, position, mismatches in matches:
                writer.writerow([kmer, transcript_id, position, mismatches])


if __name__ == '__main__':
    # Define the input query sequence and the transcriptome file path
    # query_sequence = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGC" \
    #                  "GATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCG" \
    #                  "CTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGT" \
    #                  "GAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTA" \
    #                  "TATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGAC" \
    #                  "GGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATC" \
    #                  "ACTCTCGGCATGGATGAGCTGTACAAG"  # Replace with your 700 bp sequence
    query_sequence = "GTCCAATTGTGGAGAGGTTGTTAGCTATTGTGAAGTGGTATACAGTGAAGTGTGGTAAATAAAGAGACACTTGAATGACCCTTGGTTATTGATAAGAATGTAAAGAGGGTTTATACTGTG" \
                     "TGGTAGTTTTGTTTTTTGTGGGTCAGCGTAAAAGAATTTATCCTTTGAGATGTGTTATGGAGGTAGAAAAGTGTAGGAAGTACAGTCCCTTAATATTTAGAATACAAGAGGAAATGATATAGTG" \
                     "GAGACATGTAATGTGTGAGAGAGTGTGTACGAAGTAGCTGGTCTGTGCGTGATGTGTAGGTTTATGAAATTGGAGTATTGTGGTAATTAGTAATGATGAAGGTAATGAAGTAGTGATGTTCTGG" \
                     "TTTTCGAAAGAATGGAGTAAGGGATGGTGTATGG"
    transcriptome_file_path = "gregaria_transcriptome.fna"  # Replace with your file path

    # Run the parallelized k-mer comparison with k=20 and max mismatches of 3
    kmer_matches = find_fuzzy_kmer_matches_parallel(query_sequence, transcriptome_file_path, k=21, max_mismatches=4,
                                                    max_processes=20, subonly=True)

    # Save results to CSV
    output_file_path = "kmer_matches_k21_mm4_sub_trinity.csv"
    save_matches_to_csv(kmer_matches, output_file_path)
    print(f"Results saved to {output_file_path}")

    # Display the results
    for kmer, info in kmer_matches.items():
        for transcript_id, position, mismatches in info:
            print(
                f"K-mer: {kmer} found in transcript {transcript_id} at position {position} with {mismatches} mismatches")
