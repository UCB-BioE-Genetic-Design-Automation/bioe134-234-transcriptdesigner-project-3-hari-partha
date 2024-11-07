import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript

from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker  # Import your internal RBS checker

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence using a hybrid 
    Guided Random + Sliding Window approach for optimal codon selection satisfying
    high CAI, low hairpin count, and the absence of internal promoters & forbidden sequences.
    """

    def __init__(self):

        self.codon_usage = {}
        self.rbsChooser = None

        self.codonChecker = None
        self.forbiddenSequenceChecker = None
        self.promoterChecker = None
        self.internalRBSChecker = None  # Internal RBS checker instead of RNaseE
    
    
    def initiate(self):
        """
        Initialization method for RBSChooser and checkers.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.forbiddenSequenceChecker = ForbiddenSequenceChecker()
        self.promoterChecker = PromoterChecker()
        self.codonChecker = CodonChecker()
        self.internalRBSChecker = InternalRBSChecker()  # Use InternalRBSChecker 
        
        # Initialize checkers
        self.forbiddenSequenceChecker.initiate()
        self.promoterChecker.initiate()
        self.codonChecker.initiate()
        self.internalRBSChecker.initiate()

        # Load codon usage data from the provided text file
        self.codon_usage = self.load_codon_usage('/Users/haripartha/Documents/bioe134-234-transcriptdesigner-project-3-hari-partha/genedesign/data/codon_usage.txt')

    
    def load_codon_usage(self, filepath: str) -> dict:
        """
        Parses a codon usage file and returns a dictionary mapping amino acids to their codons and frequencies.
        
        Parameters:
            filepath (str): Path to the codon usage file.
        
        Returns:
            dict: A dictionary where keys are amino acids and values are lists of tuples (codon, frequency).
        """
        codon_usage = {}

        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
            
                if len(parts) >= 3:
                    codon = parts[0].strip()  # Codon (e.g., TTT)
                    aa = parts[1].strip()     # Amino acid (e.g., F)
                    frequency = float(parts[2])  # Frequency (e.g., 0.58)

                    # Add codon and frequency to the corresponding amino acid entry
                    if aa not in codon_usage:
                        codon_usage[aa] = []
                    codon_usage[aa].append((codon, frequency))
        
        return codon_usage

    
    def guided_random_codon(self, aa: str) -> str:
        """
        Selects a codon for an amino acid using guided random selection based on frequency.
        
        Parameters:
            aa (str): Amino acid single-letter code
        
        Returns:
            str: Selected codon
        """
        codons = self.codon_usage.get(aa)
        
        if not codons:
            raise ValueError(f"No codons available for amino acid {aa}")
        
        # Normalize frequencies and select based on random value
        total_freq = sum(freq for _, freq in codons)
        rand_val = random.uniform(0, total_freq)
        
        cumulative_freq = 0
        for codon, freq in codons:
            cumulative_freq += freq
            if rand_val <= cumulative_freq:
                return codon
        
        return codons[-1][0]  # Fallback to last option
   
    
    def score_candidates(self, candidates):
        """
        Scores candidate solutions based on various criteria like forbidden sequences,
        secondary structure formation, internal RBS sites, etc.
        
        Parameters:
            candidates (List[List[str]]): List of candidate solutions (codons).
        
        Returns:
            List[str]: Best candidate solution that passes the most checks.
        """
        
        scored_candidates = []

        # Define weights for each checker 
        weights = {
            "forbidden": 6,
            "hairpin": 4,
            "promoter": 1,
            "internal_rbs": 1,
            "codon_usage": 4
        }

        for candidate in candidates:
            dna_seq = ''.join(candidate)
            score = 0

            # Run all checkers on the candidate DNA sequence and accumulate scores
            if self.forbiddenSequenceChecker.run(dna_seq):
                score += weights["forbidden"]

            if hairpin_checker(dna_seq):
                score += weights["hairpin"]

            if self.promoterChecker.run(dna_seq):
                score += weights["promoter"]

            if self.internalRBSChecker.run(dna_seq):
                score += weights["internal_rbs"]

            if self.codonChecker.run(candidate):
                score += weights["codon_usage"]

            # Append the candidate and its score to the list
            scored_candidates.append((candidate, score))

        # Sort candidates by their scores in descending order
        scored_candidates.sort(key=lambda x: x[1], reverse=True)

        # Return the candidate with the highest score
        best_candidate = scored_candidates[0][0] if scored_candidates else candidates[0]
        
        return best_candidate
    
    
    def validate_window(self, candidate):
        """
        Validates a candidate solution by running it through all checkers.
        
        Parameters:
            candidate (List[str]): A list of codons representing a potential solution.
        
        Returns:
            bool: True if the candidate passes all checks; False otherwise.
        """
        
        dna_seq = ''.join(candidate)

        # Run all checkers on the candidate DNA sequence
        
        if not self.forbiddenSequenceChecker.run(dna_seq):
            return False

        if not hairpin_checker(dna_seq):
            return False

        if not self.promoterChecker.run(dna_seq):
            return False

        if not self.internalRBSChecker.run(dna_seq):
            return False

        if not self.codonChecker.run(candidate):
            return False

        # If all checks pass
        return True
    
    
    def sliding_window_optimization(self, peptide: str) -> str:
        """
        Optimizes peptide translation using sliding window approach combined with guided random selection.
        
        Parameters:
            peptide (str): The protein sequence
        
        Returns:
            str: Optimized DNA coding sequence including dynamically growing preamble.
        """
        
        preamble = []
        cds = []
        window_size = 3

        for i in range(0, len(peptide), window_size):
            window_peptide = peptide[i:i + window_size]
            downstream_peptide = peptide[i + window_size:i + window_size + 6]

            # Generate multiple possible solutions using guided random selection
            candidate_codons = [
                [self.guided_random_codon(aa) for aa in window_peptide]
                for _ in range(10)
            ]

            # Validate each candidate until we find one that passes all checks
            best_candidate = None
            for candidate in candidate_codons:
                if self.validate_window(candidate):
                    best_candidate = candidate
                    break
            
            # If no valid candidates are found after validation retries, get candidate with highest score
            if best_candidate is None:
                best_candidate = self.score_candidates(candidate_codons)

            # Retain only middle part of this candidate (for overlap), or all if it's at the end of the sequence
            if len(window_peptide) == window_size:
                cds.extend(best_candidate[:window_size])
            else:
                cds.extend(best_candidate[:len(window_peptide)])

            # Update preamble with newly optimized middle part of this candidate 
            preamble += best_candidate[:window_size]

        return ''.join(cds)
    
    
    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA using hybrid algorithm and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with selected RBS and translated codons.
        """
        
        # Optimize CDS using sliding window + guided random approach
        cds_sequence = self.sliding_window_optimization(peptide)

        # Append stop codon (TAA)
        cds_sequence += "TAA"

        # Choose an RBS using RBSChooser while ignoring specified options
        selected_rbs = self.rbsChooser.run(cds_sequence, ignores)

        # Return transcript object with selected RBS and translated CDS as a list of codons
        return Transcript(selected_rbs, peptide, cds_sequence.split())

if __name__ == "__main__":
    peptide = "MYPFIRTARMTV"

    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)

    print(transcript)
