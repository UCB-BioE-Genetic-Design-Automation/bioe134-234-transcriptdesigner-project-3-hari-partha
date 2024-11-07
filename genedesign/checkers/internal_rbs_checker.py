class InternalRBSChecker:
    """
    Identifies potential internal ribosome binding sites (RBS), such as Shine-Dalgarno sequences 
    followed by a start codon (ATG, GTG, TTG), within a DNA sequence to prevent unintended translation initiation.
    """

    def __init__(self):
        self.shine_dalgarno_motifs = []
        self.start_codons = []

    def initiate(self):
        """
        Sets up the internal RBS checker by initializing Shine-Dalgarno motifs and start codons.
        """
        self.shine_dalgarno_motifs = ["AGGAGG", "GGAGG"]  # Common Shine-Dalgarno sequences
        self.start_codons = ["ATG", "GTG", "TTG"]

    def run(self, dna_sequence):
        """
        Checks for internal RBS (Shine-Dalgarno sequence + start codon) in a given DNA sequence.

        Parameters:
            dna_sequence (str): DNA sequence to check for internal RBS.

        Returns:
            tuple: (bool, str or None)
                - True and None if no internal RBS is found.
                - False and the problematic sequence if an internal RBS is detected.
        """
        dna_sequence = dna_sequence.upper()
        
        # Search for each Shine-Dalgarno motif in the sequence
        for motif in self.shine_dalgarno_motifs:
            position = dna_sequence.find(motif)
            
            # If motif is found, check for a start codon within 5-10 bases downstream
            while position != -1:
                # Define the search window for the start codon (5-10 bases downstream of the Shine-Dalgarno motif)
                start_search_position = position + len(motif) + 5
                end_search_position = position + len(motif) + 10
                downstream_sequence = dna_sequence[start_search_position:end_search_position]
                
                for codon in self.start_codons:
                    if codon in downstream_sequence:
                        # Internal RBS found, return False with the problematic sequence
                        return False, dna_sequence[position:start_search_position + 3]
                
                # Continue searching for additional Shine-Dalgarno motifs in the sequence
                position = dna_sequence.find(motif, position + 1)

        # If no internal RBS is detected, return True
        return True, None

if __name__ == "__main__":
    checker = InternalRBSChecker()
    checker.initiate()  # Initialize the checker

    # Example test sequences
    result = checker.run("AAAGGAGGTAGGGGTGATGAAA")
    print(result)  # Expected output: False with problematic sequence

    result = checker.run("TTTCCCGTGGGCACTGAGCACTG")
    print(result)  # Expected output: True (no internal RBS detected)
