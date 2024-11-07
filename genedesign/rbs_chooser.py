from genedesign.models.rbs_option import RBSOption
from genedesign.seq_utils.calc_edit_distance import calculate_edit_distance
from genedesign.seq_utils.hairpin_counter import hairpin_counter
from genedesign.seq_utils.Translate import Translate
import pandas as pd  # type:ignore

df = pd.read_csv('/Users/haripartha/Documents/bioe134-234-transcriptdesigner-project-3-hari-partha/genedesign/data/top_5_percent_gene_data')

from typing import Set

class RBSChooser:
    """
    A class to choose the best RBS for a given CDS sequence.
    """

    rbs_options: Set[RBSOption] = set()
    translator: Translate = Translate()

    def initiate(self) -> None:
        """
        Initialization method for RBSChooser.
        """
        self.translator.initiate()

        # Populate RBS options from the provided dataset
        for _, value in df.iterrows():
            gene_name = value['gene']
            utr = value['UTR']
            cds = value['CDS']
            first_six_aas = self.translator.run(cds)[:6]

            rbs_option = RBSOption(
                utr=utr,
                cds=cds,
                gene_name=gene_name,
                first_six_aas=first_six_aas
            )
            self.rbs_options.add(rbs_option)

    def run(self, cds: str, ignores: Set[RBSOption]) -> RBSOption:
        """
        Executes the RBS selection process for the given CDS.

        Parameters:
        - cds (str): The coding sequence to pair with an RBS.
        - ignores (Set[RBSOption]): A set of RBSOption instances to ignore during selection.

        Returns:
        - RBSOption: The selected RBSOption that best pairs with the given CDS.
        """
        # Exclude RBS options in the ignore set
        valid_rbs_options = [rbs for rbs in self.rbs_options if rbs not in ignores]

        # Check if valid options remain after exclusion
        if not valid_rbs_options:
            raise ValueError("No valid RBS options remaining after exclusion.")

        best_rbs = None
        best_score = float("inf")  # Lower scores are better
        fallback_rbs = None
        fallback_score = float("inf")  # Fallback option with lowest peptide edit distance

        for rbs in valid_rbs_options:
            # Combine UTR and CDS sequences to evaluate hairpin structure
            combined_sequence = rbs.utr + cds
            hairpin_count = hairpin_counter(combined_sequence)[0]

            # Calculate peptide similarity (edit distance) with the translated CDS
            translated_input_peptide = self.translator.run(cds)[:6]
            peptide_edit_distance = calculate_edit_distance(translated_input_peptide, rbs.first_six_aas)

            # Fallback: Track the RBS with the lowest peptide edit distance for use if no RBS meets all criteria
            if peptide_edit_distance < fallback_score:
                fallback_score = peptide_edit_distance
                fallback_rbs = rbs

            # Primary scoring criteria: lower hairpin count and closest peptide match
            if hairpin_count <= 4 and peptide_edit_distance < best_score:
                best_score = peptide_edit_distance
                best_rbs = rbs

        # Return the best RBS if found; otherwise, use the fallback RBS with the lowest peptide edit distance
        return best_rbs if best_rbs is not None else fallback_rbs
