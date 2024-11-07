from genedesign.models.rbs_option import RBSOption
from genedesign.seq_utils.calc_edit_distance import calculate_edit_distance
from genedesign.seq_utils.hairpin_counter import hairpin_counter
from genedesign.seq_utils.Translate import Translate
import pandas as pd #type:ignore

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

        for key, value in df.iterrows(): #what to do with the top 5 percent gene data items
          gene_name = value['gene']
          utr = value['UTR']
          cds = value['CDS']
          first_six_aas = self.translator.run(cds)[:6]

          rbs_option = RBSOption(
                utr = utr,
                cds= cds,
                gene_name= gene_name,
                first_six_aas = first_six_aas
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
        # Step 1: Exclude the RBS options in the ignores set
        valid_rbs_options = [rbs for rbs in self.rbs_options if rbs not in ignores]

        # Check if valid options remain after exclusion
        if not valid_rbs_options:
            return None  # No valid RBS options left

        best_rbs = None
        best_score = float("inf")  # Lower scores are better

        for rbs in valid_rbs_options:
            # Step 2: Check for secondary structure occlusion (hairpin count)
            combined_sequence = rbs.utr + cds  # UTR + input CDS
            hairpin_count = hairpin_counter(combined_sequence)[0]

            # Step 3: Compare peptide similarity (edit distance)
            translated_input_peptide = self.translator.run(cds)[:6]  # Translate first 6 amino acids
            peptide_edit_distance = calculate_edit_distance(translated_input_peptide, rbs.first_six_aas)

            # Combine both criteria: fewer hairpins and smaller peptide distance
            if hairpin_count <= 4: # four was the smallest amount of hairpin counts I could identify
              if peptide_edit_distance < best_score: # I want to sort by the smallest peptide edit distance
                best_score = peptide_edit_distance
                best_rbs = rbs
            
        # if best_rbs is None:
        #     print(f"Warning: No valid RBS found for peptide {peptide}. Defaulting to empty RBSOption.")
        #     selectedRBS = RBSOption("", cds)  # Default to an empty RBS sequence if none is found

        return best_rbs

# Instantiate RBSChooser and initiate
chooser = RBSChooser()
chooser.initiate()

# Test CDS
cds = "ATGGCTAGCAAATACGATTTTACAATATAA"

# Initialize empty ignores set
ignores = set()

# First run, no ignored RBS options
selected_rbs_1 = chooser.run(cds, ignores)
print(f"Selected RBS: {selected_rbs_1}")

# Repeat to confirm determinacy
selected_rbs_2 = chooser.run(cds, ignores)
print(f"Selected RBS: {selected_rbs_2}")

# Add the returned RBSOption to ignores and run again, confirm a different selected RBS
ignores.add(selected_rbs_1)
selected_rbs_3 = chooser.run(cds, ignores)
print(f"Selected RBS after ignoring {selected_rbs_1}: {selected_rbs_3}")

