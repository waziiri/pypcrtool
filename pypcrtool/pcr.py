import matplotlib.pyplot as plt
import numpy as np

class InSilicoPCR:
    """
    A class to perform in silico PCR (Polymerase Chain Reaction) and visualize the results on a simulated gel electrophoresis.

    Attributes:
        forward_primer (str): The forward primer sequence.
        reverse_primer (str): The reverse primer sequence.
        sequence_file (str): The file containing the DNA sequence.
        sequence (str): The DNA sequence loaded from the file.
        forward_mismatch_tolerance (int): The allowed number of mismatches for the forward primer.
        reverse_mismatch_tolerance (int): The allowed number of mismatches for the reverse primer.
    """

    def __init__(self, forward_primer, reverse_primer, sequence_file, forward_mismatch_tolerance=0, reverse_mismatch_tolerance=0):
        """
        Initializes the InSilicoPCR object with primer sequences, mismatch tolerances, and DNA sequence file.
        
        Args:
            forward_primer (str): The forward primer sequence.
            reverse_primer (str): The reverse primer sequence.
            sequence_file (str): The file containing the DNA sequence.
            forward_mismatch_tolerance (int): Allowed mismatches for the forward primer.
            reverse_mismatch_tolerance (int): Allowed mismatches for the reverse primer.
        """
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.sequence_file = sequence_file
        self.sequence = self.load_sequence()
        self.forward_mismatch_tolerance = forward_mismatch_tolerance
        self.reverse_mismatch_tolerance = reverse_mismatch_tolerance
        
    def load_sequence(self):
        """
        Loads the DNA sequence from the file.

        Returns:
            str: The loaded DNA sequence.
        """
        sequence = ""
        with open(self.sequence_file, "r") as file:
            for line in file:
                if not line.startswith(">"):
                    sequence += line.strip()
        return sequence
    
    def reverse_complement(self, seq):
        """
        Computes the reverse complement of a DNA sequence.

        Args:
            seq (str): The DNA sequence.

        Returns:
            str: The reverse complement of the sequence.
        """
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement[base] for base in reversed(seq))
    
    def find_approximate_matches(self, primer, sequence, max_mismatches):
        """
        Finds positions in the sequence where the primer can bind with allowed mismatches.

        Args:
            primer (str): The primer sequence.
            sequence (str): The DNA sequence.
            max_mismatches (int): The allowed number of mismatches.

        Returns:
            list: List of positions where the primer can bind.
        """
        matches = []
        primer_len = len(primer)
        
        for i in range(len(sequence) - primer_len + 1):
            mismatches = sum(1 for a, b in zip(primer, sequence[i:i+primer_len]) if a != b)
            if mismatches <= max_mismatches:
                matches.append(i)
        
        return matches
    
    def find_primer_binding_sites(self):
        """
        Finds the binding sites of forward and reverse primers with allowed mismatches.

        Returns:
            tuple: Lists of binding sites for forward and reverse primers.
        """
        forward_binding_sites = self.find_approximate_matches(self.forward_primer, self.sequence, self.forward_mismatch_tolerance)
        reverse_primer_rc = self.reverse_complement(self.reverse_primer)
        reverse_binding_sites = self.find_approximate_matches(reverse_primer_rc, self.sequence, self.reverse_mismatch_tolerance)
        return forward_binding_sites, reverse_binding_sites
    
    def perform_pcr(self):
        """
        Simulates the PCR process by identifying product sequences formed by primer binding.

        Returns:
            list: List of tuples containing product sequence, size, start, and end positions.
        """
        forward_binding_sites, reverse_binding_sites = self.find_primer_binding_sites()
        products = []
        
        for f_site in forward_binding_sites:
            for r_site in reverse_binding_sites:
                if f_site < r_site:
                    product_seq = self.sequence[f_site:r_site + len(self.reverse_primer)]
                    products.append((product_seq, r_site + len(self.reverse_primer) - f_site, f_site, r_site + len(self.reverse_primer)))
        
        return products
    
    def check_primer_specificity(self):
        """
        Checks if the primers are specific (bind to one site only).

        Prints:
            Primer specificity information.
        """
        forward_occurrences = len(self.find_approximate_matches(self.forward_primer, self.sequence, self.forward_mismatch_tolerance))
        reverse_occurrences = len(self.find_approximate_matches(self.reverse_complement(self.reverse_primer), self.sequence, self.reverse_mismatch_tolerance))
        
        if forward_occurrences != 1 or reverse_occurrences != 1:
            print("Primers are not specific!")
        else:
            print("Primers are specific.")
    
    def print_products(self, products):
        """
        Prints the PCR product sequences to the console.

        Args:
            products (list): List of PCR products.
        """
        for i, (product, size, start, end) in enumerate(products):
            product_id = f"product_{i+1}"
            header = f">{product_id} size={size}bp start={start} end={end}"
            print(header)
            print(product)
        print(f"\nTotal PCR products: {len(products)}")

    def save_products(self, products, output_file):
        """
        Saves the PCR product sequences to a specified file.

        Args:
            products (list): List of PCR products.
            output_file (str): The file to save the PCR products.
        """
        with open(output_file, "w") as file:
            for i, (product, size, start, end) in enumerate(products):
                product_id = f"product_{i+1}"
                header = f">{product_id} size={size}bp start={start} end={end}"
                file.write(f"{header}\n")
                file.write(f"{product}\n")
        print(f"\nProducts saved to {output_file}")

    def visualize_gel(self, products, save_path=None):
        """
        Creates a simulated gel electrophoresis image of the PCR products.

        Args:
            products (list): List of PCR products.
            save_path (str, optional): Path to save the gel image.
        """
        sizes = [size for _, size, _, _ in products]
        sizes.sort()
        
        fig, ax = plt.subplots(figsize=(5, 10))

        ax.set_xlim(0, 100)
        max_size = max(sizes) if sizes else 1000
        ax.set_ylim(0, max_size + 100)
        ax.set_xticks([10, 30, 50, 70])
        ax.set_xticklabels(['M', '1', '2', '3'])
        ax.xaxis.tick_top()
        ax.set_xlabel('')
        ax.set_ylabel('Size (bp)')
        ax.set_title('Gel Electrophoresis of PCR Products')

        # Set y-ticks and labels to represent sizes from bottom to top
        tick_positions = np.linspace(0, max_size + 100, 10)
        ax.set_yticks(tick_positions)
        ax.set_yticklabels([f'{int(pos)}' for pos in tick_positions])

        # Create a gradient background to simulate a gel
        gradient = np.linspace(1, 0, 256)
        gradient = np.vstack((gradient, gradient))
        ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap('gray'), extent=[0, 100, 0, max_size + 100])

        # Plot bands in lane 2 only
        lane_position = 50
        for size in sizes:
            ax.plot([lane_position - 5, lane_position + 5], [size, size], color='black', lw=2)
            ax.text(lane_position - 15, size, f'{size} bp', va='center', ha='right', fontsize=8, color='black')

        if save_path:
            plt.savefig(save_path)
            print(f"Gel image saved to {save_path}")
        plt.show()

# Usage Example:
if __name__ == "__main__":
    forward_primer = "TCGAGAGGAACAGCCAAACT"
    reverse_primer = "TTCCTCATGTCCAGGTCCTC"
    sequence_file = "DV1_NSP5.fasta"  # Replace with your sequence file
    forward_mismatch_tolerance = 0  # Set mismatch tolerance for forward primer (0, 1, or 2)
    reverse_mismatch_tolerance = 0  # Set mismatch tolerance for reverse primer (0, 1, or 2)

    pcr_simulator = InSilicoPCR(forward_primer, reverse_primer, sequence_file, forward_mismatch_tolerance, reverse_mismatch_tolerance)
    pcr_simulator.check_primer_specificity()
    products = pcr_simulator.perform_pcr()
    
    if products:
        pcr_simulator.print_products(products)
        pcr_simulator.save_products(products, "pcr_products.fasta")
        pcr_simulator.visualize_gel(products, save_path="gel_image.png")
    else:
        print("No PCR products found.")
