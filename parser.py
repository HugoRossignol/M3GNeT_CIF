import numpy as np
import re

class CIFParser:
    def __init__(self, cif_file):
        """
        Initialize a CIFParser instance.

        Args:
            cif_file (str): The path to the CIF file to be parsed.
        """
        self.symmetry_operations = []
        self.atomic_positions = {}
        self.lattice_parameters = []
        self.pos_info = []
        self.index = 0

        with open(cif_file, 'r') as file:
            lines = file.read().splitlines()

        pos_block = False
        symmetry_block = False
        for line in lines:
            if line.startswith("_space_group_symop_operation_xyz") or line.startswith("_symmetry_equiv_pos_as_xyz"):
                symmetry_block = True
                continue
            if line.startswith("_atom_site_label"):
                pos_block = True
            if line.startswith("loop_") or line == "":
                symmetry_block = False
                pos_block = False
            if pos_block:
                if "_" in line:
                    self.parse_pos_info(line)
                    self.index += 1
                else:
                    self.parse_pos(line)
            if symmetry_block:
                self.parse_symmetry(line)
            elif line.startswith("_cell_length_a"):
                self.lattice_parameters.append(float(line.split()[-1].split("(")[0]))
            elif line.startswith("_cell_length_b"):
                self.lattice_parameters.append(float(line.split()[-1].split("(")[0]))
            elif line.startswith("_cell_length_c"):
                self.lattice_parameters.append(float(line.split()[-1].split("(")[0]))
            elif line.startswith("_cell_angle_alpha"):
                self.lattice_parameters.append(float(line.split()[-1].split("(")[0]))
            elif line.startswith("_cell_angle_beta"):
                self.lattice_parameters.append(float(line.split()[-1].split("(")[0]))
            elif line.startswith("_cell_angle_gamma"):
                self.lattice_parameters.append(float(line.split()[-1].split("(")[0]))
        self.apply_symmetries()

    def parse_symmetry(self, line):
        """
        Parse symmetry operations from CIF file and add them to the list of symmetry operations.

        Args:
            line (str): A line from the CIF file.
        """

        # Use a regular expression to match and extract the coordinates
        match = re.search(r'(?:(^\d+\s)|(?<=\s))?[\'"]?(-?[a-zA-Z\d./+*()-]+,\s?-?[a-zA-Z\d./+*()-]+,\s?-?[a-zA-Z\d./+*()-]+)[\'"]?',line)

        if match:
            self.symmetry_operations.append(match.group(2).replace("'","")) 
        else:
            raise ValueError("Positions in CIF cannot be parsed.")

    def parse_pos_info(self, line):
        """
        Parse position information from CIF file.

        Args:
            line (str): A line from the CIF file.
        """
        if line.startswith("_atom_site_label") or line.startswith("_atom_site_fract"):
            self.pos_info.append(self.index)

    def parse_pos(self, line):
        """
        Parse atomic positions from CIF file and store them in a dictionary.

        Args:
            line (str): A line from the CIF file.
        """
        line_list = line.split()
        symbol, x, y, z = [line_list[i].split("(")[0] for i in self.pos_info]
        
        match = re.match(r'^[A-Za-z]+', symbol)

        if match:
            symbol = match.group()
        else:
            raise ElementFormatError("Element is not formatted properly")


        if symbol in self.atomic_positions:
            self.atomic_positions[symbol] = np.vstack((self.atomic_positions[symbol], np.asarray([float(x), float(y), float(z)])))
        else:
            self.atomic_positions[symbol] = np.asarray([[float(x), float(y), float(z)]])

    def apply_symmetries(self):
        """
        Apply symmetry operations to generate unique atomic positions.
        """
        self.unique_positions = {}
        for symbol, positions in self.atomic_positions.items():
            self.unique_positions[symbol] = []
            for symm_op in self.symmetry_operations:
                for position in positions:
                    new_pos = self.apply_symmetry_operator(symm_op, position)
                    if not any(np.allclose(new_pos, existing_pos) for existing_pos in self.unique_positions[symbol]):
                        self.unique_positions[symbol].append(new_pos)

        for s in self.unique_positions:
            self.unique_positions[s] = np.asarray(self.unique_positions[s])

    def apply_symmetry_operator(self, symm_op, position):
        """
        Apply a symmetry operator to a given position.

        Args:
            symm_op (str): The symmetry operation to apply.
            position (numpy.ndarray): The position to which the symmetry operation should be applied.

        Returns:
            numpy.ndarray: The new position after applying the symmetry operation.
        """
        x, y, z = position
        new_x = eval(symm_op.split(',')[0]) % 1
        new_y = eval(symm_op.split(',')[1]) % 1
        new_z = eval(symm_op.split(',')[2]) % 1
        return np.array([new_x, new_y, new_z])

# Example usage:
if __name__ == "__main__":
    cif_parser = CIFParser("CIF/Test/Si.cif")
    lattice_parameters = cif_parser.lattice_parameters

    print("Lattice Parameters:")
    print(lattice_parameters)

    print("Unique Positions:")
    for symbol, positions in cif_parser.unique_positions.items():
        print(f"Element: {symbol}")
        print(positions)
