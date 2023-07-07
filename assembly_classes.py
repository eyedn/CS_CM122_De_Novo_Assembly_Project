###############################################################################
#   Aydin Karatas
#   CS CM122 Project 3
#   project3_classes.py
###############################################################################

class Read:
    def __init__(self, label: str) -> None:
        self.label = label
        self.sequence = ""

    def add_to_sequence(self, sequence: str) -> None:
        self.sequence += sequence