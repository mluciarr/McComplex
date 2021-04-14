import os
import argparse

class InputError(FileNotFoundError):
    """Checks if the input variable is a directory."""

    def __init__(self, input_dir):
        self.input_dir = input_dir

    def __str__(self):
        return "%s is not a directory. Please, try again with a proper directory." % (self.input_dir)

