"""
Test script for the PS-TEROS output formatter.

This script demonstrates the functionality of the output formatter
by creating a sample output file with mock calculation results.
"""

from teros.utils.output_formatter import create_output_file


def test_output_formatter():
    """Test the output formatter with sample data."""
    # Sample results data that might come from a real calculation
    sample_results = {
        "Formation Enthalpy": (2.345678, "eV/atom"),
        "Surface Energy (100)": (0.123456, "J/m²"),
        "Surface Energy (110)": (0.234567, "J/m²"),
        "Surface Energy (111)": (0.345678, "J/m²"),
        "Oxygen Chemical Potential": (-4.567890, "eV"),
        "Temperature": (298.15, "K"),
        "Pressure": (1.0, "atm"),
    }
    
    # Create the output file
    create_output_file(sample_results, "test_output.txt")
    print("Test output file 'test_output.txt' created successfully!")


if __name__ == "__main__":
    test_output_formatter()
