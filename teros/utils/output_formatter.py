"""
Output formatter for PS-TEROS calculations using Rich library.

This module provides functions to create formatted output files
with a professional appearance including headers, tables, and results.
"""

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich import box
import datetime
import teros


def generate_front_page():
    """
    Create a formatted front page with PS-TEROS logo, version info, and description.
    
    Returns:
        str: Formatted text representation of the front page
    """
    # Create a console that records output
    console = Console(record=True, width=80)
    
    # Create the title with ASCII art style using Rich Text
    title = Text("PS-TEROS", style="bold blue")
    subtitle = Text("Ab Initio Atomistic Thermodynamics for Surface Energetics", 
                   style="italic bright_blue")
    version = Text(f"Version {teros.__version__}", style="dim")
    description = Text(
        "A package for calculating surface Gibbs free energy using ab initio atomistic\n"
        "thermodynamics with AiiDA Workgraph framework.", 
        style="white"
    )
    
    # Create the main panel with all information
    panel_content = Text()
    panel_content.append(title, style="bold blue")
    panel_content.append("\n")
    panel_content.append(subtitle, style="italic bright_blue")
    panel_content.append("\n")
    panel_content.append(version, style="dim")
    panel_content.append("\n\n")
    panel_content.append(description, style="white")
    
    # Add timestamp
    timestamp = Text(f"\n\nGenerated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", 
                    style="dim")
    panel_content.append(timestamp)
    
    # Print the panel to the console
    console.print(Panel(
        panel_content,
        title="[bold]PS-TEROS Output Report[/bold]",
        subtitle="[italic]Surface Thermodynamics Calculation Results[/italic]",
        border_style="bright_blue",
        box=box.ROUNDED,
        padding=(1, 2)
    ))
    
    # Return the recorded text
    return console.export_text()


def create_results_table(results_data):
    """
    Create a formatted table of results.
    
    Args:
        results_data (dict): Dictionary containing calculation results
        
    Returns:
        Table: Rich Table object with formatted results
    """
    table = Table(title="Calculation Results", box=box.ROUNDED, style="blue")
    
    # Add columns
    table.add_column("Property", style="cyan", no_wrap=True)
    table.add_column("Value", style="magenta", justify="right")
    table.add_column("Units", style="green")
    
    # Add data rows
    for property_name, (value, units) in results_data.items():
        table.add_row(property_name, f"{value:.6f}", units)
    
    return table


def create_output_file(results_data, filename="PS_TEROS_OUTPUT.txt"):
    """
    Create a formatted output file with front page and results.
    
    Args:
        results_data (dict): Dictionary containing calculation results
                            Format: {'property_name': (value, 'units')}
        filename (str): Name of the output file
    """
    # Create the output file
    with open(filename, "w") as f:
        # Create a console that writes to the file
        console = Console(file=f, width=80)
        
        # Write the front page
        front_page_text = generate_front_page()
        f.write(front_page_text)
        f.write("\n\n")
        
        # Write the results table
        table = create_results_table(results_data)
        console.print(table)
        
        # Add completion message
        console.print("\n[bold green]Calculation completed successfully![/]")
        console.print("[dim]End of report[/dim]")


# Example usage function
def example_usage():
    """Example of how to use the output formatter"""
    # Sample results data
    sample_results = {
        "Formation Enthalpy": (2.345678, "eV/atom"),
        "Surface Energy (100)": (0.123456, "J/m²"),
        "Surface Energy (110)": (0.234567, "J/m²"),
        "Surface Energy (111)": (0.345678, "J/m²"),
        "Oxygen Chemical Potential": (-4.567890, "eV"),
    }
    
    # Create the output file
    create_output_file(sample_results, "example_output.txt")
    print(f"Example output file 'example_output.txt' created successfully!")


if __name__ == "__main__":
    example_usage()
