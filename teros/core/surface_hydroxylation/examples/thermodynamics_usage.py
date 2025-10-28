"""
Usage examples for JANAF thermodynamics database.

This script demonstrates how to use the JanafDatabase class for
chemical potential corrections in surface energy calculations.
"""

from teros.core.surface_hydroxylation import JanafDatabase


def example_basic_usage():
    """Basic usage: get corrections at room temperature."""
    print("=" * 60)
    print("Example 1: Basic Usage")
    print("=" * 60)

    db = JanafDatabase()

    # Get correction at room temperature (298 K, 1 bar)
    mu_h2o = db.get_mu_correction('H2O', T=298)
    mu_h2 = db.get_mu_correction('H2', T=298)
    mu_o2 = db.get_mu_correction('O2', T=298)

    print(f"\nChemical potential corrections at 298 K, 1 bar:")
    print(f"  μ(H2O) = {mu_h2o:.4f} eV")
    print(f"  μ(H2)  = {mu_h2:.4f} eV")
    print(f"  μ(O2)  = {mu_o2:.4f} eV")

    # For atomic oxygen (not molecular)
    mu_O_atom = mu_o2 / 2.0
    print(f"  μ(O atom) = {mu_O_atom:.4f} eV")


def example_temperature_dependence():
    """Demonstrate temperature dependence."""
    print("\n" + "=" * 60)
    print("Example 2: Temperature Dependence")
    print("=" * 60)

    db = JanafDatabase()

    temps = [298, 350, 500, 1000]
    print(f"\nμ(H2O) at different temperatures:")

    for T in temps:
        mu = db.get_mu_correction('H2O', T=T)
        print(f"  {T:4d} K: {mu:.4f} eV")


def example_pressure_correction():
    """Demonstrate pressure corrections."""
    print("\n" + "=" * 60)
    print("Example 3: Pressure Corrections")
    print("=" * 60)

    db = JanafDatabase()

    T = 350  # Near 100°C
    pressures = [0.1, 0.5, 1.0, 2.0]

    print(f"\nμ(H2O) at {T} K, different pressures:")

    for P in pressures:
        mu = db.get_mu_correction('H2O', T=T, P=P)
        print(f"  {P:4.1f} bar: {mu:.4f} eV")


def example_surface_energy_calculation():
    """Example: Calculate hydroxylation energy."""
    print("\n" + "=" * 60)
    print("Example 4: Surface Energy Calculation")
    print("=" * 60)

    db = JanafDatabase()

    # Example energies from VASP (dummy values)
    E_pristine = -1500.0  # eV
    E_hydroxylated = -1495.0  # eV
    n_OH = 2  # Added 2 OH groups

    # Experimental conditions
    T = 298  # K
    P_H2O = 0.023  # bar (room temperature water vapor)

    # Chemical potentials
    mu_h2o = db.get_mu_correction('H2O', T=T, P=P_H2O)
    mu_h2 = db.get_mu_correction('H2', T=T, P=1.0)

    # Formation reaction: Surface + n*H2O -> Surface-nOH + n/2*H2
    # ΔE = E_hydroxylated - E_pristine - n*μ(H2O) + (n/2)*μ(H2)
    delta_E = E_hydroxylated - E_pristine - n_OH * mu_h2o + (n_OH / 2) * mu_h2

    print(f"\nHydroxylation calculation:")
    print(f"  Temperature: {T} K")
    print(f"  P(H2O): {P_H2O} bar")
    print(f"  Number of OH groups: {n_OH}")
    print(f"  μ(H2O): {mu_h2o:.4f} eV")
    print(f"  μ(H2): {mu_h2:.4f} eV")
    print(f"  ΔE_formation: {delta_E:.4f} eV")


def example_raw_data_access():
    """Demonstrate raw data access for verification."""
    print("\n" + "=" * 60)
    print("Example 5: Raw Data Access")
    print("=" * 60)

    db = JanafDatabase()

    raw = db.get_raw_data('H2O', T=298)

    print(f"\nRaw JANAF data for H2O at 298 K:")
    print(f"  H(T) - H(0K): {raw['H']:.3f} kJ/mol")
    print(f"  S(T): {raw['S']:.3f} J/(mol·K)")
    print(f"  Δμ⁰: {raw['delta_mu']:.4f} eV")

    # Show available temperatures
    temps = db.list_temperatures()
    print(f"\nAvailable temperature points: {len(temps)}")
    print(f"  Range: {temps[0]}-{temps[-1]} K")
    print(f"  Step: {temps[1] - temps[0]} K")


if __name__ == '__main__':
    example_basic_usage()
    example_temperature_dependence()
    example_pressure_correction()
    example_surface_energy_calculation()
    example_raw_data_access()

    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)
