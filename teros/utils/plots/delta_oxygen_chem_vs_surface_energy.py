import matplotlib.pyplot as plt
from aiida import load_profile
from aiida.orm import load_node

# Load AiiDA profile (skip if running in verdi shell)
try:
    load_profile()
except Exception:
    pass

MAIN_NODE_PK = 16275  # Change if needed

# Load main node
main_node = load_node(MAIN_NODE_PK)

# Output keys for surface terminations
surface_keys = ['s_0', 's_1', 's_2', 's_3']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

plt.figure(figsize=(8, 6))

for key, color in zip(surface_keys, colors):
    output_dict = main_node.outputs[key].get_dict()
    gamma_dict = output_dict['gamma_values_fixed_metals']
    # Convert keys (mu_O) to float and sort
    mu_o = [float(mu) for mu in gamma_dict.keys()]
    gamma = [gamma_dict[k] for k in gamma_dict.keys()]
    # Sort by mu_O
    mu_o, gamma = zip(*sorted(zip(mu_o, gamma)))
    plt.plot(mu_o, gamma, marker='o', label=key, color=color)

plt.xlabel(r'$\Delta \mu_O$ (eV)')
plt.ylabel(r'$\gamma$ (eV/$\mathrm{\AA}^2$)')
plt.title('Surface Energies vs. Oxygen Chemical Potential')
plt.legend(title='Surface')
plt.tight_layout()
plt.show()