import os
import re

ROOT_DIR = os.path.expanduser("~/git/quantum-lego")

# Order matters: more specific patterns first
REPLACEMENTS = [
    # Specific function extraction from slabs
    (r"from teros\.core\.slabs import extract_total_energy", "from quantum_lego.core.common.utils import extract_total_energy"),
    
    # Internal dependencies -> quantum_lego.core.common.*
    (r"from teros\.core\.utils", "from quantum_lego.core.common.utils"),
    (r"import teros\.core\.utils", "import quantum_lego.core.common.utils"),
    
    (r"from teros\.core\.constants", "from quantum_lego.core.common.constants"),
    (r"import teros\.core\.constants", "import quantum_lego.core.common.constants"),
    
    (r"from teros\.core\.fixed_atoms", "from quantum_lego.core.common.fixed_atoms"),
    (r"import teros\.core\.fixed_atoms", "import quantum_lego.core.common.fixed_atoms"),
    
    (r"from teros\.core\.aimd_functions", "from quantum_lego.core.common.aimd_functions"),
    (r"import teros\.core\.aimd_functions", "import quantum_lego.core.common.aimd_functions"),
    
    (r"from teros\.core\.aimd", "from quantum_lego.core.common.aimd"),
    (r"import teros\.core\.aimd", "import quantum_lego.core.common.aimd"),
    
    (r"from teros\.core\.convergence", "from quantum_lego.core.common.convergence"),
    (r"import teros\.core\.convergence", "import quantum_lego.core.common.convergence"),
    
    (r"from teros\.core\.u_calculation", "from quantum_lego.core.common.u_calculation"),
    (r"import teros\.core\.u_calculation", "import quantum_lego.core.common.u_calculation"),

    # General package rename
    (r"from teros\.core\.lego", "from quantum_lego"),
    (r"import teros\.core\.lego", "import quantum_lego"),
    
    # Internal package references within the package itself (if any absolute imports remain)
    # This might overlap with the above, but let's be safe. 
    # Since we are moving teros/core/lego -> quantum_lego/core, 
    # imports like 'from . import ...' are fine.
    # But 'from teros.core.lego.bricks' should become 'from quantum_lego.core.bricks'
]

def process_file(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError:
        print(f"Skipping binary file: {filepath}")
        return
    
    new_content = content
    for pattern, replacement in REPLACEMENTS:
        new_content = re.sub(pattern, replacement, new_content)
        
    if new_content != content:
        print(f"Updating {filepath}")
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)

def main():
    print(f"Scanning {ROOT_DIR}...")
    for root, dirs, files in os.walk(ROOT_DIR):
        if '.git' in dirs:
            dirs.remove('.git')
        if '__pycache__' in dirs:
            dirs.remove('__pycache__')
            
        for file in files:
            if file.endswith('.py') or file.endswith('.md'):
                process_file(os.path.join(root, file))

if __name__ == "__main__":
    main()
