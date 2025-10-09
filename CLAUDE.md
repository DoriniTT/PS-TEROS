# Working with PS-TEROS and AiiDA-WorkGraph

This is a comprehensive guide for developing and working with the **PS-TEROS** code using **AiiDA** and **AiiDA-WorkGraph**.

---

## Quick Start Commands

* **Before you begin**

Make sure that you are in the 'psteros' profile in AiiDA by doing 'verdi profile set-default psteros', then do 'verdi status' to check if everything is working fine.
If the daemon is not running, you can start it by doing 'verdi daemon start'.

* **How the PS-TEROS is organized**

Folders:
* teros: Where all the code lives.
  * core: Where all the main python code for development are.
  * experimental: Here new features and experiments are made.
* examples: Exemples folder that test the full development PS-TEROS code (core files).
* legacy: Old files for backup
* docs: Documentation of the PS-TEROS code.

* **After every modification**

  ```bash
  verdi daemon restart
  ```

* **Python binary**

  ```bash
  source ~/envs/psteros/bin/activate && /home/thiagotd/envs/psteros/bin/python
  ```

* **Launch workflow**

  ```bash
  source ~/envs/psteros/bin/activate && /home/thiagotd/envs/psteros/bin/python /home/thiagotd/git/PS-TEROS/examples/vasp/update_psteros/psteros_vasp.py
  ```

* **Wait for results**: Use `sleep 15` (or more) before analyzing nodes

* **Analyze results**

  ```bash
  verdi process show <PK>
  verdi process report <PK>
  ```

* **Clear Python cache**

  ```bash
  find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
  ```

---

## Documentation

When developing a new feature to the code, please create the documentation in the folder where the test script is in the example, not in the root git directory.

---