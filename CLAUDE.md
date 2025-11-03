# Working with PS-TEROS and AiiDA-WorkGraph

This is a comprehensive guide for developing and working with the **PS-TEROS** code using **AiiDA** and **AiiDA-WorkGraph**.

---

## Quick Start Commands

* **Before you begin**

Make sure that you are in the 'presto' profile in AiiDA by doing 'verdi profile set-default presto', then do 'verdi status' to check if everything is working fine.
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

  **When running tests**

  Use the VASP6.5.0@cluster02 code. It uses a cluster with 24 processors which have these information:
-----------------------  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PK                       122432
UUID                     699ca05a-6901-429f-a5ee-6fe21edd9db7
Type                     core.code.installed
Pk                       122432
Uuid                     699ca05a-6901-429f-a5ee-6fe21edd9db7
Node type                data.core.code.installed.InstalledCode.
Process type
Repository metadata      {}
Ctime                    2025-11-03 07:48:43.086994-03:00
Mtime                    2025-11-03 07:48:43.120913-03:00
Label                    VASP6.5.0
Description              VASP 6.5.0
Attributes               {'with_mpi': None, 'append_text': '', 'input_plugin': 'vasp.vasp', 'prepend_text': '. /home/umsl02/Programs/spack/share/spack/setup-env.sh', 'use_double_quotes': False, 'filepath_executable': '/usr/local/bin/vasp_std', 'wrap_cmdline_params': False}
Extras                   {'hidden': False, '_aiida_hash': 'a0bfbc12e319d46a026a97fd3f0e9216ac7b6bc6febb106eb6c8f273cb93c576'}
Computer                 cluster02 (umsl02), pk: 95
User                     aiida@localhost
Source
Default calc job plugin  vasp.vasp
Use double quotes        False
With mpi
Prepend text             . /home/umsl02/Programs/spack/share/spack/setup-env.sh
Append text
Filepath executable      /usr/local/bin/vasp_std
-----------------------  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  - It do not have a queue.
  - Use the num_cores_per_machine=24, like this:
  options = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 24, 
    },
}

* **Python binary**

  ```bash
  source ~/envs/aiida/bin/activate && /home/thiagotd/envs/aiida/bin/python
  ```

* **Launch workflow**

  ```bash
  source ~/envs/aiida/bin/activate && /home/thiagotd/envs/aiida/bin/python /home/thiagotd/git/PS-TEROS/examples/vasp/update_psteros/psteros_vasp.py
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

## GitHub Workflow

After you successfully implement and test a new feature, you may add, commit, and merge it into the **`develop`** branch. Below is the recommended workflow.

---

### Repository Structure and Branch Management

The main development branch is located at:

```
/home/thiagotd/git/PS-TEROS
```

For each new feature or update to the PS-TEROS code, create a dedicated folder in the /home/thiagotd/git/PS-TEROS/teros/experimental path. This ensures that each feature is isolated and can be tested independently before merging.

You will perform all feature development in this directory.

---

### Testing Before Merging

Every feature must be **fully tested** before merging into `develop`.
To do this, create a **new folder** inside:

```
/home/thiagotd/git/PS-TEROS/examples/
```

Use this folder to run complete test cases that validate your new functionality. You may also create small scripts to test specific parts of the code, but final validation must be performed using **production-like parameters**, similar to those used in the existing example folders.

Your implementation is considered correct and ready to merge when:

* The main node returns a `[0]`, indicating successful completion.
* The feature you implemented appears correctly in the main node output and behaves as expected.
* The calculations may take some time, so use 'sleep' bash command to wait until you clearly see the modification that we are implemeted work correctly when checking the calculation with 'verdi process show <PK>' or 'verdi process report <PK>'.

---

### Committing and Pushing Changes

Once the above criteria are met, you can commit and push your changes:

```bash
git add -A
git commit -m "Describe the feature or fix clearly."
git pull
git push
```

Then return to the `develop` branch:

```bash
cd /home/thiagotd/git/PS-TEROS
```
---

## Documentation

When developing a new feature to the code, please create the documentation in the folder where the test script is in the example, not in the root git directory.

---
