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

  **When running tests**

  Use the VASP-6.4.1@cluster02 code. It uses a cluster with 24 processors which have these information:
PK                       42941
UUID                     6bb827ad-8615-430c-8231-9267c68d67ff
Type                     core.code.installed
Pk                       42941
Uuid                     6bb827ad-8615-430c-8231-9267c68d67ff
Node type                data.core.code.installed.InstalledCode.
Label                    VASP-6.4.1
Description              VASP 6.4.1
Attributes               {'with_mpi': None, 'append_text': '', 'input_plugin': 'vasp.vasp', 'prepend_text': 'module load vasp\n', 'use_double_quotes': False, 'filepath_executable': '/usr/sw/vasp.6.4.1/bin/vasp_std', 'wrap_cmdline_params': False}
Computer                 cluster06 (umsl06), pk: 6
Default calc job plugin  vasp.vasp
Prepend text             module load vasp
Append text
Filepath executable      /usr/sw/vasp.6.4.1/bin/vasp_std

  - It do not have a queue.
  - Use the num_cores_per_machine=24, like this:
  options = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 24, 
    },
    # Optional settings:
    # 'prepend_text': 'module load vasp/6.4.2',  # Commands before VASP
    # 'custom_scheduler_commands': '#SBATCH --constraint=haswell',
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

For each new feature or update to the PS-TEROS code, create a dedicated worktree branch. This ensures that each feature is isolated and can be tested independently before merging.

```bash
git worktree add ../worktree/PS-TEROS/<feature_name> <feature_branch>
```

Then move into the new branch directory:

```bash
cd ../worktree/PS-TEROS/<feature_name>
```

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

And merge your feature branch:

```bash
git merge <feature_branch>
```

---

### Handling Merge Conflicts

If conflicts occur, **stop immediately** so they can be reviewed.
You may resolve simple conflicts yourself (e.g., cache files, `__init__.py`, etc.), but for more complex cases, wait for a manual review.

Do **not delete** the `<feature_branch>` after merging.
A final validation will be performed on the `develop` branch to ensure that the integration works as intended.

---

## Documentation

When developing a new feature to the code, please create the documentation in the folder where the test script is in the example, not in the root git directory.

---