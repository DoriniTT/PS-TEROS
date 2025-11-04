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

  Use the VASP6.5.1@cluster02 code. It uses a cluster with 24 processors which have these information:

PK                       124130
UUID                     2cc402ca-891e-4844-b84b-bacec1df2c0e
Type                     core.code.installed
Pk                       124130
Uuid                     2cc402ca-891e-4844-b84b-bacec1df2c0e
Node type                data.core.code.installed.InstalledCode.
Process type
Repository metadata      {}
Ctime                    2025-11-03 20:51:03.443188-03:00
Mtime                    2025-11-03 20:51:03.482991-03:00
Label                    VASP-6.5.1
Description              VASP 6.5.1
Attributes               {'with_mpi': None, 'append_text': '', 'input_plugin': 'vasp.vasp', 'prepend_text': 'source /home/umsl02/Programs/load-vasp-6.5.1.sh\n', 'use_double_quotes': False, 'filepath_executable': '/home/umsl02/Programs/spack/opt/spack/linux-skylake/vasp-6.5.1-mxw7wmvghgvg3vkknylaqrvminkn6gwj/bin/vasp_std', 'wrap_cmdline_params': False}
Extras                   {'hidden': False, '_aiida_hash': 'de588cba661946ef235a7dcc48dec3e415862c107fd8318074f63dc7152f869c'}
Computer                 cluster02 (umsl02), pk: 97
User                     aiida@localhost
Source
Default calc job plugin  vasp.vasp
Use double quotes        False
With mpi
Prepend text             source /home/umsl02/Programs/load-vasp-6.5.1.sh
Append text
Filepath executable      /home/umsl02/Programs/spack/opt/spack/linux-skylake/vasp-6.5.1-mxw7wmvghgvg3vkknylaqrvminkn6gwj/bin/vasp_std

It do not have a queue.
Use the num_cores_per_machine=24, like this:
options = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 24, 
    },
}
Keep in mind that it can only launch one calculation at a time, so you MUST set max_concurrent_jobs=1 in the workgraph.
When correcting errors, dont forget to kill the main workgraph (and all its sub nodes) node with verdi process kill <WORKGRAPH_PK>. In addition, it will not kill the calculation in the cluster02, so you need to manually run ssh cluster02 killall vasp_std. After this you may run the test script again.

* **Python binary**

  ```bash
  source ~/envs/aiida/bin/activate && /home/thiagotd/envs/aiida/bin/python
  ```

* **Launch workflow**

  ```bash
  source ~/envs/aiida/bin/activate && /home/thiagotd/envs/aiida/bin/python /home/thiagotd/git/PS-TEROS/examples/vasp/step_x_example.py
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

## Documentation

When developing a new feature to the code, please create the documentation in the folder where the test script is in the example, not in the root git directory.

---
