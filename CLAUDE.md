digraph CLAUDE {
    // Main process using proper shapes from style guide

    // WHEN: New request
    subgraph cluster_request {
        "Request arrives" [shape=ellipse];
        "Do I understand?" [shape=diamond];
        "Ask specific questions" [shape=box];
        "grep -r 'similar' ." [shape=plaintext];
        "Read similar files" [shape=box];
        "Will change >10 files?" [shape=diamond];
        "STOP: Get permission" [shape=octagon, style=filled, fillcolor=orange];
        "Begin implementation" [shape=doublecircle];

        "Request arrives" -> "Do I understand?";
        "Do I understand?" -> "Ask specific questions" [label="no"];
        "Do I understand?" -> "grep -r 'similar' ." [label="yes"];
        "Ask specific questions" -> "Do I understand?";
        "grep -r 'similar' ." -> "Read similar files";
        "Read similar files" -> "Will change >10 files?";
        "Will change >10 files?" -> "STOP: Get permission" [label="yes"];
        "Will change >10 files?" -> "Begin implementation" [label="no"];
        "STOP: Get permission" -> "Request arrives" [label="denied"];
        "STOP: Get permission" -> "Begin implementation" [label="approved"];
    }

    // WHEN: Writing code
    subgraph cluster_implement {
        "Write ONE failing test" [shape=box];
        "npm test" [shape=plaintext];
        "Test fails as expected" [shape=ellipse];
        "Write MINIMAL code" [shape=box];
        "npm test" [shape=plaintext];
        "Test passes?" [shape=diamond];
        "git add file1 file2" [shape=plaintext];
        "git commit -m 'desc'" [shape=plaintext];
        "More to implement?" [shape=diamond];
        "Implementation complete" [shape=doublecircle];

        "Begin implementation" -> "Write ONE failing test";
        "Write ONE failing test" -> "npm test";
        "npm test" -> "Test fails as expected";
        "Test fails as expected" -> "Write MINIMAL code";
        "Write MINIMAL code" -> "npm test";
        "npm test" -> "Test passes?";
        "Test passes?" -> "Write MINIMAL code" [label="no"];
        "Test passes?" -> "git add file1 file2" [label="yes"];
        "git add file1 file2" -> "git commit -m 'desc'";
        "git commit -m 'desc'" -> "More to implement?";
        "More to implement?" -> "Write ONE failing test" [label="yes"];
        "More to implement?" -> "Implementation complete" [label="no"];
    }

    // WHEN: Stuck
    subgraph cluster_stuck {
        "I'm stuck" [shape=ellipse];
        "Document the problem" [shape=box];
        "Third attempt?" [shape=diamond];
        "Say: I don't understand X" [shape=box];
        "Remove half the code" [shape=box];
        "Add debug output" [shape=box];
        "console.log(state)" [shape=plaintext];
        "Copy working example" [shape=box];

        "I'm stuck" -> "Document the problem";
        "Document the problem" -> "Third attempt?";
        "Third attempt?" -> "Say: I don't understand X" [label="yes"];
        "Third attempt?" -> "Remove half the code" [label="no"];
        "Say: I don't understand X" -> "Request arrives" [label="after help"];
        "Remove half the code" -> "Add debug output";
        "Add debug output" -> "console.log(state)";
        "console.log(state)" -> "Copy working example";
        "Copy working example" -> "I'm stuck" [label="still stuck"];
        "Copy working example" -> "Write MINIMAL code" [label="unstuck"];
    }

    // WHEN: Debugging
    subgraph cluster_debug {
        "Test is failing" [shape=ellipse];
        "Read error CAREFULLY" [shape=box];
        "Can reproduce?" [shape=diamond];
        "Simplify test case" [shape=box];
        "git diff HEAD~1" [shape=plaintext];
        "Find working example" [shape=box];
        "grep -r 'working_pattern' ." [shape=plaintext];
        "Compare differences" [shape=box];
        "Form ONE hypothesis" [shape=box];
        "Hypothesis correct?" [shape=diamond];
        "Apply minimal fix" [shape=box];

        "Test is failing" -> "Read error CAREFULLY";
        "Read error CAREFULLY" -> "Can reproduce?";
        "Can reproduce?" -> "Simplify test case" [label="no"];
        "Can reproduce?" -> "git diff HEAD~1" [label="yes"];
        "Simplify test case" -> "Can reproduce?";
        "git diff HEAD~1" -> "Find working example";
        "Find working example" -> "grep -r 'working_pattern' .";
        "grep -r 'working_pattern' ." -> "Compare differences";
        "Compare differences" -> "Form ONE hypothesis";
        "Form ONE hypothesis" -> "Hypothesis correct?";
        "Hypothesis correct?" -> "Apply minimal fix" [label="yes"];
        "Hypothesis correct?" -> "Form ONE hypothesis" [label="no"];
        "Apply minimal fix" -> "npm test";
    }

    // WHEN: Verifying complete
    subgraph cluster_verify {
        "Ready to complete" [shape=ellipse];
        "All tests pass?" [shape=diamond];
        "npm test" [shape=plaintext];
        "Style matches?" [shape=diamond];
        "npm run lint" [shape=plaintext];
        "Debug output removed?" [shape=diamond];
        "grep -r 'console.log' src/" [shape=plaintext];
        "TODOs updated?" [shape=diamond];
        "Fix issues" [shape=box];
        "Task complete" [shape=doublecircle];

        "Ready to complete" -> "All tests pass?";
        "All tests pass?" -> "npm test" [label="check"];
        "npm test" -> "Style matches?" [label="pass"];
        "npm test" -> "Fix issues" [label="fail"];
        "Style matches?" -> "npm run lint" [label="check"];
        "npm run lint" -> "Debug output removed?" [label="pass"];
        "npm run lint" -> "Fix issues" [label="fail"];
        "Debug output removed?" -> "grep -r 'console.log' src/" [label="check"];
        "grep -r 'console.log' src/" -> "TODOs updated?" [label="clean"];
        "grep -r 'console.log' src/" -> "Fix issues" [label="found"];
        "TODOs updated?" -> "Task complete" [label="yes"];
        "TODOs updated?" -> "Fix issues" [label="no"];
        "Fix issues" -> "Ready to complete";
    }

    // Critical warnings
    subgraph cluster_warnings {
        "NEVER git add -A" [shape=octagon, style=filled, fillcolor=red, fontcolor=white];
        "STOP if breaking rules" [shape=octagon, style=filled, fillcolor=red, fontcolor=white];
        "Test MUST come first" [shape=octagon, style=filled, fillcolor=red, fontcolor=white];
        "Say I don't understand" [shape=octagon, style=filled, fillcolor=orange];
    }

    // Process transitions (dotted = cross-process)
    "npm test" -> "I'm stuck" [label="confused", style=dotted];
    "npm test" -> "Test is failing" [label="fails", style=dotted];
    "Task complete" -> "Request arrives" [label="next", style=dotted];
    "Implementation complete" -> "Ready to complete" [style=dotted];
}

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