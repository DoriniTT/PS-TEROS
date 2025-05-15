# How to create a documentation

This directory contains the source files for the documentation of TEROS.

## Building the documentation

To build the documentation locally:

1. Install the documentation requirements:
   ```
   pip install -e .[docs]
   ```

2. Change to the docs directory:
   ```
   cd docs
   ```

3. Build the HTML documentation:
   ```
   make html
   ```

4. View the documentation:
   ```
   open build/html/index.html
   ```

## Documentation structure

- `source/index.rst`: Main entry point
- `source/installation.rst`: Installation instructions
- `source/tutorial.rst`: Tutorial for new users
- `source/theory.rst`: Theoretical background
- `source/api.rst`: API reference documentation
- `source/examples.rst`: Example usage scripts
- `source/contributing.rst`: Guidelines for contributors
- `source/authors.rst`: List of authors and contributors
- `source/history.rst`: Release history and changelog

## Adding new documentation

1. Create a new `.rst` file in the `source` directory
2. Add the file to the table of contents in `index.rst`
3. Build the documentation to test your changes

## Publishing

Documentation is automatically built and published on [Read the Docs](https://readthedocs.org/) 
when changes are pushed to the main branch of the repository.
