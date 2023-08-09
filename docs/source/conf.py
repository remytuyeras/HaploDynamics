# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'HaploDynamics'
copyright = '2023, Remy Tuyeras'
author = 'Remy Tuyeras'

release = '0.2-beta.4'
version = '0.2-beta.4'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']
