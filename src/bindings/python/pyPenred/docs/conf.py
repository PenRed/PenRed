# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pyPenred'
copyright = '2025, penred colaboration'
author = 'penred colaboration'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'sphinx_copybutton',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"

# Logo configuration
html_logo = "_static/logo.png"
html_title = "pyPenred API"

html_static_path = ['_static']

autodoc_default_options = {
    'imported-members': True,
    'special-members': '__init__',
    'show-inheritance': True,
    'undoc-members': True,
}

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}

# Better type hint formatting
autodoc_typehints = "description"
autodoc_typehints_format = "short"
python_use_unqualified_type_names = True
typehints_use_signature = True  # More compact display
typehints_use_signature_return = True  # For return type
add_function_parentheses = True  # Ensures parentheses are shown

# Add these Napoleon settings (below your existing napoleon config):
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_ivar = True

# Configure copybutton for code blocks
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_remove_prompts = True

html_css_files = ['custom.css']
html_js_files = ['custom.js']
