# -*- coding: utf-8 -*-

import datetime
import importlib
import inspect
import os
import sys

pardir = os.path.abspath(os.pardir)
sys.path.insert(1, pardir)
# so existing cache directory is used when src is imported
os.environ['WORKDIR'] = pardir

project = 'hic-param-est'
version = release = ''
author = 'Jonah Bernhard'
copyright = '{} {}'.format(datetime.date.today().year, author)

source_suffix = '.rst'
master_doc = 'index'

templates_path = ['_templates']
html_static_path = ['_static']
exclude_patterns = ['_build']

default_role = 'py:obj'
pygments_style = 'sphinx'

html_theme = 'sphinx_rtd_theme'
html_title = ''
html_domain_indices = False
html_use_index = False

html_context = dict(show_source=False)

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.extlinks',
    'sphinx.ext.linkcode',
]

autodoc_member_order = 'bysource'

github_url = 'https://github.com/jbernhard/hic-param-est/blob/master/'

extlinks = {
    'ghlink': (github_url + '%s', '')
}

def linkcode_resolve(domain, info):
    """
    Get the github URL (with line numbers) for a code object.

    """
    if domain != 'py':
        return

    try:
        obj = importlib.import_module(info['module'])
        for name in info['fullname'].split('.'):
            obj = getattr(obj, name)
        lines, startline = inspect.getsourcelines(obj)
        sourcefile = os.path.basename(inspect.getsourcefile(obj))
    except Exception:
        return

    return '{}src/{}#L{}-L{}'.format(
        github_url, sourcefile,
        startline, startline + len(lines) - 1
    )
