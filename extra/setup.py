try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'PTR extraction pipeline',
    'author': 'Daniel Hermansson',
    'url': '',
    'download_url': '',
    'author_email': 'hedani@chalmers.se',
    'version': '0.1',
    'install_requires': ['jinja2'],
    'scripts': [''],
    'name': 'ptr_pipeline'
}

setup(**config)