from distutils.core import setup
  
setup(
    name='psi4fockci',
    version='0.1.0',
    author='Shannon E. Houck',
    author_email='shouck@vt.edu',
    packages=['psi4fockci', 'psi4fockci.test'],
    scripts=[],
    url='',
    description='Psi4 RAS-SF-IP/EA code.',
    long_description=open('README.md').read(),
)
