from setuptools import setup, find_packages

setup(name='modelseedpy_ext',
      version='0.0.9',
      description='ModelSEEDpy External Tools',
      url='https://github.com/Fxe/ModelSEEDpy-Ext',
      author='Filipe Liu',
      author_email='fliu@anl.gov',
      license='MIT',
      packages=['modelseedpy_ext'],
      install_requires=[
          # "modelseedpy >= 1.0.0", # when available in pypi
          "pandas >= 1.0.0",
          "networkx >= 2.4",
          "cobra >= 0.17.1"
      ],
      zip_safe=True)
