from setuptools import setup, find_packages

setup(name='modelseedpy_ext',
      version='0.0.9',
      description='ModelSEEDpy External Tools',
      url='https://github.com/Fxe/ModelSEEDpy-Ext',
      author='Filipe Liu',
      author_email='fliu@anl.gov',
      license='MIT',
      packages=find_packages(),
      package_data={
          "modelseedpy_ext": ["profiler/*"],
      },
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Intended Audience :: Science/Research",
          "Operating System :: OS Independent",
          "Programming Language :: Python :: 3.8",
          "Programming Language :: Python :: 3.9",
          "Programming Language :: Python :: 3.10",
          "Programming Language :: Python :: 3.11",
          "Natural Language :: English",
      ],
      install_requires=[
          "modelseedpy >= 0.3.0",
          "requests",
          "pandas >= 1.0.0",
          "networkx >= 2.4",
          "cobra >= 0.17.1",
          "pymongo",
          "numpy",
          "pyarrow",
          "lxml",
          "pyArango"
      ],
      zip_safe=True)
