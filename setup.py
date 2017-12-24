#!/usr/bin/env python
install_requires = ['numpy','xarray']
tests_require=['nose','coveralls']
# %%
from setuptools import setup,find_packages

setup(name='AeroPlanets',
      packages=find_packages(),
      author='Michael Hirsch, Ph.D.',
      description='1-D ionosphere particle penetration simulatio',
      url='https://github.com/scivision/aeroplanets',
      version='0.2.0',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'License :: OSI Approved :: GPL License',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3',
      ],
      install_requires=install_requires,
      python_requires='>=3.6',
      tests_require=tests_require,
      extras_require={'tests':tests_require,
                      'plot':['matplotlib','seaborn',]},
	  )

