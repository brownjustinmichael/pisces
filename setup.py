from setuptools import setup
import os

setup(name='pisces_utils',
      version='0.1',
      description='A set of utilities for analyzing and running the hydrocode PISCES',
      url='https://github.com/brownjustinmichael/pisces',
      author='Justin Brown',
      author_email='jumbrown@ucsc.edu',
      license='MIT',
      packages=['pisces_utils'],
      install_requires=["sqlalchemy","celery","netCDF4","pyyaml","psycopg2"],
      zip_safe=False)