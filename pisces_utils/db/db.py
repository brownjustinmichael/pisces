import sys
import os
import hashlib
import getpass
import json
import yaml
import subprocess

import numpy as np
import numpy.ma as ma
import netCDF4

import sqlalchemy
from sqlalchemy.sql.expression import func
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm.exc import NoResultFound

import pisces_utils.config as config

# The numpy types are not recognized by psycopg2, so we need to register them
from psycopg2.extensions import register_adapter, AsIs
def adapt_numpy_float64(numpy_float64):
  return AsIs(numpy_float64)
register_adapter(np.float64, adapt_numpy_float64)

def adapt_numpy_int32(numpy_int32):
  return AsIs(numpy_int32)
register_adapter(np.int32, adapt_numpy_int32)

# Create the declarative base from which the entries will derive
Base=declarative_base()

# See if a password is defined; if so, use it
# Otherwise, request the password at the terminal
# This doesn't seem to work inside SublimeText, so the password can be specified manually by setting password= here
try:
    engine = sqlalchemy.create_engine ('postgresql://justinbrown:%s@loki.ucsc.edu:5432/piscesdb' % password, echo = False)
except NameError:
    password = getpass.getpass ()
    engine = sqlalchemy.create_engine ('postgresql://justinbrown:%s@loki.ucsc.edu:5432/piscesdb' % password, echo = False)

# Reflect the metadata from the existing database
Base.metadata.reflect(bind=engine)
Session=sqlalchemy.orm.sessionmaker(bind=engine)

# This dictionary allows easy translation of types into types that sql will recognize
strtypes={int: "integer", 
          float: "double precision", 
          str: "text", 
          np.float64: "double precision", 
          np.int32: "integer", 
          bool: "bool", 
          type(None): "text"}
