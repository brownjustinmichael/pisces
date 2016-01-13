import sys
import os
import hashlib
import getpass
import json
import yaml
import subprocess

import numpy as np
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

class SimulationEntry(object):
    if "simulations" not in Base.metadata.tables:
        # If not, create one with columns id and hash
        Table=type("SimulationTable", (Base,), {"__table__": sqlalchemy.Table("simulations", Base.metadata, sqlalchemy.Column("id", sqlalchemy.Integer, primary_key=True))})
        Base.metadata.tables["simulations"].create(engine)
    else:
        Table=type("SimulationTable", (Base,), {"__table__": sqlalchemy.Table("simulations", Base.metadata, autoload=True, autoload_with=engine, extend_existing=True)})

    def __init__(self, entry=None, *args, **kwargs):
        self._entry=entry
        if self._entry is None:
            self._entry=self.Table(*args, **kwargs)
            self.add_params(**kwargs)

    @property
    def entry(self):
        """
        Return the entry associated with the current simulation. This property first checks that the database entry has the correct table to be in the database and if not, copies the contents of the previous entry into the correct table.
        """
        if not isinstance(self._entry, self.Table):
            tmp=self._entry
            self._entry=self.Table()
            for column in tmp.__table__.columns:
                try:
                    setattr(self._entry, column.key, getattr(tmp, column.key))
                except AttributeError:
                    pass
        return self._entry

    @classmethod
    def add_columns(cls, session=None, **kwargs):
        """
        Add columns to the underlying database. kwargs should be a dictionary where the keys are the parameter names and the values are default values for the column.
        """
        changed=False
        for arg in kwargs:
            try:
                getattr(cls.Table, arg)
            except AttributeError:
                changed=True
                if session is not None and session.dirty:
                    raise RuntimeError("Dirty session")
                conn=engine.connect()
                print("TABLE CHANGE: ADDING COLUMN %s TO SIMULATIONS" % arg)
                conn.execute("alter table simulations add column %s %s null" % (arg, strtypes[type(kwargs[arg])]))
                conn.close()
        if changed:
            Base.metadata.reflect(bind=engine)
            cls.Table=type("SimulationTable", (Base,), {"__table__": sqlalchemy.Table("simulations", Base.metadata, autoload=True, autoload_with=engine, extend_existing=True)})

    def steps(self, session = None):
        if session is None:
            session=Session.object_session(self.entry)
        return session.query(StepEntry.Table).filter(StepEntry.Table.simulation==self.entry)

    def add_params(self, **kwargs):
        self.add_columns(**kwargs)
        for param in kwargs:
            setattr(self.entry, param, kwargs[param])

    @classmethod
    def from_config(cls, session=None, **kwargs):
        """
        Generate a simulation entry from a YAML string of parameters. If the simulation already exists in the database, return that instead
        """
        # Parse the input parameters
        translation=config.process(kwargs)
        cls.add_columns(session=session, **translation)

        # No acceptable simulation has been found, so we generate one from scratch
        return cls(kwargs)

    def same_sub(self, query, key, tolerance=1.e-4):
        for column in SimulationEntry.Table.__table__.columns:
            if column.key.startswith(key):
                value=getattr(self.entry, column.key)
                if type(value) != float:
                    query=query.filter(getattr(SimulationEntry.Table, column.key)==value)
                else:
                    query=query.filter(func.abs(getattr(SimulationEntry.Table, column.key)-value)<tolerance)
        return query

    def clone(self):
        sim=SimulationEntry(default_file="")
        for column in SimulationEntry.Table.__table__.columns:
            setattr(sim.entry, column.key, getattr(self.entry, column.key))
        return sim

class StepEntry(object):
    """
    This class produces a class named StepEntry that represents the steps table in the database
    When MetaStepEntry() is called, it automatically reflects off the existing database and replaces the StepEntry class of this module with the newly-constructed class
    """
    if "steps" not in Base.metadata.tables:
        # If not, create one with columns id and hash
        Table=type("StepTable", (Base,), {"__table__": 
            sqlalchemy.Table("steps", Base.metadata, sqlalchemy.Column("id", sqlalchemy.Integer, primary_key=True), 
                sqlalchemy.Column("file", sqlalchemy.String), 
                sqlalchemy.Column("line", sqlalchemy.Integer), 
                sqlalchemy.Column("step", sqlalchemy.Integer), 
                sqlalchemy.Column("simulation_id", sqlalchemy.Integer, sqlalchemy.ForeignKey("simulations.id")), 
                sqlalchemy.UniqueConstraint("file", "line", name="unique_file_line")), 
            "simulation": sqlalchemy.orm.relationship(SimulationEntry.Table)})
        Base.metadata.tables["steps"].create(engine)
    else:
        Table=type("StepTable", (Base,), {"__table__": sqlalchemy.Table("steps", Base.metadata, autoload=True, autoload_with=engine, extend_existing=True), "simulation": sqlalchemy.orm.relationship(SimulationEntry.Table)})

    def __init__(self, entry=None, *args, **kwargs):
        if entry:
            self._entry=entry
        else:
            self._entry=self.Table(*args, **kwargs)

    @property
    def entry(self):
        if not isinstance(self._entry, self.Table):
            tmp=self._entry
            self._entry=self.Table()
            for column in tmp.__table__.columns:
                try:
                    setattr(self._entry, column.key, getattr(tmp, column.key))
                except AttributeError:
                    pass
        return self._entry

    @classmethod
    def from_file (cls, session, file_name, sim=None, reset=False):
        """
        Generate a new set of StepEntries from a cdf file.
        If reset is set, any existing entries with that file will be deleted.
        """
        # Gather the data from file
        data=netCDF4.Dataset(file_name)

        # Iterate over the data, checking that all scalars are present as columns in the database, adding those not present
        times=0
        changed=False
        for variable in data.variables:
            if len(data[variable].shape) == 1:
                times=data[variable].shape[0]
                try:
                    getattr(cls.Table, variable)
                except AttributeError:
                    changed=True
                    if session.dirty:
                        raise RuntimeError("Dirty session")
                    conn=engine.connect()
                    print("TABLE CHANGE: ADDING COLUMN %s TO STEPS" % variable)
                    conn.execute("alter table steps add column %s %s null" % (variable, strtypes[type(data[variable][0])]))
                    conn.close()
        if changed:
            cls.Table=type("StepTable", (Base,), {"__table__": sqlalchemy.Table("steps", Base.metadata, autoload=True, autoload_with=engine, extend_existing=True), "simulation": sqlalchemy.orm.relationship(SimulationEntry.Table)})

        # Generate a new simulation from the parameters read in from that file
        if sim is None:
            sim=SimulationEntry.from_params(session=session, **{attr: data.getncattr(attr) for attr in data.ncattrs()})
        else:
            for attr in data.ncattrs():
                if attr != "params":
                    setattr(sim.entry, attr, data.getncattr(attr))
                else:
                    sim.add_params(data.getncattr("params"))

        # # If reset is set, delete all entries in the database matching the current file
        if reset:
            for entry in session.query(StepEntry.Table).filter(StepEntry.Table.file==os.path.abspath(file_name)).all():
                session.delete(entry)

        # Irrelevant, but the code breaks without this line
        session.flush ()
        # session.query(StepEntry).filter(StepEntry.file==os.path.abspath(file_name)).count()

        # Add the lines from the code as new entries
        entries=[]
        for time in range(times):
            entries.append(StepEntry.Table())
            entries[-1].file=os.path.abspath(file_name)
            entries[-1].line=time
            entries[-1].simulation_id=sim.entry.id
            for variable in data.variables:
                if len(data[variable].dimensions) == 1:
                    setattr(entries[-1], variable, data[variable][time])
            session.add(entries[-1])

        # session.flush()

        return entries

# sqlalchemy.orm.mapper(SimulationEntry, SimulationEntry.Table.__table__)
# sqlalchemy.orm.mapper(StepEntry, StepEntry.Table.__table__)
