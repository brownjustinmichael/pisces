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

password="occuthe9"

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

    def __init__(self, entry=None, default_file="../src/defaults.yaml", *args, **kwargs):
        self._entry=entry
        if self._entry is None:
            self._entry=self.Table(*args, **kwargs)
            self.add_params(default_file)

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
    def _process(cls, value, dictionary=None, current_key=""):
        """
        Process a YAML nested dictionary into a flat dictionary with different levels becoming underscore-separated keys
        """
        if dictionary is None:
            dictionary={}
        if isinstance(value, dict):
            for key in value:
                if current_key:
                    composite_key="__".join([current_key, key])
                else:
                    composite_key=key
                cls._process(value[key], dictionary, composite_key)
        else:
            dictionary[current_key]=value
        return dictionary

    @classmethod
    def add_columns(cls, session=None, **kwargs):
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

    # def __getattr__(self, attr):
    #     try:
    #         return getattr(self.entry, attr)
    #     except AttributeError:
    #         pass

    #     if not isinstance(self.entry, self.Table):
    #         tmp=self.entry
    #         self.entry=self.Table()
    #         for column in tmp.__table__.columns:
    #             try:
    #                 setattr(self.entry, column.key, getattr(tmp, column.key))
    #             except AttributeError:
    #                 pass
    #     return getattr(self.entry, attr)
      
    def steps(self, session = None):
        if session is None:
            session=Session.object_session(self.entry)
        return session.query(StepEntry.Table).filter(StepEntry.Table.simulation==self.entry)

    def add_params(self, param_file):
        if param_file:
            try:
                translation=self._process(yaml.load(open(param_file)))
            except OSError:
                translation=self._process(yaml.load(param_file))
            self.add_columns(**translation)
            for param in translation:
                setattr(self.entry, param, translation[param])

    @classmethod
    def from_params(cls, *args, default_file="../src/defaults.yaml", session=None, **kwargs):
        """
        Generate a simulation entry from a YAML string of parameters. If the simulation already exists in the database, return that instead
        """
        # Parse the input parameters
        translation={}
        for params in (default_file,) + args:
            if not params:
                continue
            try:
                params=open(params)
            except IOError:
                pass
            translation=cls._process(yaml.load(params), translation)
            for kw in kwargs:
                translation[kw]=kwargs[kw]

            cls.add_columns(session=session, **translation)

        # No acceptable simulation has been found, so we generate one from scratch
        entry=cls(default_file="")
        for param in translation:
            setattr(entry.entry, param, translation[param])
        return entry

    @classmethod
    def _unprocess(cls, value, dictionary, keys):
        if len(keys) == 1:
            if value is not None:
                dictionary[keys[0]]=value
            else:
                return
        else:
            if keys[0] not in dictionary:
                dictionary[keys[0]]={}
            cls._unprocess(value, dictionary[keys[0]], keys[1:])

    def to_file(self, file_name="config.yaml", **kwargs):
        result={}
        for column in self.Table.__table__.columns:
            if column.key != "id":
                if column.key not in kwargs:
                    self._unprocess(getattr(self.entry, column.key), result, column.key.split("__"))
                else:
                    self._unprocess(getattr(kwargs, column.key), result, column.key.split("__"))
        f=open(file_name, "w")
        yaml.dump(result, f, default_flow_style=False)

    def resume(self, execute="../../run/pisces", cwd=None, reset=False, debug=None, session=None, **kwargs):
        if session is None:
            session=Session.object_session(self.entry)
        # final=self.steps(session).order_by(StepEntry.Table.step.desc()).first()

        # if final is not None:
        #     kwargs["output__number"]=self.entry.output__number + 1
        #     kwargs["input__file"]=self.entry.dump__file
        #     kwargs["input__directory"]=self.entry.dump__directory
        # if final is None or (final.t < self.entry.time__stop and final.step < self.entry.time__steps):
        currentdir=os.getcwd()
        if cwd:
            os.chdir(cwd)
        else:
            os.chdir(self.entry.cwd)
        self.to_file(**kwargs)

        try:
            np = str(self.entry.np if self.entry.np else 1), init
        except AttributeError:
            np = "1"

        exit=subprocess.call(["mpiexec", "-np", np, execute] + ([] if debug is None else ["-D%i" % debug]))
        if exit != 0:
            raise RuntimeError("The subprocess exited with a nonzero exit code of %i" % exit)
        new=[]

        for i in range(int (np)):
            new+=StepEntry.from_file(session, os.path.join(self.entry.root, self.entry.output__directory, ((self.entry.output__stat__file % i) % self.entry.output__number) + ".cdf"), sim=self, reset=reset)
        os.chdir(currentdir)
        return new

    def run(self, init="../../run/isces_init", cwd=None, debug=None, session=None, **kwargs):
        if session is None:
            session=Session.object_session(self.entry)

        # try:
        #     self.entry=session.query(SimulationEntry).filter(SimulationEntry.hash==self.calculate_hash()).one()
        # except NoResultFound:
        print ("DEBUG IS ", debug, [] if debug is None else ["-D%i" % debug])
        currentdir=os.getcwd()
        if cwd:
            os.chdir(cwd)
        else:
            os.chdir(self.entry.cwd)
        self.to_file()
        os.makedirs(os.path.join(self.entry.root, self.entry.input__directory), exist_ok=True)
        os.makedirs(os.path.join(self.entry.root, self.entry.output__directory), exist_ok=True)
        try:
            np = str(self.entry.np if self.entry.np else 1)
        except AttributeError:
            np = "1"
        subprocess.call(["mpiexec", "-np", np, init] + ([] if debug is None else ["-D%i" % debug]))
        os.chdir(currentdir)

        return self.resume(cwd=cwd, debug=debug, session=session, **kwargs)

    def same_sub(self, query, key, tolerance=1.e-4):
        for column in SimulationEntry.Table.__table__.columns:
            if column.key.startswith(key):
                value=getattr(self.entry, column.key)
                print (column.key, value)
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
