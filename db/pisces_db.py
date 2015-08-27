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

# This dictionary allows easy translation of types into types that sql will recognize
strtypes={int: "integer", 
          float: "real", 
          str: "text", 
          np.float64: "real", 
          np.int32: "integer", 
          bool: "bool", 
          type(None): "text"}

class BaseSimulationEntry(object):
    @classmethod
    def Factory(cls):
        """
        This function produces a class named SimulationEntry that represents the simulations table in the database
        When SimulationEntryFactory() is called, it automatically reflects off the existing database and replaces the SimulationEntry class of this module with the newly-constructed class
        """
        @property
        def steps(self):
            """
            Return the steps associated with this simulation
            """
            return Session.object_session(self).query(StepEntry).filter(StepEntry.simulation==self)

        # Extend the current table using the one in the database
        tab=sqlalchemy.Table("simulations", Base.metadata, autoload=True, autoload_with=engine, extend_existing=True)

        # Construct a new SimulationEntry class
        newcls=type("SimulationEntry", (BaseSimulationEntry, Base), {"__table__": tab, "steps": steps})

        # Set the module level SimulationEntry to the new class
        setattr(sys.modules[__name__], "SimulationEntry", newcls)
        
        # Reset the StepEntry metaclass to correctly associate the new class with the foreign key
        BaseStepEntry.Factory()
        return newcls

    def __new__(cls, default_file="../src/defaults.yaml", *args, **kwargs):
        if default_file:
            return cls.from_params(default_file)
        else:
            return object.__new__(cls, *args, **kwargs)

    def __init__(self, default_file="../src/defaults.yaml", *args, **kwargs):
        Base.__init__(self, *args, **kwargs)

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

    def to_file(self, file_name="config.yaml"):
        result={}
        for column in self.__table__.columns:
            if column.key != "id" and column.key != "hash":
                self._unprocess(getattr(self, column.key), result, column.key.split("__"))
        f=open(file_name, "w")
        yaml.dump(result, f, default_flow_style=False)

    @classmethod
    def from_params(cls, params, session=None, **kwargs):
        """
        Generate a simulation entry from a YAML string of parameters. If the simulation already exists in the database, return that instead
        """
        # Parse the input parameters
        try:
            params=open(params)
        except IOError:
            pass
        translation=cls._process(yaml.load(params))
        for kw in kwargs:
            translation[kw]=kwargs[kw]

        # Check that each parameter is a column in the database; if not, add it
        changed=False
        for param in translation:
            try:
                getattr(cls, param)
            except AttributeError:
                changed=True
                conn=engine.connect()
                print(param,translation[param])
                conn.execute("alter table simulations add column %s %s null" % (param, strtypes[type(translation[param])]))
                conn.close()
        if changed:
            cls=cls.Factory()

        # Generate a hash of the parameters and look up this hash in the database
        hash_value=cls._hash_dictionary(translation)
        if session:
            try:
                entry=session.query(SimulationEntry).filter(SimulationEntry.hash==hash_value).one()
                return entry
            except NoResultFound:
                pass

        # No acceptable simulation has been found, so we generate one from scratch
        entry=cls(default_file="")
        for param in translation:
            setattr(entry, param, translation[param])
        return entry

    @staticmethod
    def _hash_dictionary(dictionary):
        copy=dictionary.copy()
        for key in dictionary:
            if key.startswith("output") or key.startswith("input") or key.startswith("time") or key == "time__steps" or key == "time__stop":
                copy.pop(key)
        return hashlib.md5(json.dumps(copy, sort_keys=True).encode("utf-8")).hexdigest()

    def calculate_hash(self):
        dictionary = {}
        for column in self.__table__.columns:
            if column.key != "id" and column.key != "hash" and not column.key.startswith("output") and not column.key == "time__steps" and not column.key == "time__stop" and not column.key.startswith("input"):
                dictionary[column.key]=getattr(self, column.key)
        print (dictionary)
        return self._hash_dictionary(dictionary)

    def __eq__(self, other):
        return self.calculate_hash() == other.calculate_hash()

    def add_params(self, param_file):
        translation=self._process(yaml.load(open(param_file)))
        print(translation)
        for param in translation:
            try:
                getattr(self, param)
                print ("PARAM ", param, " EXISTS")
            except AttributeError:
                conn=engine.connect()
                print(param,translation[param])
                conn.execute("alter table simulations add column %s %s null" % (param, strtypes[type(translation[param])]))
                conn.close()
        result=self.__class__()
        for column in self.__table__.columns:
            setattr(result, column.key, getattr(self, column.key))
        for param in translation:
            setattr(result, param, translation[param])
        return result

    def resume(self, execute="../../run/pisces", np=1, cwd=None, reset=False, debug=None, **kwargs):
        final=self.steps.order_by(StepEntry.step.desc()).first()
        print("Final is ", final)
        print(self.id)
        if final is not None:
            self.output__number+=1
            self.input__file=self.dump__file
            self.input__directory=self.dump__directory
        if final is None or (final.t < self.time__stop and final.step < self.time__steps):
            currentdir=os.getcwd()
            if cwd:
                os.chdir(cwd)
            else:
                os.chdir(self.cwd)
            self.to_file()
            subprocess.call(["mpiexec", "-np", str(np), execute] + ([] if debug is None else ["-D%i" % debug]))
            session=Session.object_session(self)
            for i in range(self.np):
                print ("Adding")
                new=StepEntry.from_file(session, os.path.join(self.root, self.output__directory, ((self.output__stat__file % i) % self.output__number) + ".cdf"), sim=self, reset=reset)
            os.chdir(currentdir)

    def run(self, init="../../run/isces_init", np=1, cwd=None, debug=None, **kwargs):
        session=Session.object_session(self)
        this=self
        try:
            print(self.calculate_hash())
            this=session.query(SimulationEntry).filter(SimulationEntry.hash==self.calculate_hash()).one()
            print("FOUND PREVIOUS")
        except NoResultFound:
            this.np=np
            currentdir=os.getcwd()
            if cwd:
                os.chdir(cwd)
            else:
                os.chdir(this.cwd)
            this.to_file()
            os.makedirs(os.path.join(this.root, this.input__directory), exist_ok=True)
            os.makedirs(os.path.join(this.root, this.output__directory), exist_ok=True)
            subprocess.call(["mpiexec", "-np", str(np), init] + ([] if debug is None else ["-D%i" % debug]))
            os.chdir(currentdir)

        print("RUN HERE")
        this.resume(np=np, cwd=cwd, debug=debug, **kwargs)

class BaseStepEntry(object):
    """
    This class produces a class named StepEntry that represents the steps table in the database
    When MetaStepEntry() is called, it automatically reflects off the existing database and replaces the StepEntry class of this module with the newly-constructed class
    """
    @classmethod
    def Factory(cls):
        """
        This class produces a class named StepEntry that represents the steps table in the database
        When MetaStepEntry() is called, it automatically reflects off the existing database and replaces the StepEntry class of this module with the newly-constructed class
        """
        # Extend the current table by reflecting off the database
        tab=sqlalchemy.Table("steps", Base.metadata, autoload=True, autoload_with=engine, extend_existing=True)

        # Construct the new StepEntry class with the new table and reset the relationship with SimulationEntry
        newcls=type("StepEntry", (BaseStepEntry, Base), {"__table__": tab, "simulation": relationship(SimulationEntry)})

        # Set the module definition of StepEntry to the new class
        setattr(sys.modules[__name__], "StepEntry", newcls)
        return newcls

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
                    getattr(cls, variable)
                except AttributeError:
                    changed=True
                    conn=engine.connect()
                    conn.execute("alter table steps add column %s %s null" % (variable, strtypes[type(data[variable][0])]))
                    conn.close()
        if changed:
            cls=cls.Factory()

        # Generate a new simulation from the parameters read in from that file
        if sim is None:
            sim=SimulationEntry.from_params(session=session, **{attr: data.getncattr(attr) for attr in data.ncattrs()})
        else:
            for attr in data.ncattrs():
                if attr != "params":
                    setattr(sim, attr, data.getncattr(attr))
                else:
                    sim=sim.add_params(data.getncattr("params"))

        session.add(sim)

        # If reset is set, delete all entries in the database matching the current file
        if reset:
            for entry in session.query(StepEntry).filter(StepEntry.file==os.path.abspath(file_name)).all():
                session.delete(entry)

        # Irrelevant, but the code breaks without this line
        session.flush ()
        # session.query(StepEntry).filter(StepEntry.file==os.path.abspath(file_name)).count()

        # Add the lines from the code as new entries
        entries=[]
        for time in range(times):
            entries.append(cls())
            entries[-1].file=os.path.abspath(file_name)
            entries[-1].fid=time
            entries[-1].simulation=sim
            for variable in data.variables:
                if len(data[variable].dimensions) == 1:
                    setattr(entries[-1], variable, data[variable][time])
            session.add(entries[-1])

        session.flush()

        return entries

# Check if the simulation table exists in the database
if "simulations" not in Base.metadata.tables:
    # If not, create one with columns id and hash
    sqlalchemy.Table("simulations", Base.metadata, sqlalchemy.Column("id", sqlalchemy.Integer, primary_key=True), sqlalchemy.Column("hash", sqlalchemy.String, unique=True, onupdate=BaseSimulationEntry.calculate_hash))
    Base.metadata.tables["simulations"].create(engine)

# Check if the steps database exists in the database
if "steps" not in Base.metadata.tables:
    # If not, create one with columns id, file, fid (fileline id), and simulation_id, and set that file and fid should be unique together
    sqlalchemy.Table("steps", Base.metadata, sqlalchemy.Column("id", sqlalchemy.Integer, primary_key=True), sqlalchemy.Column("file", sqlalchemy.String), sqlalchemy.Column("fid", sqlalchemy.Integer), sqlalchemy.Column("simulation_id", sqlalchemy.ForeignKey("simulations.id"), nullable=False), sqlalchemy.Column("step", sqlalchemy.Integer), sqlalchemy.UniqueConstraint("file", "fid", name="unique_file_fid"))
    Base.metadata.tables["steps"].create(engine)

# Define SimulationEntry and StepEntry here for the first time
SimulationEntry=BaseSimulationEntry.Factory()
StepEntry=BaseStepEntry.Factory()

# An example use
Session=sqlalchemy.orm.sessionmaker(bind=engine)
session=Session()
# new=StepEntry.from_file (session, "../sims/pisces_test/output/stat_00_00.cdf")
sim=SimulationEntry()
sim=sim.add_params("../sims/pisces_test/config_copy.yaml")
for column in sim.__table__.columns:
    print (column.key, getattr (sim, column.key))
session.add(sim)
sim.dump__every=100
sim.time__steps=200
sim.run(cwd="../sims/pisces_test", reset=True, debug=True)
session.commit()
print("HASH IS", sim.hash)
sim.time__steps+=100
sim.run(cwd="../sims/pisces_test", reset=True)
# # session.rollback()
session.commit()

# print(session.query(SimulationEntry).first().steps)

# print (new[0].file)
