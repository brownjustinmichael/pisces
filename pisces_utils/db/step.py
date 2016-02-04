import os
import datetime

import numpy.ma as ma
import sqlalchemy
import netCDF4

import pisces_utils.config as config
from .db import Base, engine, strtypes
from .simulation import SimulationEntry

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
                sqlalchemy.Column("date", sqlalchemy.DateTime), 
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
    def from_file (cls, session, file_name, sim=None, reset=True):
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
            var = ma.filled(data[variable], 0.0)
            if len(var.shape) == 1:
                times=var.shape[0]
                try:
                    getattr(cls.Table, variable)
                except AttributeError:
                    changed=True
                    if session.dirty:
                        raise RuntimeError("Dirty session")
                    conn=engine.connect()
                    print("TABLE CHANGE: ADDING COLUMN %s TO STEPS" % variable)
                    conn.execute("alter table steps add column %s %s null" % (variable, strtypes[type(var[0])]))
                    conn.close()
        if changed:
            cls.Table=type("StepTable", (Base,), {"__table__": sqlalchemy.Table("steps", Base.metadata, autoload=True, autoload_with=engine, extend_existing=True), "simulation": sqlalchemy.orm.relationship(SimulationEntry.Table)})

        # Generate a new simulation from the parameters read in from that file
        if sim is None:
            sim=SimulationEntry.from_config(session=session, **{attr: data.getncattr(attr) for attr in data.ncattrs()})
        else:
            for attr in data.ncattrs():
                if attr != "params":
                    setattr(sim, attr, data.getncattr(attr))
                else:
                    configuration = config.Configuration(data.getncattr("params"))
                    for param in configuration:
                        setattr(sim, param, configuration[param])

        date = datetime.datetime.fromtimestamp (os.path.getmtime(file_name))
        # If reset is set, delete all entries in the database matching the current file
        previous = session.query(cls).filter(cls.file == os.path.abspath(file_name)).order_by(cls.line).all()
        if len(previous) > 0:
            if date > previous[0].date:
                if reset:
                    for prev in previous:
                        session.delete(prev)
            else:
                return previous

        # Irrelevant, but the code breaks without this line
        session.flush()
        # session.query(StepEntry).filter(StepEntry.file==os.path.abspath(file_name)).count()

        # Add the lines from the code as new entries
        entries=[]
        for time in range(times):
            # Otherwise, make a new entry and set its basic values
            entries.append(StepEntry.Table())
            entries[-1].file = os.path.abspath(file_name)
            entries[-1].line = time
            entries[-1].simulation_id = sim.id
            entries[-1].date = date

            for variable in data.variables:
                var = data[variable][time]
                var = ma.filled(var, 0.0).tolist()
                if len(data[variable].dimensions) == 1:
                    setattr(entries[-1], variable, var)
            session.add(entries[-1])

        # session.flush()

        print(entries)

        return entries

sqlalchemy.orm.mapper(StepEntry, StepEntry.Table.__table__)