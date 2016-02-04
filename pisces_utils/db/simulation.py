import sqlalchemy

from .db import Base, engine, strtypes
from pisces_utils import config

class SimulationEntry(object):
    if "simulations" not in Base.metadata.tables:
        # If not, create one with columns id and hash
        Table=type("SimulationTable", (Base,), {"__table__": sqlalchemy.Table("simulations", Base.metadata, sqlalchemy.Column("id", sqlalchemy.Integer, primary_key=True))})
        Base.metadata.tables["simulations"].create(engine)
    else:
        Table=type("SimulationTable", (Base,), {"__table__": sqlalchemy.Table("simulations", Base.metadata, autoload=True, autoload_with=engine, extend_existing=True)})

    def __init__(self, entry=None, *args, **kwargs):
        if entry:
            self._entry=entry
        else:
            self._entry=self.Table(*args, **config.process(kwargs))

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

        remove = []
        for key in translation:
            if key[:5] == "input":
                remove.append(key)
            if key == "output__number":
                remove.append(key)

        for key in remove:
            translation.pop(key)
            
        try:
            entry = session.query(cls).filter_by(**translation).one()
            return entry
        except sqlalchemy.orm.exc.NoResultFound:
            pass

        # No acceptable simulation has been found, so we generate one from scratch
        return cls.Table(**config.process(kwargs))

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


sqlalchemy.orm.mapper(SimulationEntry, SimulationEntry.Table.__table__)
