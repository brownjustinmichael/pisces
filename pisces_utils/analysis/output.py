import os
from subprocess import call

import netCDF4 as nc
import numpy as np
import numpy.ma as ma
from numpy import meshgrid
from scipy.interpolate import LinearNDInterpolator
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt

class Output(nc.Dataset):
    """An output object for analysis"""
    def __init__(self, file, *args, **kwargs):
        super(Output, self).__init__(file, *args, **kwargs)

    @classmethod
    def from_simulation_entry(cls, simulation_entry, **kwargs):
        files = simulation_entry.get_files(**kwargs)
        outputs = []
        for file in files:
            outputs.append(cls(file))
        return outputs

    def __len__(self):
        time_key = list(self.dimensions.keys())[0]
        return len(self.dimensions[time_key])

    def norm(self, other, variable, self_index=-1, other_index=-1, n=2, positions=("x", "z")):
        """
        Calculate the Ln norm of the file with another. If the two are not the same size, self will be interpolated via nearest neighbors to the points of other.

        :type other: :class:`Output`
        :param other: The other output object from which the difference should be calculated
        :type variable: :class:`str`
        :param variable: The variable key to use for the norm calculation
        :type index: :class:`int`
        :param index: The index to use for the norm calculation
        :type n: :class:`int`
        :param n: The order of the norm
        :type positions: :class:`tuple` of :class:`str`
        :param positions: A tuples of variable keys to use as positional arguments to interpolate if necessary

        :return: The relative Ln norm difference between the two outputs
        """
        try:
            diff = (self.variables[variable][self_index] - other.variables[variable][other_index])
        except ValueError:
            flattened = [self.variables[var][self_index].flatten() for var in positions]
            points = [i for i in zip(*flattened)]
            values = self.variables[variable][self_index].flatten()

            inter_flattened = [other.variables[var][other_index].flatten() for var in positions]
            inter_points = [i for i in zip(*inter_flattened)]

            interpolator = LinearNDInterpolator(points, values)
            new_values = np.array([interpolator(point) for point in inter_points])
            diff = (new_values - other.variables[variable][other_index].flatten())

        scale = np.sum(other.variables[variable][other_index])
        return (np.sum(diff ** n)) ** (1.0 / n) / scale

    def plot(self, ckey, xkey="x", ykey="z", ax=None, index=-1, **kwargs):
        x = ma.getdata(self[xkey][index])
        y = ma.getdata(self[ykey][index])
        c = ma.getdata(self[ckey][index])

        if ax is None:
            fig, ax = plt.subplots(1, 1)
        
        p = ax.pcolormesh(x.T, y.T, c.T, **kwargs)

        ax.set_xlim((x [0,0], x [-1,-1]))
        ax.set_ylim((y [0,0], y [-1,-1]))
        
        ax.set_xlabel("$%s$" % xkey)
        ax.set_ylabel("$%s$" % ykey)
        ax.set_title("$%s$" % ckey)
        ax.set_aspect(1)
        
        return p

    def quiver(self, ukey="x_velocity", vkey="z_velocity", xkey="x", ykey="z", dx=0.1, dy=0.1, ax=None, index=-1, **kwargs):
        x = ma.getdata(self[xkey][index])
        y = ma.getdata(self[ykey][index])

        if ax is None:
            fig, ax = plt.subplots(1, 1)
        
        X, Y = meshgrid (np.arange (np.min (x), np.max (x), dx), np.arange (np.min (y), np.max (y), dy))

        u = ma.getdata(self.variables[ukey][index])
        w = ma.getdata(self.variables[vkey][index])
        U = griddata (x.flatten (), y.flatten (), u.flatten (), X, Y, interp = "linear")
        W = griddata (x.flatten (), y.flatten (), w.flatten (), X, Y, interp = "linear")    
        
        v = ax.quiver (X, Y, U, W, **kwargs)

        ax.set_xlim ((x [0,0], x [-1,-1]))
        ax.set_ylim ((y [0,0], y [-1,-1]))
        
        ax.set_xlabel ("$%s$" % xkey)
        ax.set_ylabel ("$%s$" % ykey)
        ax.set_aspect(1)
        
        return v

    def contour(self, ckey, levels=10, xkey="x", ykey="z", ax=None, index=-1, **kwargs):
        x = ma.getdata(self[xkey][index])
        y = ma.getdata(self[ykey][index])
        c = ma.getdata(self[ckey][index])

        if ax is None:
            fig, ax = plt.subplots(1, 1)
        
        p = ax.contour(x.T, y.T, c.T, levels, **kwargs)

        ax.set_xlim ((x [0,0], x [-1,-1]))
        ax.set_ylim ((y [0,0], y [-1,-1]))
        
        ax.set_xlabel ("$%s$" % xkey)
        ax.set_ylabel ("$%s$" % ykey)
        ax.set_title ("$%s$" % ckey)
        ax.set_aspect(1)
        
        return p

    def record(self, index=-1):
        return Output.Record(self, index)

    class Record:
        def __init__(self, output, index=-1):
            self.output = output
            self.index = index

        def contour(self, *args, **kwargs):
            return self.output.contour(*args, index=self.index, **kwargs)

        def plot(self, *args, **kwargs):
            return self.output.plot(*args, index=self.index, **kwargs)

        def quiver(self, *args, **kwargs):
            return self.output.quiver(*args, index=self.index, **kwargs)

class Movie:
    def __init__(self, movie_file="movie.mp4", file_format="tmp_%04d.png", frame_rate=5, remove_files=True):
        self.movie_file = movie_file
        self.file_format = file_format
        self.files = []
        self.frame_rate = frame_rate
        self.remove_files = remove_files
        if self.movie_file.split(".")[-1] == "mp4":
            self.vcodec = "mpeg4"
        else:
            raise RuntimeError("Unknown movie extension %s" % self.movie_file.split(".")[-1])

    def plots(self, *args, **kwargs):
        i = 0
        while(True):
            fig, ax = plt.subplots(*args, **kwargs)
            yield fig, ax
            fig.savefig(self.file_format % i)
            self.files.append(self.file_format % i)
            plt.close()
            i += 1

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        call(["ffmpeg", "-r", str(self.frame_rate), "-i", self.file_format, "-vcodec", str(self.vcodec), "-y", self.movie_file])
        if self.remove_files:
            for file in self.files:
                os.remove(file)

