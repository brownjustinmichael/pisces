import time
import subprocess

from db.launch import Launcher, CodeRegistry, LauncherRegistry

class TorqueLauncher(Launcher):
    """
    A class designed to launch code execution through Torque.

    Since Torque is its own queuing system, it doesn't make sense 
    """
    def __init__(self, code, sub_launcher="Launcher"):
        super().__init__(code)
        self.code = code
        self.sub_launcher = sub_launcher

    def _launch(self, *args, **kwargs):
        cwd = os.getcwd()
        os.chdir(self.code.wd)

        batch_file = open("batch.pbs" % guess, "w")

        batch_file.write("#PBS -S /bin/bash\n")
        batch_file.write("#PBS -q normal\n")
        batch_file.write("#PBS -N %s\n" % commandRoot)

        batch_file.write("#PBS -l nodes=%d:ppn=%d\n" % (self.code.np, self.code.threads))
        batch_file.write("#PBS -l walltime=%02i:00:00\n" % hours)
        batch_file.write("cd $PBS_O_WORKDIR\n")
        batch_file.write("cp $PBS_NODEFILE .\n")
        batch_file.write("export HOSTFILE=hostfile\n")
        batch_file.write("export I_MPI_PIN_DOMAIN=omp")
        batch_file.write("export KMP_AFFINITY=compact")
        
        batch_file.write("module load python\n")
        batch_file.write("module switch python python/3.4.1\n")
        batch_file.write("python3 host_rewrite.py $PBS_NODEFILE --ppn %d > $HOSTFILE\n" %(self.code.threads))
        batch_file.write("export OMP_NUM_THREADS=%d\n" % self.code.threads)

        batch_file.write(" ".join(self.code.call()))
        batch_file.write("\n")
        batch_file.close()
        
        p = subprocess.Popen(["qsub", "batch.pbs"], stdout=subprocess.PIPE)
        output = p.communicate()
        print(output.split(".")[0])
        self.process = output.split(".")[0]

        os.chdir(cwd)

    def _wait(self, *args, **kwargs):
        output = "False"
        while output != "":
            p = subprocess.Popen(["qstat", "|", "grep", self.process], stdout=subprocess.PIPE)
            output = p.communicate()
            time.sleep (5)