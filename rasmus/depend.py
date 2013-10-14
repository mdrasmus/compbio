import datetime
import os
import sys
import time
import threading



# states for jobs
STATUS_UNDONE   = "incomplete"
STATUS_PENDING  = "pending"
STATUS_RUNNING  = "running"
STATUS_DONE     = "done"
STATUS_ERROR    = "error"

VALID_STATUS = {
    STATUS_UNDONE  :1,
    STATUS_PENDING :1,
    STATUS_RUNNING :1,
    STATUS_DONE    :1,
    STATUS_ERROR   :1
}


DEFAULT_DISPATCH = "bash $SCRIPT"
BASH_DISPATCH = "bash $SCRIPT"
LSF_DISPATCH = "bsub -o $STATUSDIR/$JOBNAME.output -K bash $SCRIPT"


# autodetect dispatch
def getDefaultDispatch():
    platform = getPlatform()
    
    if platform == "lsf":
        return LSF_DISPATCH
    else:
        return DEFAULT_DISPATCH

def getDefaultBackground():
    platform = getPlatform()
    
    if platform == "lsf":
        return True
    else:
        return False



class PipelineException (Exception):
    def __init__(self, msg):
        Exception.__init__(self, "pipeline: " + msg)
    


class Job:
    def __init__(self, name, task, depends=[], 
                 background=True, dispatch=DEFAULT_DISPATCH):
        
        self.pid = -1
        self.name = name
        self.task = task
        self.background = background        
        self.status = STATUS_UNDONE
        self.msg = ""
        self.parents = []
        self.children = []
        self.waitSet = {}
        self.dispatch = dispatch
        self.subjobs = []

        
        # add dependencies
        for dep in depends:
            self.addDep(dep)        
        
        # determine task
        if isinstance(task, str):
            self.tasktype = "shell"
        else:
            self.tasktype = "function"
    
    def undo(self):
        self.status = STATUS_UNDONE
        for child in self.children:
            child.undo()
    
    def addDep(self, dep):
        self.parents.append(dep)
        dep.children.append(self)
    
    def wait(self, job):
        self.waitSet[job] = 1

    def unwait(self, job):
        if job in self.waitSet:
            del self.waitSet[job]
        
    def setWaits(self):
        for parent in self.parents:
            if parent.status != STATUS_DONE:
                self.wait(parent)
            
            if parent.status == STATUS_ERROR:
                parent.raiseError()
    
    def notifyChildren(self):
       for child in self.children:
           child.unwait(self)
    
    def isWaiting(self):
        return len(self.waitSet) > 0
    
    
    def raiseError(self):
        raise PipelineException("job '%s' has error" % self.name)

 



class Pipeline:
    def __init__(self, 
                 statusDir="pipeline", 
                 background=None,
                 dispatch=None):
        
        self.statusDir = statusDir
        self.jobs = {}      # set of all jobs registered with Pipeline
        self.pending = {}   # set of all jobs that are ready to run
        self.pids = {}      # set of all process ids currently running
        self.isInit = False
        self.testing = False
        self.logOutput = None
        self.maxNumProc = 40 # maximum allowed number of processes
        self.needReset = False
        
        if background == None:
            self.background = getDefaultBackground()
        else:
            self.background = background
        
        if dispatch == None:
            self.dispatch = getDefaultDispatch()
        else:
            self.dispatch = dispatch

        
    
    
    def init(self):
        # set all job states to UNDONE
        if self.needReset:
            for job in self.jobs.itervalues():
                filename = self.getJobStatusFile(job)
                if os.path.exists(filename):
                    os.remove(filename)
            self.needReset = False
        
        # read in the status of all jobs
        self.readStatus(False)
        
        if not self.isInit:
            # clear any pending jobs
            self.pending = {}
            
            for job in self.jobs.itervalues():
                # job that were running are now back to pending
                if job.status == STATUS_RUNNING or \
                   job.status == STATUS_PENDING or \
                   job.status == STATUS_ERROR:
                    self.writeJobStatus(job, STATUS_UNDONE)
            
            self.isInit = True
    
    def reinit(self):
        self.isInit = False
        self.init()
    
    def reset(self):
        self.needReset = True
    
    def enableTesting(self, enable=True):
        self.testing = enable

    def setLogOutput(self, out=sys.stdout):
        self.logOutput = out
    
    def setStatusDir(self, statusDir):
        assert not self.isInit
        self.statusDir = statusDir
    
    def ensureStatusDir(self):
        if not os.path.exists(self.statusDir):
            os.mkdir(self.statusDir)
    
    def setMaxNumProc(self, nproc):
        self.maxNumProc = nproc
    
    
    def getJobStatusFile(self, job):
        return os.path.join(self.statusDir, job.name + ".status")
    
    def getJobOutputFile(self, job):
        return os.path.join(self.statusDir, job.name + ".output")
        
    def getJobErrorFile(self, job):
        return os.path.join(self.statusDir, job.name + ".error")
    
    def getJobScriptFile(self, job):
        return os.path.join(self.statusDir, job.name + ".script")
    
    
    def readStatus(self, retry=True):
        """Read in status information for all jobs"""
        
        for job in self.jobs.itervalues():
            self.readJobStatus(job, retry)
    
    def writeStatus(self):
        """Write status information for all jobs"""
        
        for job in self.jobs.values():
            self.writeJobStatus(job)
    
    
    def readJobStatus(self, job, retry=True):
        filename = self.getJobStatusFile(job)
        
        while True:
            if os.path.exists(filename):
                infile = file(filename, "r")
                job.status = infile.read().rstrip()
            else:
                job.status = STATUS_UNDONE
            
            # handle the case where the status file is only partially written
            if job.status in VALID_STATUS or not retry:
                break
            else:
                # wait for a little 
                time.sleep(.05)
        
        if job.status not in VALID_STATUS:
            self.writeJobStatus(job, STATUS_UNDONE)
    
    
    def writeJobStatus(self, job, status=None):
        self.ensureStatusDir()
        
        if status != None:
            job.status = status
        
        out = file(self.getJobStatusFile(job), "w")
        out.write(job.status)
        out.close()
    
    
    def log(self, *text):
        if self.logOutput:
            self.logOutput.write("pipeline: " + " ".join(text) + "\n")
    
    
    def add(self, name, task, depends=[], background=None, dispatch=None):
        # set defaults
        if background == None:
            background = self.background
        if dispatch == None:
            dispatch = self.dispatch
    
        parents = []
        for dep in depends:
            try:
                parents.append(self.jobs[dep])
            except KeyError:
                raise PipelineException("unknown job '%s'" % dep)
    
        self.jobs[name] = Job(name, task, parents, background, dispatch)
        
        return name
    
    
    def addGroup(self, name, subjobnames, depends=[], 
                 background=None,
                 dispatch=None):
        # set defaults
        if background == None:
            background = self.background
        if dispatch == None:
            dispatch = self.dispatch
        
        parents = []
        for dep in depends:
            try:
                parents.append(self.jobs[dep])
            except KeyError:
                raise PipelineException("unknown job '%s'" % dep)
        
        tasks = []
        for jobname in subjobnames:
            try:
                job = self.jobs[jobname]
            except KeyError:
                raise PipelineException("unknown job '%s'" % jobname)
            
            assert isinstance(job.task, str), \
                "subjob '%s' task must be a command line string" % jobname
            tasks.append(job.task)
        task = " && ".join(map(lambda t: "( %s )" % t, tasks))
        
        self.jobs[name] = Job(name, task, parents, background, dispatch)
        self.jobs[name].subjobs = map(lambda x: self.jobs[x], subjobnames)
        
        return name
    
    
    def addGroups(self, name, subjobnames, size=1, depends=[], 
                  background=None, dispatch=None):
        # set defaults
        if background == None:
            background = self.background
        if dispatch == None:
            dispatch = self.dispatch
        
        groups = []
        j = 1
        for i in xrange(0, len(subjobnames), size):
            groups.append(self.addGroup("%s%d" % (name, j),
                                        subjobnames[i:i+size],
                                        depends,
                                        background,
                                        dispatch))
            j += 1
        return groups
    
    
    def run(self, name):
        self.init()
        
        if name in self.jobs:
            self.runJob(self.jobs[name])
        else:
            raise PipelineException("unknown job '%s'" % name)
    
    
    def undo(self, name):
        if name in self.jobs:
            self.undoJob(self.jobs[name])
        else:
            raise PipelineException("unknown job '%s'" % name)
    
    
    def addPending(self, job):
        self.writeJobStatus(job, STATUS_PENDING)
        self.pending[job] = 1
    
    
    def removePending(self, job, status=None):
        if job in self.pending:
            del self.pending[job]
            
            if status != None:
                job.status = status
                self.writeJobStatus(job, job.status)
        
        assert job.status != STATUS_PENDING
    
    
    def finishJob(self, job):
        self.log("%s: END" % job.name)
        
        self.removePending(job)
        if job.pid != -1:
            retpid, retcode = os.waitpid(job.pid, 0)
            self.finishPid(job.pid)
        job.notifyChildren()
    
    def finishPid(self, pid):
        job = self.pids[pid]
        del self.pids[pid]
        job.pid = -1

    
    def runJob(self, job): 
        # return immediately if job is already done
        if job.status == STATUS_DONE:
            self.finishJob(job)
            return STATUS_DONE
        
        # do not job if it has an error
        elif job.status == STATUS_ERROR:
            job.raiseError()
        
        # do not run job if it is already running
        elif job.status == STATUS_RUNNING:
            return STATUS_RUNNING
        
        # ensure these are the only valid job states
        assert job.status == STATUS_UNDONE or \
               job.status == STATUS_PENDING, \
               "unknown job status '%s'" % job.status
        
        
        
        # determine which jobs to wait for
        job.setWaits()
        
        # add job to pending jobs
        if job not in self.pending:
            self.addPending(job)
            
            # make sure all parents are run first
            for parent in job.parents:
                self.runJob(parent)
        
        # run job if it is waiting for no one
        # and number of processes is less than max allowed
        if not job.isWaiting() and \
           len(self.pids) < self.maxNumProc:
            self.execJob(job)
        
        # return job status
        return job.status
    
    
    def execJob(self, job):
        self.log("%s: BEGIN" % job.name)
    
        # mark job as running
        self.writeJobStatus(job, STATUS_RUNNING)
    
        if job.tasktype == "function":
            # currently functions can not be backgrounded
            #
            
            if self.testing:
                print "* running job '%s' (python function)\n" % job.name
                self.writeJobStatus(job, STATUS_DONE)
                return
            
            # run task in main thread (BLOCK)
            if job.task():
                self.writeJobStatus(job, STATUS_DONE)
            else:
                self.writeJobStatus(job, STATUS_ERROR)
            
        elif job.tasktype == "shell":
            if self.testing:
                print "* running job '%s':\n%s\n" % (job.name, job.task)
                self.writeJobStatus(job, STATUS_DONE)
                return
            
            
            if job.background:
                # run task in separate process and dont block
                
                # save task into script file
                script = self.getJobScriptFile(job)
                out = file(script, "w")
                out.write(job.task)
                out.close()
                
                # expand dispatch
                dispatch = job.dispatch
                dispatch = dispatch.replace("$JOBNAME", job.name)
                dispatch = dispatch.replace("$SCRIPT", script)
                dispatch = dispatch.replace("$STATUSDIR", self.statusDir)
                
                statusfile = self.getJobStatusFile(job)
                job.pid = os.spawnlp(os.P_NOWAIT, "bash", "bash", "-c", 
                    """( %s ) && 
                       (echo done > '%s') ||
                       (echo error > '%s')
                    """ % 
                    (   
                        dispatch,
                        statusfile,
                        statusfile
                    ))

                # save process id
                self.pids[job.pid] = job
            else:
                # run task in separate process and wait for it
                
                if os.system(job.task) == 0:
                    self.writeJobStatus(job, STATUS_DONE)
                else:
                    self.writeJobStatus(job, STATUS_ERROR)
            
        else:
            raise PipelineException("unknown tasktype '%s'" % job.tasktype)
    
    
    def undoJob(self, job):
        """Make job and all depending jobs as 'incomplete'"""
        
        self.writeJobStatus(job, STATUS_UNDONE)
        for child in job.children:
            self.undoJob(child)
    
    
    def process(self, poll=False):
        self.init()
        
        while len(self.pending) > 0:
            # try to run all pending jobs
            for job in self.pending.keys():
                self.runJob(job)
        
            # wait for jobs to finish
            while len(self.pids) > 0:
                if not poll:
                    pid, retcode = os.waitpid(0, 0)                
                else:
                    # do not hang if only polling
                    pid, retcode = os.waitpid(0, os.WNOHANG)
                
                if pid in self.pids:
                    job = self.pids[pid]
                    self.finishPid(pid)
                    self.readJobStatus(job)

                    # these are the only allowed states
                    assert job.status == STATUS_DONE or \
                           job.status == STATUS_ERROR

                    if job.status == STATUS_DONE:
                        self.finishJob(job)

                        # try execute children
                        for child in job.children:
                            if child in self.pending:
                                self.runJob(child)

                    elif job.status == STATUS_ERROR:
                        job.raiseError()
                    
                    
                    # try to run all pending jobs
                    for job in self.pending.keys():
                        self.runJob(job)
                else:
                    print >>sys.stderr, "ERROR: unknown job with pid %d reporting" % \
                                            pid
                
                # do not loop if only polling
                if poll:
                    return


def hasLsf():
    """Returns True only if LSF is available"""

    ret = os.system("which bsub > /dev/null")
    return (ret == 0)

PLATFORM = None

def getPlatform():
    global PLATFORM

    if PLATFORM == None:
        if hasLsf():
            PLATFORM = "lsf"
        else:
            PLATFORM = "none"
    return PLATFORM


def submit(cmd):
    if getPlatform() == "lsf":
        return "bsub -o /dev/null -K " + cmd
    else:
        return cmd



# testing
if __name__ == "__main__":
    pipeline = Pipeline()
    
    pipeline.add("job1", "echo job1")
    pipeline.add("job2", "echo job2 && sleep 3")
    pipeline.add("job3", "echo job3", ["job1", "job2"])
    pipeline.add("job4", "cat nosuchthing", ["job3"])
    pipeline.add("job5", "echo job5", ["job1", "job3"])
    
    #pipeline.undo("job1")
    
    pipeline.run("job5")
    
    pipeline.process()

