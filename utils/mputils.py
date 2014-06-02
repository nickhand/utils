import multiprocessing
from multiprocessing.queues import Queue as mpQueue
import logging, os, sys
import tempfile, datetime
from utils import utilities
import os

progressLoaded = True
try:
    from utils.utilities import initializeProgressBar
except ImportError:
    progressLoaded = False

#-------------------------------------------------------------------------------
class Counter(object):
    
    def __init__(self, init_val=0):
        self.lock = multiprocessing.Lock()
        self.val = multiprocessing.Value('i', init_val)
    
    def increment(self):
        with self.lock:
            self.val.value += 1
    
    def decrement(self):
        with self.lock:
            self.val.value -= 1

    def value(self):
        with self.lock:
            return self.val.value

#-------------------------------------------------------------------------------
class Queue(mpQueue):
    """
    Queue-like class with overflow limits
    """
    MAX_QUEUE_SIZE = 30000
    
    def __init__(self, *args, **kwargs):
        
        super(Queue, self).__init__(*args, **kwargs)
        
        self.queue_size = Counter()
        self.overflow = []
        
    #end __init__
    
    #---------------------------------------------------------------------------
    @property 
    def size(self):
        return self.queue_size.value() + self.overflow_size
        
    @property
    def overflow_size(self):
        return len(self.overflow)
    
    #---------------------------------------------------------------------------
    def get(self, **kwargs):
        """
        Dequeue and return an object
        """
        if self.empty():
            raise multiprocessing.queues.Empty("Cannot dequeue from an empty Queue object")
        
        if self.queue_size.value() > 0:
            self.queue_size.decrement()
            return super(Queue, self).get(**kwargs)
        else:
            return self.overflow.pop()
    #end get
    
    #---------------------------------------------------------------------------
    def empty(self): 
        return (self.size == 0)
    
    #end empty
    
    #---------------------------------------------------------------------------
    def put(self, obj, **kwargs):
        """
        Enqueue onto the queue
        """
        if self.queue_size.value() + 1 <= Queue.MAX_QUEUE_SIZE:
            super(Queue, self).put(obj, **kwargs)
            
            print "queue size = ", self.queue_size.value()
            self.queue_size.increment()
        else:
            self.overflow.append(obj)
    #end put
    
    #---------------------------------------------------------------------------
#endclass Queue

#-------------------------------------------------------------------------------
class worker(multiprocessing.Process):
    """
    Worker class that a dequeues a task from an input queue, perform a specified
    computation, and store the results until the input queue is empty
    """
    
    def __init__(self, task_queue, result_queue, except_event, pbar=None):
        
        # initialize a new Process for each worker
        multiprocessing.Process.__init__(self)
        
        
        # save the task and results queue
        self.task_queue   = task_queue 
        self.result_queue = result_queue 
        
        # handle the progress bar
        if pbar is not None:
            self.pbar = pbar
        else:
            self.pbar = None

        # handle an exception
        self.exception = except_event
    #end __init__
    
    #---------------------------------------------------------------------------
    def run(self):
        """
        Start the worker class doing the tasks until there are none left
        """
        # pull tasks until there are none left and we don't exit
        while not self.exception.is_set():
            
            # dequeue the next task
            next_task = self.task_queue.get()
            print "task queue size = ", self.task_queue.size, os.getpid()
            
            # task == None means we should exit
            if next_task is None:
                break
            
            # try to update the progress bar
            if self.pbar is not None:
                try: 
                    self.pbar.update(next_task.num+1)
                except:
                    self.exception.set()
                    raise Exception("Exception event in multiprocessing")
                    
            # try to do the work
            try:  
                answer = next_task()
                self.result_queue.put(answer)
                print "results queue size = ", self.result_queue.queue_size.value(), os.getpid()
                
            # set the exception event so main process knows to exit, 
            # and then raise the exception
            except:
                self.exception.set()
                raise Exception("Error trying to perform task.")
    #end run
#endclass worker    

#-------------------------------------------------------------------------------
class task(object):
    """
    A class representing a 'task' where a specified computation is performed
    """
    def __init__(self, function, *args, **kwargs):
        
        self.func = function
        self.args = args
        self.num = kwargs.get('num', 0)
        
        
    def __call__(self):
        
        return self.func(*self.args)
#endclass task

#-------------------------------------------------------------------------------
class mp_master(object):
    """
    A class to control a multiprocessing job 
    """
    
    def __init__(self, nprocs, njobs, progress=True, log=True):
        """
        Initialize the input/output queues and make the workers
        """
        
        # set up the results queue
        self.results = Queue()
        
        # set up the tasks queue
        self.tasks = Queue()

        # hold the dequeued results
        self.deqd_results = []
        
        self.log = log        
        if self.log:
        
            # redirect stderr to a file
            self.temp_stderr = tempfile.TemporaryFile()
            sys.stderr = self.temp_stderr
        
            # make a unique file name for std out
            fileName, extension = os.path.splitext(os.path.basename(sys.argv[0]))
            time_stamp = utilities.get_timestamp()
            self.stdout = open(os.getcwd() + os.sep + "%s.%s.log" \
                                                %(fileName, time_stamp), 'w')
            sys.stdout = self.stdout
        
            # set up the logger to log to sys.stderr
            self.logger = multiprocessing.log_to_stderr()
            self.logger.setLevel(logging.INFO)
        
        # if we want a progress bar
        if progress and progressLoaded:
            bar = initializeProgressBar(njobs, fd=sys.__stderr__)
        else:
            bar = None
        
        # create an exception event
        self.exception = multiprocessing.Event()
        
        # start a worker for each cpu available
        print 'creating %d workers' %nprocs
        self.workers = [ worker(self.tasks, self.results, self.exception, pbar=bar) for i in range(nprocs) ]
        
    #end __init__

    #---------------------------------------------------------------------------
    def run(self):
        """
        Start the workers and do the work
        """
        
        # make sure to catch exceptions
        try: 
            # start the work
            for w in self.workers:
                w.start()
            
            # add a poison pill for each worker
            for i in range(len(self.workers)):
                self.enqueue(None)
                
            # while processes still alive, dequeue the results and store
            while any([w.is_alive() for w in self.workers]):
                while self.more_results():
                    self.deqd_results.append(self.results.get())
                
            print "len of deqd results = ", len(self.deqd_results)
            print "len of results = ", self.results.qsize()
            
            # if exception, raise
            if self.exception.is_set():
                raise Exception("Exception event in multiprocessing")
        except:
            
            # close all the workers gracefully
            for w in self.workers:
               w.terminate()
               w.join()
            raise Exception("Exception event in multiprocessing")
             
        finally: 
            
            # append the temp stderr to stdout file
            if self.log:
                self.stdout.write('%s\n' %('-'*100))
                self.temp_stderr.seek(0)
                self.stdout.write(self.temp_stderr.read())
                self.stdout.write('%s\n' %('-'*100))
            
            # summary
            self.info()
    #end run
        
    #---------------------------------------------------------------------------
    def enqueue(self, task):
        """
        Enqueue a task onto the tasks queue
        """
        self.tasks.put(task)
    #end enqueue
    
    #---------------------------------------------------------------------------
    def more_results(self):
        """
        Check if there are any results to dequeue
        """
        return (self.results.size > 0)
    #end more_results
    
    #---------------------------------------------------------------------------
    def dequeue(self):
        """
        Dequeue all the results
        """
        while self.more_results():
            self.deqd_results.append(self.results.get())
                
        return self.deqd_results
    #end dequeue
    
    #---------------------------------------------------------------------------
    def info(self):
        """
        Summarize process info
        """
        # print out exit codes
        for w in self.workers:
            sys.stdout.write("exit code for Process '%s' is %s\n" % (w.name, w.exitcode))
            
        # print out finish time
        now = datetime.datetime.now()
        sys.stdout.write("job finished at %s\n\n" %str(now))
    #end info
#endclass mp_master

#-------------------------------------------------------------------------------

        
    
