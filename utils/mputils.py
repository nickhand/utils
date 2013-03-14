import multiprocessing as mp
import logging, os, sys
progressLoaded = True
try:
    from utils.utilities import initializeProgressBar
except:
    progressLoaded = False

class worker(mp.Process):
    """
    @brief worker class that a dequeues a task from an input queue, perform a specified
    computation, and store the results until the input queue is empty
    """
    
    def __init__(self, task_queue, result_queue, pbar=None):
        
        # initialize a new Process for each worker
        mp.Process.__init__(self) 
        
        # save the task and results queue
        self.task_queue   = task_queue 
        self.result_queue = result_queue 
        
        if pbar is not None:
            self.pbar = pbar
        else:
            self.pbar = None
            
        return
        

    def run(self):
        """
        @brief start the worker class doing the tasks until there
        are none left
        """
        
        # pull tasks until there are none left
        while True:
            # dequeue the next task
            next_task = self.task_queue.get()
            
            # task = None means this worker is finished
            if next_task is None:
                # make sure we tell the queue we finished the task
                self.task_queue.task_done()
                break
                
            # try to update the progress bar
            if self.pbar is not None:
                try: 
                    self.pbar.update(next_task.num+1)
                except:
                    raise
                
            # do the work by calling the task    
            answer = next_task()
            
            # store the answer
            self.result_queue.put(answer)
            
            # make sure we tell the queue we finished the task
            self.task_queue.task_done()
        
        return 0
    
class task(object):
    """
    @brief a class representing a 'task' where a specified computation
    is performed
    """
    
    def __init__(self, function, *args, **kwargs):
        
        self.func = function
        self.args = args
        self.num = kwargs.get('num', 0)
        
        
    def __call__(self):
        
        # call the function with the arguments
        ans = self.func(*self.args)
            
        return ans



class mp_master(object):
    """
    @brief a class to control a multiprocessing job 
    """
    
    def __init__(self, nprocs, njobs, progress=True):
        """
        @brief initialize the input/output queues and make the workers
        """
        
        # set up the queues
        self.results = mp.Queue()
        self.tasks = mp.JoinableQueue()
        
        # set up the logger to log to sys.stderr
        self.logger = mp.log_to_stderr()
        self.logger.setLevel(logging.INFO)
        
        # if we want a progress bar
        if progress and progressLoaded:
            bar = initializeProgressBar(njobs)
        else:
            bar = None
        
        # start a worker for each cpu available
        print 'creating %d workers' % nprocs
        self.workers = [ worker(self.tasks, self.results, pbar=bar) for i in range(nprocs) ]
        
        return
    
    def enqueue(self, task):
        """
        @brief enqueue a task onto the tasks queue
        """
        
        self.tasks.put(task)
        
        return
        
    def run(self):
        """
        @brief start the workers and do the work
        """
        
        try: 
            # start the work
            for w in self.workers:
                w.start()
            
            # add a poison pill for each worker
            for i in range(len(self.workers)):
                self.tasks.put(None)
            
            # wait for all tasks to finish
            self.tasks.join()
            
        except KeyboardInterrupt:
            
            print "caught keyboard interrupt..", mp.active_children()
            for w in self.workers:
               
               w.terminate()
               w.join()
            
        
        return
        
    def dequeue(self):
        """
        @brief dequeue the results, if available, else None
        """
        
        return self.results.get()

    def more_results(self):
        """
        @brief return True if there are more results to dequeue
        """
        
        return not self.results.empty()
    
