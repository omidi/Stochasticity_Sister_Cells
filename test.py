import threading
import time
import numpy as np

a = np.zeros(10)

class myThread (threading.Thread):
    def __init__(self, threadID, index, a):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.index = index
        self.a = a
    def run(self):
        # Get lock to synchronize threads
        threadLock.acquire()
        self.a[self.index] = 2**self.index
        # Free lock to release next thread
        threadLock.release()


threadLock = threading.Lock()
threads = []

# Create new threads
for i in xrange(len(a)):
    threads.append(myThread(i, i, a))
    threads[i].start()

# Wait for all threads to complete
for t in threads:
    t.join()

print a