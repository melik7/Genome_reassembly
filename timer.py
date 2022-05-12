import time
# This class allows us to monitor the time spent by a function
class Timer():
    # Function called right before the execution of the function
    def __enter__(self):
        self.t1 = time.perf_counter()
        return self
    # Function called right after the function
    def __exit__(self, type, value, traceback):
        self.t2 = time.perf_counter()
        self.t = self.t2 - self.t1
    # Function that prints on the shell the time spent by the instructions
    def print(self, template: str = "{}"):
        print(template.format(round(self.t, 2)))

