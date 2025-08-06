import time

class Data:
    def __init__(self):
        self.funcs = {}
        self.calls = {}
    
    def show(self):
        print("TIME:", self.funcs)
        print("NUMBER OF CALLS:", self.calls, flush=True)


data = Data()


def get_time(func):
    def inside(*args):
        t0 = time.time()
        res = func(*args)
        if not func.__name__ in data.funcs:
            data.funcs[func.__name__] = 0
            data.calls[func.__name__] = 0
        data.funcs[func.__name__] += time.time() - t0
        data.calls[func.__name__] += 1
        return res
    return inside


@get_time
def test():
    time.sleep(1)

if __name__ == "__main__":
    test()
    print(data.funcs)