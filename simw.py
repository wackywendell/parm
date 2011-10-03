from sim import *

def itern(obj, n):
    for i in range(n):
        yield obj.get(i)

def stliter(obj):
    b = obj.begin()
    while b != obj.end():
        yield b.next()
