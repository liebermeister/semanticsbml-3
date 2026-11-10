'''
custom file class with read progress
'''
import os,stat
from . import stderr

class CustomFile:

        def __init__(self,fn,mode='rw'):
                self.f = open(fn, mode)
                self.num_lines = reduce(lambda x,y: x+1, open(fn), 0)
                #self.size = int(os.stat(fn)[stat.ST_SIZE])
                self.cur_line=0
                #self.cur_size=0
                stderr.downloadProgress(self.num_lines,0,True,url=fn)

        def read(self,val):
                #self.cur_size+=val
                #print 'read R ',self.size,self.cur_size
                #stderr.downloadProgress(self.size,self.cur_size,False,False)
                return self.f.read(val)

        def readline(self):
                self.cur_line+=1
                if self.cur_line%100000==0:stderr.downloadProgress(self.num_lines,self.cur_line,False,False)
                return self.f.readline()

        def close(self):
                stderr.downloadProgress(0,0,False,True) 
                self.f.close()
