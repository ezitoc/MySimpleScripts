#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : abstractClasses.py
#
#* Purpose :
#
#* Creation Date : 18-02-2013
#
#* Last Modified : Mon 18 Feb 2013 10:31:18 AM ART
#
#* Created By :  Ezequiel Castillo
#
#_._._._._._._._._._._._._._._._._._._._._.

class AbstractClass(object):

    def first(self):
       raise NotImplementedError

    def second(self):
       raise NotImplementedError

class ConcreteClass(AbstractClass):

    def first(self):
       print self.__class__.__name__

if __name__ == "__main__":

    obj = ConcreteClass()

    try:
       obj.first()
       obj.second()
    except NotImplementedError:
       print "Not implemented.  Good to know."
