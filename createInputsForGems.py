#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : createInputsForGems.py
#
#* Purpose :
#
#* Creation Date : 29-01-2013
#
#* Last Modified : Mon 04 Feb 2013 04:45:41 PM ART
#
#* Created By :  Ezequiel Castillo
#
#_._._._._._._._._._._._._._._._._._._._._.

import os
import sys
import shutil

def pairwise(iterable):
    itnext = iter(iterable).next
    while True:
        yield itnext( ), itnext( )

def dictFromSequence(seq):
    return dict(pairwise(seq))

def replace_words(text, word_dic):
    """
    take a text and replace words that match a key in a dictionary with
    the associated value, return the changed text
    """
    rc = re.compile('|'.join(map(re.escape, word_dic)))
    def translate(match):
        return word_dic[match.group(0)]
    return rc.sub(translate, text)

def multiple_subst(infile, replace_pattern):
    # read the file
    fin = open(infile, "r")
    str1 = fin.read()
    fin.close()
    # add 'replace pattern' : 'replacement' to dictionary
    dict_make = dictFromSequence(replace_pattern)
    # call the function and get the changed text
    str2 = replace_words(str1, dict_make)
    # write changed text back out
    fout = open("run-"+infile, "w")
    fout.write(str2)
    fout.close()


class Base(object):

    def __init__(self, projectName='project', modFiles=None, filesToCopy=[]):
        self.projectName = projectName
        self.rootDir = os.path.abspath('.')
        self.ignoreFiles = modFiles
        self.filesToCopy = filesToCopy

    def selectFiles(self):
        """Return a list of files that won't be modified"""
        if not self.filesToCopy:
            #Copy all files by default.
            filesList = [f for f in os.listdir(self.rootDir) if
                         os.path.isfile(os.path.join(self.rootDir, f))]
        else:
            filesList = self.filesToCopy
        if self.ignoreFiles:
            filesListMod = [f for f in filesList if f not in self.ignoreFiles]
            return filesListMod
        else:
            return filesList

    def createFolders(self, suffixNo=None, suffixName=None, foldersName="run"):
        """Create subfolders containing all files declared with the filesToCopy
        option (all files included, otherwise) and ignoring those declared at the
        modFiles option"""

        if suffixNo and suffixName:
            print "ERROR: Two suffix types declared. Declare only one."
            sys.exit(1)

        if suffixNo or suffixName:
            if suffixNo:
                self.suffixList = range(suffixNo)
            else:
                self.suffixList = suffixName
        else:
            print "ERROR: No suffix declared."
            sys.exit(1)

        self.foldersName = foldersName
        # CWD = CWD/projectName
        if not os.path.isdir(self.projectName):
            os.mkdir(self.projectName)
        else:
            print 'WARN: folder "%s" already exists.' % (self.projectName)
        os.chdir(self.projectName)
        for suffix in self.suffixList:
            self.folderName = '%s_%s' % (self.foldersName, suffix)
        # CWD = CWD/projectName/folderName
            os.mkdir(self.folderName)
            os.chdir(self.folderName)
            copySelectedFiles = self.selectFiles()
            for file in copySelectedFiles:
                shutil.copy(os.path.join(self.rootDir, file), os.curdir)
        # CWD = CWD/projectName
            os.chdir(os.pardir)

        os.chdir(self.rootDir)


if __name__ == "__main__":
    a = Base(modFiles='ispymodule.py')
    a.createFolders(suffixNo=9)
    #a.replacePattern
