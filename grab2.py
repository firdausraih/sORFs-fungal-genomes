#grep cluster yg ada 3 3 result atau 2 result sofwares#
#usage:python grab2.py <query file> > <output.file>
#!/usr/bin/python
import sys
import os
import glob
import re

def main():
    _inputfile = open(sys.argv[1], "r")

    # flags
    flag = 0
    getigr_seq = []
    getge_seq = []
    sf_seq = []
    f_seq = []

    for read in _inputfile:

        if flag is 1:
            if re.search(r'>getigr', read):
                getigr_seq.append(read)

            elif re.search(r'>getge', read, re.M|re.I):
                getge_seq.append(read)

            elif re.search(r'>sf', read, re.M|re.I):
                sf_seq.append(read)

            elif re.search(r'>Cluster', read, re.M|re.I):
                flag = 0
                if (len(getigr_seq) > 0 and len(getge_seq) > 0 and len(sf_seq) > 0):
                    f_seq.extend(getigr_seq)
                    f_seq.extend(getge_seq)
                    f_seq.extend(sf_seq)
                    for i in f_seq:
                        print i.rstrip()

                getigr_seq = []
                getge_seq = []
                sf_seq = []
                f_seq = []

        if flag is 0:
            if re.search(r'>Cluster', read, re.M|re.I):
                flag = 1
                f_seq.append(read)

if __name__ == "__main__":
    main()


