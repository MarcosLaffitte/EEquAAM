#!/usr/bin/env python3

import sys

inputFileName = sys.argv[1]
outputFileName = sys.argv[2]

with open(inputFileName, mode='r') as inp:
    with open(outputFileName, mode='w') as out:
        count = 0
        for i in inp:
            print(i)
            print(count)
            if count == 4:
                out.write(i)
                out.write('\n')
                count = 0
            else: 
                out.write(i)
                count = count +1
                
