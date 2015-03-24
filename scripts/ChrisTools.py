#!/usr/bin/python

import sys

def parse_fa(file):
    x = []
    id = ''
    s = ''

    with open(file) as f:
        while True:
            l = f.readline()
            if not l: break

            if l[0] == '>':
                if id != '':
                    x.append({ 'id': id.strip(), 'seq': s.strip() })
                
                # Start a new sequence
                id = l

                # Clear the old sequence
                s = ''
            else:
                s += l

        x.append({ 'id': id.strip(), 'seq': s.strip() })

    return x

def parse_bed(file):
    x = []
    with open(file) as f:
        while True:
            l = f.readline()
            if not l: break

            tokens = l.split('\t')
            if len(tokens) == 1:
                continue

            x.append({
                       "start": int(tokens[1]),
                       "end"  : int(tokens[2]),
                       "name" : tokens[3]
                     })

    return x

def parse_csv(file):
    x = []
    with open(file) as f:
        while True:
            l = f.readline()
            if not l: break

            tokens = l.split('\t')
            x.append({ "src": tokens[0].strip(), "dst": tokens[1].strip() })

    return x

def replace(x, b, y, output):
    # Assume x is a chromosome file...
    o = x[0]['seq']

    for i in range(len(b)):
        name = b[i]["name"]
        s = ''
        print("----> 1: Search for " + name)

        # Get the sequence to replace from y
        for j in range(len(y)):
            if y[j]['id'][1:] == name:
                s = y[j]['seq']
                print("----> 2: Found: " + s + " to replace, pos is: " + str(b[i]["start"]))
                break

        if s == '':
            print("Unable to find a matching sequence to replace for " + name)
            return

        # The position to replace
        p = b[i]["start"]

        t1 = o[:p]
        t2 = s
        t3 = o[p+len(t2):]

        print("----> 3: Src: " + o)
        print("----> 4: Dst: " + t1 + "___" + t2 + "___" + t3)
        print("")

        o = t1 + t2 + t3

    f = open(output, 'w')
    f.write(x[0]['id']+'\n')
    f.write(o)

def subtitute(x, y, output):
    # Assume x is a chromosome file...
    o = x[0]['seq']

    for i in range(len(y)):
        src = y[i]['src']
        dst = y[i]['dst']

        print("----> 1: " + src )
        print("----> 2: " + dst )

        o = o.replace(src, dst)

    f = open(output, 'w')
    f.write(x[0]['id']+'\n')
    f.write(o)

if __name__ == '__main__':

    print('Argument List:', str(sys.argv))

    if (sys.argv[1] == 'replace'):
        x = parse_fa(sys.argv[2])
        y = parse_bed(sys.argv[3])
        z = parse_fa(sys.argv[4])
        replace(x, y, z, sys.argv[5])
    elif (sys.argv[1] == 'subtitute'):
        #x = parse_fa("C://Sources//QA//Scripts//sub_tab.fa")
        #y = parse_csv("C://Sources//QA//Scripts//col.csv")
        #subtitute(x, y, "C://Sources//QA//Scripts//output_sub.txt")
        x = parse_fa(sys.argv[2])
        y = parse_csv(sys.argv[3])
        subtitute(x, y, sys.argv[4])
