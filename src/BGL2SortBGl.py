import sys, os, re
import time

p = re.compile(r'^(\w+)\s+(\w+)') # Frist two columns (ex. 'P pedigree' or 'rs969931 29602876', ... )


def BGL2SortBGL(sort_markers, inbgl, outbgl):


    selectMarkers = {}
    first_time = True
    conut_num=1


    with open(sort_markers) as fin, open(outbgl, 'w') as of:
        for l in fin:
            c = l.split()
            selectMarkers[c[0]] = True
        for t in selectMarkers:
            with open(inbgl) as pf:
                for l in pf:
                    c = l.split()
                    header = c[:2]
                    data = c[2:]
                    if header[0] != "M" and first_time:
                        newdata = []
                        for j in range(len(data)):
                            newdata.append(data[j])
                        of.write(' '.join(header + newdata) + '\n')
                    if t == header[1]:
                        newdata = []
                        for j in range(len(data)):
                            newdata.append(data[j])
                        conut_num = conut_num+1
                        of.write(' '.join(header + newdata) + '\n')
                        continue

                first_time = False


    return outbgl


def BGL2SortBGL_WS(sort_markers, inbgl, outbgl):



    ### Obtaining marker labels from sorted *.markers file(`sort_markers`)

    selectMarkers = []
    with open(sort_markers) as fin:
        selectMarkers = [p.match(string=l).group(1) for l in fin]


    ### Proprocessing input beagle file(`inbgl`) using enough memory
    IN_BGL_header = {}
    IN_BGL_body = {}

    with open(inbgl) as pf:
        for l in pf:

            m = p.match(string=l)

            if m.group(1) != 'M':
                IN_BGL_header[m.group(2)] = l
            else:
                IN_BGL_body[m.group(2)] = l


    # Header check.
    # count = 0
    #
    # for k, v in IN_BGL_header.items():
    #     print("{} : {}".format(k, v))
    #     count += 1
    #
    #     if count > 5: break;
    #
    # print("Memory : %0.2f(Mb)"%(sys.getsizeof(IN_BGL_header)/1024**2))

    # Body check.
    # count = 0
    #
    # for k, v in IN_BGL_body.items():
    #     print("{} : {}".format(k, v))
    #     count += 1
    #
    #     if count > 5: break;
    #
    # print("Memory : %0.2f(Mb)"%(sys.getsizeof(IN_BGL_body)/1024**2))


    ### Sorting the order of `inbgl` file's markers(rows) based on that of sorted markers(`sort_markers` -> `selectMarkers`)
    with open(outbgl, 'w') as of:

        # Writing Header
        # print("Writing Header")
        of.writelines((v for v in IN_BGL_header.values()))

        # Writing Body
        # print("Writing Body")
        of.writelines((IN_BGL_body[item] for item in selectMarkers))


    return outbgl





if __name__ == '__main__':

    """
    [BGL2SortBGL.py]
    
    The module to sort the rows of given Input Beagle file(*.bgl(.phsed)) based on the order of given Input Marker file
    (*.markers)
    
    """

    [sort_markers, inbgl, outbgl] = sys.argv[1:4]

    ### <Test> ###

    # # (2019. 05. 22.) First Testing (due to too much slow result)
    # sort_markers = '/Users/wansun/Git_Projects/CookHLA/tests/modules/BGL2SortBGl/input/Exon234_Panel_test.markers'
    # inbgl = '/Users/wansun/Git_Projects/CookHLA/tests/modules/BGL2SortBGl/input/T1DGC_REF.STEP2_exon234.bgl.phased'
    # outbgl = '/Users/wansun/Git_Projects/CookHLA/tests/modules/BGL2SortBGl/Exon234_Panel_test.sorted.generator2.bgl.phased'

    # start = time.time()

    BGL2SortBGL_WS(sort_markers, inbgl, outbgl)

    # end = time.time()

    # print("Time : %0.2f"%((end - start)/60))