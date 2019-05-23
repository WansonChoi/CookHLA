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


    selectMarkers = []
    first_time = True
    conut_num=1

    ### Obtaining marker labels from sorted *.markers file(`sort_markers`)
    with open(sort_markers) as fin:
        selectMarkers = [p.match(string=l).group(1) for l in fin]


    ### Sorting the order of `inbgl` file's markers based on that of sorted markers(`sort_markers` -> `selectMarkers`)
    with open(outbgl, 'w') as of:

        for t in selectMarkers:
            with open(inbgl) as pf:
                for l in pf:

                    # c = l.split()
                    # header = c[:2]
                    # data = c[2:]

                    m = p.match(string=l)
                    header = [m.group(1), m.group(2)]

                    if header[0] != "M" and first_time:
                        of.write(l)
                        # print(header)
                    if t == header[1]:
                        of.write(l)
                        # print(header)
                        conut_num = conut_num+1
                        break

                first_time = False



    return outbgl


# def GEN_MakeOutBGL(selectMarkers, inbgl, first_time, count_num):
#
#     for t in selectMarkers:
#         with open(inbgl) as pf:
#             for l in pf:
#
#                 # c = l.split()
#                 # header = c[:2]
#                 # data = c[2:]
#
#                 m = p.match(string=l)
#                 header = [m.group(1), m.group(2)]
#
#                 if header[0] != "M" and first_time:
#                     # newdata = []
#                     # for j in range(len(data)):
#                     #     newdata.append(data[j])
#                     of.write(l)
#                     print(header)
#                 if t == header[1]:
#                     # newdata = []
#                     # for j in range(len(data)):
#                     #     newdata.append(data[j])
#                     of.write(l)
#                     print(header)
#                     conut_num = conut_num + 1
#                     break
#
#             first_time = False




if __name__ == '__main__':

    """
    [BGL2SortBGL.py]
    
    
    """

    [sort_markers, inbgl, outbgl] = sys.argv[1:4]

    # ### <Test> ###
    #
    # # (2019. 05. 22.) First Testing (due to too much slow result)
    # sort_markers = '/Users/wansun/Git_Projects/CookHLA/tests/modules/BGL2SortBGl/input/Exon234_Panel_test.markers'
    # inbgl = '/Users/wansun/Git_Projects/CookHLA/tests/modules/BGL2SortBGl/input/T1DGC_REF.STEP2_exon234.bgl.phased'
    # outbgl = '/Users/wansun/Git_Projects/CookHLA/tests/modules/BGL2SortBGl/Exon234_Panel_test.sorted.bgl.phased'

    start = time.time()

    BGL2SortBGL_WS(sort_markers, inbgl, outbgl)

    end = time.time()

    # print("Time : %0.2f"%((end - start)/60))