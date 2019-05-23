import sys, os


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
    # outbgl = '/Users/wansun/Git_Projects/CookHLA/tests/modules/BGL2SortBGl/Exon234_Panel_test'

    BGL2SortBGL(sort_markers, inbgl, outbgl)