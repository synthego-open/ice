import numpy as np
import pandas as pd

def needleman_wunsch_matrix_pointers(score_matrix):
    """
    fill in the DP matrix according to the Needleman-Wunsch algorithm.
    Returns the matrix of scores and the matrix of pointers
    """
    ptr_dict = {'DIAG': 0, 'UP': 1, 'LEFT': 2}

    #     match =  1  # match score
    #     mismatch = -1  # mismatch penalty
    #     indel = -1 # indel penalty
    n, m = score_matrix.shape

    s = np.zeros((n + 1, m + 1))  # DP matrix
    ptr = np.zeros((n + 1, m + 1), dtype=int)  # matrix of pointers

    ##### INITIALIZE SCORING MATRIX (base case) #####
    s[:n, :m] = score_matrix

    #     for i in range(1, n+1) :
    #         s[i,0] = indel * i
    #     for j in range(1, m+1):
    #         s[0,j] = indel * j

    ########## INITIALIZE TRACEBACK MATRIX ##########

    # Tag first row by LEFT, indicating initial "-"s
    ptr[0, 1:] = ptr_dict['LEFT']

    # Tag first column by UP, indicating initial "-"s
    ptr[1:, 0] = ptr_dict['UP']

    #####################################################

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # match
            # relative scores
            relative_scores = {
                'DIAG': s[i - 1, j - 1],
                'UP': s[i - 1, j],
                'LEFT': s[i, j - 1]}

            optimal = max(relative_scores.keys(), key=(lambda key: relative_scores[key]))
            ptr[i, j] = ptr_dict[optimal]

            s[i, j] += relative_scores[optimal]

    return s, ptr


def needleman_wunsch_trace(score_matrix, ptr, ind,control_calls):
    #### TRACE BEST PATH TO GET ALIGNMENT ####
    ptr_dict = {'DIAG': 0, 'UP': 1, 'LEFT': 2}

    n, m = score_matrix.shape
    i = ind[0]
    j = ind[1]
    curr = ptr[i, j]

    x_array = []
    y_array = []
    aligned_control=''

    while (i > 0 or j > 0):

        if curr == ptr_dict['DIAG']:
            aligned_control += control_calls[i - 1]
            x_array.append(i - 1)
            y_array.append(j - 1)
            i -= 1
            j -= 1
        elif curr == ptr_dict['LEFT']:
            aligned_control += '*'
            x_array.append(i)
            y_array.append(j - 1)
            j -= 1
        elif curr == ptr_dict['UP']:
            aligned_control += '-'
            x_array.append(i - 1)
            y_array.append(j)
            i -= 1

        curr = ptr[i, j]

    return pd.DataFrame({'x': np.asarray(x_array), 'y': np.asarray(y_array)}),aligned_control[::-1]


def return_tracebacks(kernal_stacks,control_calls):

    tracebacks = []
    aligned_controls=[]

    for clip in np.linspace(0.4, 0.9, 20):
        clipped = (np.median(np.asarray(kernal_stacks), axis=0) > clip).astype('int')
        clipped[clipped == 0] = -1

        s, ptr = needleman_wunsch_matrix_pointers(clipped)

        ind = np.unravel_index(np.argmax(s, axis=None), s.shape)
        position_dict,aligned_control = needleman_wunsch_trace(s, ptr, ind,control_calls)
        tracebacks.append(pd.DataFrame(position_dict))
        aligned_controls.append(aligned_control)
    return tracebacks,aligned_controls
