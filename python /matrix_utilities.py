from collections import namedtuple
import sys

def string_to_matrix(line):
    nums = line.split()
    row_count = int(nums[0])
    col_count = int(nums[1])
    elems = [float(n) for n in nums[2:]]
    rows = []
    for row_index in range(0,row_count):
        start = row_index*col_count
        end = start + col_count
        rows.append(elems[start:end])
    return rows

def matrix_to_string(matrix):
    if len(matrix) == 0:
        return "0 0"
    dims_info = str(len(matrix)) + " " + str(len(matrix[0]))

    row_to_str = lambda r: " ".join(map(str, r))
    return dims_info + " " + " ".join(map(row_to_str, matrix))

def transp(mat):
    return list(map(list, zip(*mat)))

def matmul(mat1, mat2):
    mat2_transp = transp(mat2)
    mat = []
    for row in mat1:
        new_row = []
        for col in mat2_transp:
            new_row.append(sum(map(lambda t: t[0]*t[1], zip(row,col))))
        mat.append(new_row)
    return mat

def times_vectors(vec1, vec2):
    return [a*b for (a,b) in zip(vec1, vec2)]

def times(mat1, mat2):
    zipped_rows = zip(mat1, mat2)
    return list(map(times_vectors, zipped_rows))

def getcol(mat, i):
    return transp(mat)[i]

def dims(matrix):
    if(len(matrix) == 0):
        return (0,0)
    return (len(matrix), len(matrix[0]))

def pretty_print(mat):
    for row in mat:
        str_row = map('{0:.5f}'.format, row)
        print(" ".join(str_row))

Model = namedtuple("Model", "a, b, pi, n")
def parse_model():
    A = string_to_matrix(sys.stdin.readline().strip())
    B = string_to_matrix(sys.stdin.readline().strip())
    pi = string_to_matrix(sys.stdin.readline().strip())[0]
    return Model(a=A, b=B, pi=pi, n=len(A))
