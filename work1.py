"""
Liav Huli 314917808
Yehuda Harush 324023332
Tamir Refael 
Sagi Lidani 211451091
"""
def get_minor(matrix, i, j):
    """
    Get the minor of a matrix by removing the ith row and jth column.

    Parameters:
    matrix (list of list of floats): The original matrix.
    i (int): Row index to exclude.
    j (int): Column index to exclude.

    Returns:
    list of list of floats: Minor of the matrix.
    """
    return [row[:j] + row[j+1:] for row in (matrix[:i] + matrix[i+1:])]

def determinant(matrix):
    """
    Calculate the determinant of a square matrix recursively.

    Parameters:
    matrix (list of list of floats): The matrix.

    Returns:
    float: Determinant of the matrix.
    """
    size = len(matrix)
    if size == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    
    det = 0
    for c in range(size):
        det += ((-1) ** c) * matrix[0][c] * determinant(get_minor(matrix, 0, c))
    return det

def multiply_matrix(matrix1, matrix2):
    """
    Perform matrix multiplication of two matrices.

    Parameters:
    matrix1 (list of lists): First matrix.
    matrix2 (list of lists): Second matrix.

    Returns:
    list of lists: Result of matrix multiplication.
    """
    rows1 = len(matrix1)
    cols1 = len(matrix1[0])
    rows2 = len(matrix2)
    cols2 = len(matrix2[0])

    if cols1 != rows2:
        raise ValueError("Matrix dimensions are not compatible for multiplication.")

    result = [[0 for _ in range(cols2)] for _ in range(rows1)]

    for i in range(rows1):
        for j in range(cols2):
            for k in range(cols1):
                result[i][j] += matrix1[i][k] * matrix2[k][j]

    return result

def find_inverse_matrix(matrix):
    """
    Find the inverse of a given 3x3 matrix using elementary row matrices.
    
    Parameters:
    matrix (list of list of floats): A 3x3 matrix represented as a list of lists.
    
    Returns:
    list of list of floats: The inverse of the input matrix if it exists.
    """
    # Initialize the inverse matrix as the identity matrix
    inverse_matrix = [[1 if i == j else 0 for j in range(3)] for i in range(3)]
    
    # Apply row operations to transform the input matrix to the identity matrix
    for i in range(3):
        if matrix[i][i] == 0:
            raise ValueError("Matrix is singular, cannot invert.")
        
        # Scale the ith row to have 1 on diagonal
        scale_factor = 1 / matrix[i][i]
        for j in range(3):
            matrix[i][j] *= scale_factor
            inverse_matrix[i][j] *= scale_factor
        
        # Eliminate non-diagonal elements
        for k in range(3):
            if i != k:
                factor = matrix[k][i]
                for j in range(3):
                    matrix[k][j] -= factor * matrix[i][j]
                    inverse_matrix[k][j] -= factor * inverse_matrix[i][j]
    
    return inverse_matrix


def inverse(matrix):
    """
    Find the inverse of a matrix and return it.

    Parameters:
    matrix (list of list of floats): The matrix to invert.

    Returns:
    list of list of floats: The inverse of the matrix if it exists.
    
    Raises:
    ValueError: If the matrix is singular (determinant is zero).
    """
    # Check that the matrix is not singular
    det = determinant(matrix)
    if det == 0:
        raise ValueError("Matrix is singular, cannot invert.")
    
    # Calculate and return the inverse matrix
    return find_inverse_matrix(matrix)

def matrix_1_norm(matrix):
    """
    Calculate the 1-norm of a given matrix.

    Parameters:
    matrix (list of lists): The matrix to calculate the 1-norm for.

    Returns:
    float: The 1-norm of the matrix.
    """
    num_cols = len(matrix[0])
    col_sums = [0] * num_cols
    
    for row in matrix:
        for col in range(num_cols):
            col_sums[col] += abs(row[col])
    
    return max(col_sums)

# Example usage:
matrix = [
    [1, -1, -2],
    [2, -3, -5],
    [-1, 3, 5]
]

try:
    # Calculate inverse matrix
    inv_matrix = inverse(matrix)
    print("Inverse matrix:")
    for row in inv_matrix:
        print(row)
    
    # Calculate condition number
    cond = matrix_1_norm(matrix) * matrix_1_norm(inv_matrix)
    print("Condition number:", cond)

except ValueError as e:
    print(e)
