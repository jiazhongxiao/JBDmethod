# IRJBD Documentation

IRJBD is an algorithm designed to compute the extreme generalized singular values and corresponding generalized singular vectors of matrix pairs.

## How to Use

1. **Prepare the Directory Structure**
   - Create a folder named `IRJBDmatrices` at the same level as the `IRJBD-Master` directory. This means that both folders (`IRJBDmatrices` and `IRJBD-Master`) should reside within the same parent directory.
   - Place your matrix files into the `IRJBDmatrices` folder. Each matrix file should be named as `xxx.mat` and contain a struct `Problem`.

2. **Matrix File Requirements**
   - The `Problem` struct should include the matrices `A` and `L`. These matrices can be loaded in MATLAB using the following commands:
     ```matlab
     S = load("../IRJBDmatrices/xxx.mat");
     A = S.Problem.A;
     L = S.Problem.L;
     ```
   - If the matrix file does not contain the matrix `L`, the program will automatically generate `L` as described in the paper.

## Usage Example

1. Download the matrix file `cat_ears_4_4.mat` from the [SuiteSparse Matrix Collection](https://suitesparse-collection-website.herokuapp.com/mat/JGD_Margulies/cat_ears_4_4.mat).
2. Place the file into the `IRJBDmatrices` folder.
3. Modify the parameters in `main.mat` as follows:
   ```matlab
   write = 1;
   targets = [5, -5];
   ks = [25, 50];
   ```
4. Run `main.mat`. The program will execute 4 times, using a maximum subspace dimension of 25 and 50 to calculate the 5 largest and smallest GSVD components of `{A, L}`, respectively. Here:
   - `A` is loaded from `cat_ears_4_4.mat`.
   - `L` is generated as described in the paper.

5. The results will be saved to a file named `results.csv`. The results should match those provided in the paper.
