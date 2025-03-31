#include "Main.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

#pragma region General 

double Random(double min, double max)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

Matrix::Matrix(int row, int col)
{
    this->row = row;
    this->col = col;
    this->matrix = new double* [row];

    for (int i = 0;i < row;i++) {
        matrix[i] = new double[col];
    }
}

//Copying constructor
Matrix::Matrix(const Matrix& M)
{
    this->row = M.row;
    this->col = M.col;
    this->matrix = new double* [M.row];

    for (int i = 0;i < row;i++) {
        matrix[i] = new double[col];
    }
    for (int i = 0;i < row;i++) {
        for (int j = 0;j < col;j++) {
            this->matrix[i][j] = M.matrix[i][j];
        }
    }
}

Matrix::~Matrix()
{
    if (this->matrix != nullptr) {
        for (int i = 0;i < row;i++) {
            delete matrix[i];
        }
        delete[]matrix;
    }
}

void Matrix::Fill()
{
    for (int i = 0;i < row;i++) {
        for (int j = 0;j < col;j++) {
            std::cout << "Input [" << i + 1 << "][" << j + 1 << "]: ";
            std::cin >> matrix[i][j];
            std::cout << std::endl;
        }
    }
}

void Matrix::FillRand()
{
    for (int i = 0;i < row;i++) {
        for (int j = 0;j < col;j++) {
            //matrix[i][j] = Random(-5, 5);
            matrix[i][j] = rand() % 10;
        }
    }
}

void Matrix::Print()
{
    for (int i = 0;i < row;i++) {
        for (int j = 0;j < col;j++) {
            std::cout << std::setw(12) << std::left << std::setprecision(4) << matrix[i][j];
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

#pragma endregion 

#pragma region Operators 

Matrix& Matrix::operator = (const Matrix& A)
{
    if (this == &A) {
        return *this;
    }
    else {
        for (int i = 0;i < A.row;i++) {
            for (int j = 0;j < A.col;j++) {
                matrix[i][j] = A.matrix[i][j];
            }
        }
        return *this;
    }
}

Matrix operator +(const Matrix& A, const Matrix& B)
{
    if (A.row != B.row || A.col != B.col) {
        throw std::invalid_argument("Matrices must have the same dimensions for addition");
    }
    Matrix C(A.row, A.col);
    for (int i = 0;i < A.row;i++) {
        for (int j = 0;j < A.col;j++) {
            C.matrix[i][j] = A.matrix[i][j] + B.matrix[i][j];
        }
    }
    return C;
}

Matrix operator -(const Matrix& A, const Matrix& B)
{
    if (A.row != B.row || A.col != B.col) {
        throw std::invalid_argument("Matrices must have the same dimensions for addition");
    }
    Matrix C(A.row, A.col);
    for (int i = 0;i < A.row;i++) {
        for (int j = 0;j < A.col;j++) {
            C.matrix[i][j] = A.matrix[i][j] - B.matrix[i][j];
        }
    }
    return C;
}

Matrix operator *(const Matrix& A, const Matrix& B)
{
    if (A.col != B.row) {
        throw std::invalid_argument("Matrices must be compatible");
    }
    Matrix C(A.row, B.col);

    for (int i = 0;i < A.row;i++) {
        for (int j = 0; j < B.col;j++) {
            C.matrix[i][j] = 0;
            for (int k = 0;k < B.row;k++) {
                C.matrix[i][j] += A.matrix[i][k] * B.matrix[k][j];
            }
        }
    }

    return C;
}

Matrix operator *(const double scalar, const Matrix& B)
{
    Matrix C(B.row, B.col);
    for (int i = 0;i < B.row;i++) {
        for (int j = 0;j < B.col;j++) {
            C.matrix[i][j] = scalar * B.matrix[i][j];
        }
    }
    return C;
}
#pragma endregion

#pragma region Functions

Matrix Matrix::Transpose() const {
    Matrix M(col, row);

    for (int i = 0;i < col; i++) {
        for (int j = 0;j < row;j++) {
            M.matrix[i][j] = matrix[j][i];
        }
    }
    return M;
}

double Matrix::Determinant(int col) const {
    if (row != col) {
        throw std::invalid_argument("Matrix must be square");
    }
    if (col == 1) {
        return matrix[0][0];
    }
    if (col == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }
    double det = 0.0;
    for (int p = 0;p < col;p++) {
        Matrix SubM(col - 1, col - 1);
        for (int i = 1;i < col;i++) {
            int subj = 0;
            for (int j = 0;j < col;j++) {
                if (j == p) continue;
                SubM.matrix[i - 1][subj++] = matrix[i][j];
            }
        }
        det += matrix[0][p] * ((p % 2 == 0 ? 1 : -1) * SubM.Determinant(col - 1));
    }
    return det;
}

double Matrix::Trace() const {
    if (row != col) {
        throw std::invalid_argument("Trace is not defined , matrix should be square");
    }
    double trace = 0.0;
    for (int i = 0;i < col;i++) {
        trace += matrix[i][i];
    }
    return trace;
}

double Matrix::Minor(int col, int indrow, int indcol) const {
    double minor = 0.0;
    Matrix Temp(col - 1, col - 1);
    int subi = 0;
    for (int i = 0;i < col;i++) {
        if (i == indrow) continue;
        int subj = 0;
        for (int j = 0;j < col;j++) {
            if (j == indcol) continue;
            Temp.matrix[subi][subj] = matrix[i][j];
            subj++;
        }
        subi++;
    }
    minor = Temp.Determinant(col - 1);
    return minor;
}

Matrix Matrix::Inverse(const Matrix& A) {
    double det = A.Determinant(col);
    if (row != col || det == 0) {
        throw std::invalid_argument("The inverse matrix is undefined");
    }
    Matrix InvMatrix(col, col);

    for (int i = 0;i < col;i++) {
        for (int j = 0;j < col;j++) {
            InvMatrix.matrix[j][i] = (1.0 / det) * (pow(-1, i + j)) * A.Minor(col, i, j);
        }
    }
    return InvMatrix;
}

double Matrix::DotProduct(const Matrix& V1, const Matrix& V2)
{
    if (V1.row != V2.row) {
        throw std::invalid_argument("Dot product is not defined");
    }
    double dotresult = 0.0;
    for (int i = 0;i < V1.row;i++) {
        dotresult += V1.matrix[i][0] * V2.matrix[i][0];
    }
    return dotresult;
}

double Matrix::Norm()const {
    double norm = 0.0;
    double sum = 0.0;
    for (int i = 0;i < row;i++) {
        for (int j = 0; j < col;j++) {
            sum += matrix[i][j] * matrix[i][j];
        }
    }
    norm = sqrt(sum);
    return norm;
}

Matrix Matrix::NormalizeVector(const Matrix& V) {
    Matrix Vn(V.row, 1);
    Vn = V;
    double vectNorm = Vn.Norm();
    if (vectNorm != 0) {
        for (int i = 0;i < Vn.row;i++) {
            Vn.matrix[i][0] = (Vn.matrix[i][0] / vectNorm);
        }
    }
    return Vn;
}

Matrix Matrix::SplitVector(int w)const {
    Matrix Vector(row, 1);
    for (int i = 0;i < row;i++) {
        Vector.matrix[i][0] = matrix[i][w];
    }
    return Vector;
}

Matrix Matrix::SplitRow(int r)const {
    Matrix SplRow(1, col);

    for (int i = 0;i < col;i++) {
        SplRow.matrix[0][i] = matrix[r][i];
    }
    return SplRow;
}

Matrix Matrix::ZeroVector(int z)const {
    Matrix Zv(z, 1);

    for (int i = 0;i < z;i++) {
        Zv.matrix[i][0] = 0;
    }
    return Zv;
}

Matrix Matrix::IdentityMatrix(int s) const {
    Matrix I(s, s);

    for (int i = 0;i < s;i++) {
        for (int j = 0;j < s; j++) {
            (i == j) ? I.matrix[i][i] = 1 : I.matrix[i][j] = 0;
        }
    }
    return I;
}

bool   Matrix::IsSpecial(const Matrix& A) {
    // This is a special kind of matrices that has only the non-zero antidiagonal 
    int counter = 0;
    double sum = 0.0;
    for (int i = 0;i < A.col;i++) {
        for (int j = 0;j < A.col;j++) {
            if ((i + j) == A.col - 1) {
                sum += A.matrix[i][j];
            }
            if ((i + j) != A.col - 1 && A.matrix[i][j] == 0) {
                counter++;
            }
        }
    }

    if (counter == (A.row * A.col - A.col) && sum != 0) {

        return true;
    }
    return false;
}
#pragma endregion 

#pragma region SLE
Matrix Matrix::REF(Matrix& RA) {
    Matrix A(RA.row, RA.col);
    A = RA;
    Matrix ref(A.row, A.col);

    double t;
    int flag = std::min(A.row, A.col);

    for (int i = 0;i < flag;i++) {
        int subj = 0;
        if (abs(A.matrix[i][i]) < 1e-9) {
            A.matrix[i][i] = 0;
            // Looking for non-zero pivot (row)
            for (int k = i + 1;k < A.row;k++) {
                if (abs(A.matrix[k][i]) < 1e-9) {
                    //swap row i and row k
                    for (int l = 0;l < A.col;l++) {
                        t = A.matrix[i][l];
                        A.matrix[i][l] = A.matrix[k][l];
                        A.matrix[k][l] = t;
                    }
                    break;
                }
            }
            if (abs(A.matrix[i][i]) < 1e-9) continue;
        }

        for (int k = i + 1; k < A.row;k++) {
            double coef1 = A.matrix[k][i] / A.matrix[i][i];
            for (int j = i;j < A.col;j++) {
                A.matrix[k][j] -= coef1 * A.matrix[i][j];
            }
        }

    }

    for (int i = flag - 1;i > -1;i--) {
        for (int k = i - 1; k > -1;k--) {
            if (abs(A.matrix[i][i]) < 1e-9) continue;
            double coef2 = A.matrix[k][i] / A.matrix[i][i];
            for (int j = i;j > -1;j--) {
                A.matrix[k][j] -= coef2 * A.matrix[i][j];
            }
        }
    }

    for (int i = 0;i < A.row;i++) {
        for (int j = 0;j < A.col;j++) {
            if (abs(A.matrix[i][j]) < 1e-9) A.matrix[i][j] = 0;
            ref.matrix[i][j] = A.matrix[i][j];
        }
    }
    return ref;
}

int    Matrix::Rank(Matrix& A) {
    Matrix T(A.row, A.col);
    T = REF(A);
    int rank = T.row;
    int counter = 0;
    //Checking
    for (int i = T.row - 1;i > -1;i--) {
        for (int j = T.col - 1; j > -1;j--) {
            if (abs(T.matrix[i][j]) < 1e-9) { counter++; }
            if (abs(T.matrix[i][j]) > 1e-9) { break; }
        }
        if (counter == col) rank -= 1;
        counter = 0;
    }
    return rank;
}

Matrix Matrix::HSLE(Matrix& A) {
    Matrix Solution(A.col, 1);

    Matrix N(A.row, A.col);

    N = REF(A);

    // Checking artificial rank
    int rk = N.row;
    int count = 0;
    for (int i = N.row - 1;i > -1;i--) {
        for (int j = N.col - 1; j > -1;j--) {
            if (abs(N.matrix[i][j]) < 1e-9) { count++; }
            if (abs(N.matrix[i][j]) > 1e-9) { break; }
        }
        if (count == N.col) rk -= 1;
        count = 0;
    }

    if ((N.row > N.col) && rk == N.row) {
        throw std::invalid_argument("There is no solutions");
    }
    if (rk == N.col) {
        for (int j = 0;j < N.col;j++) {
            Solution.matrix[j][0] = 0;
        }
        return Solution;
    }
    for (int i = 0;i < N.row;i++) {
        double sum = 0.0;
        for (int j = 0;j < N.col;j++) {
            if (N.matrix[i][j] != 0) {
                for (int k = j + 1;k < N.col;k++) {
                    sum += N.matrix[i][k] / N.matrix[i][j];
                }
                (sum != 0) ? Solution.matrix[i][0] = -sum : Solution.matrix[i][0] = 0;
                break;
            }
        }
    }
    if (rk < N.col) {
        for (int i = rk;i < N.col;i++) {
            Solution.matrix[i][0] = 1;
        }
    }
    return Solution;
}

Matrix Matrix::NHSLE(Matrix& A, Matrix& b) {
    if (A.row != b.row || b.col != 1) {
        throw std::invalid_argument("Dimensions are incorrect");
    }

    Matrix Ab(A.row, A.col + 1);

    for (int i = 0;i < A.row;i++) {
        for (int j = 0;j < A.col + 1;j++) {
            if (j == A.col) {
                Ab.matrix[i][j] = b.matrix[i][0];
            }
            else {
                Ab.matrix[i][j] = A.matrix[i][j];
            }

        }
    }
    Ab = Ab.REF(Ab);

    int rkA = A.Rank(A);
    int rkAb = Ab.Rank(Ab);

    if (A.row > A.col && rkAb == A.row) {
        throw std::invalid_argument("There is no solutions");
    }
    if (A.row > A.col && rkA == rkAb) {
        throw std::invalid_argument("There is linear depended rows ,please reduce them and try again ");
    }
    if (rkA < rkAb) {
        throw std::invalid_argument("There is no solutions");
    }
    if (rkA == rkAb && A.row == A.col) {
        Matrix Solution(A.row, 1);

        for (int i = 0; i < Ab.row; i++) {
            if (Ab.matrix[i][i] != 0) {
                //throw std::invalid_argument("Division by zero");
                (Ab.matrix[i][A.col] != 0) ? Solution.matrix[i][0] = -(Ab.matrix[i][A.col] / Ab.matrix[i][i]) : Solution.matrix[i][0] = 0;

            }
        }
        return Solution;
    }
    if (rkA == rkAb && A.row < A.col) {
        Matrix Solution(A.col, Ab.col - rkA);

        //Part Solution
        for (int i = 0;i < A.col;i++) {
            if (i < rkAb && Ab.matrix[i][i] != 0) {
                //throw std::invalid_argument("Division by zero");
                (Ab.matrix[0][Ab.col - 1] != 0) ? Solution.matrix[i][0] = -(Ab.matrix[0][Ab.col - 1] / Ab.matrix[i][i]) : Solution.matrix[i][0] = 0;
            }
            else {
                Solution.matrix[i][0] = 0;
            }
        }

        // Homogeneous SLE solution
        for (int j = 1;j < Ab.col - rkA;j++) {
            for (int i = 0;i < A.col;i++) {
                if (i < rkAb && Ab.matrix[i][i] != 0) {
                    //throw std::invalid_argument("Division by zero");
                    (Ab.matrix[i][rkA + j] != 0) ? Solution.matrix[i][j] = -(Ab.matrix[i][rkA + j] / Ab.matrix[i][i]) : Solution.matrix[i][j] = 0;
                }
                else {
                    (i == rkA + j - 1) ? Solution.matrix[i][j] = 1 : Solution.matrix[i][j] = 0;
                }
            }
        }

        return Solution;
    }
}

#pragma endregion

#pragma region QRdecomposition

Matrix Matrix::Qmatrix(const Matrix& A) {
    Matrix Q(A.row, A.col);

    for (int j = 0;j < A.col;j++) {
        Matrix Aj = A.SplitVector(j);
        for (int i = 0;i < A.row;i++) {
            Q.matrix[i][j] = Aj.matrix[i][0];
        }
        //Substraction
        for (int k = 0; k < j;k++) {
            Matrix Qk(A.row, 1);

            Qk = Q.SplitVector(k);
            double norm = DotProduct(Qk, Qk);
            for (int i = 0;i < A.row;i++) {
                if (norm != 0) {
                    Q.matrix[i][j] -= (DotProduct(Qk, Aj) / norm) * Q.matrix[i][k];
                }
            }
        }
        for (int j = 0;j < col;j++) {
            Matrix Qj = Q.SplitVector(j);
            double normj = DotProduct(Qj, Qj);
            if (normj != 0) {
                for (int i = 0;i < row;i++) {
                    Q.matrix[i][j] = NormalizeVector(Qj).matrix[i][0];
                }
            }
        }
    }
    return Q;
}

Matrix Matrix::Rmatrix(const Matrix& Q, const Matrix& A) {
    Matrix R(A.row, A.row);

    for (int j = 0;j < col;j++) {
        Matrix aj = A.SplitVector(j);
        for (int i = 0;i < row;i++) {
            Matrix ei = Q.SplitVector(i);
            R.matrix[i][j] = DotProduct(ei, aj);
        }
    }
    return R;
}

Matrix Matrix::QRalgBasic(const Matrix& A)
{
    if (IsSpecial(A)) {
        throw std::invalid_argument("The matrix is special");
    }
    Matrix CurrentA(A.row, A.col);

    CurrentA = A; // A(0) = A

    Matrix Q(A.row, A.col);
    Matrix R(A.row, A.row);

    int counter = 0;

    while (counter < 1000) {
        Q = CurrentA.Qmatrix(CurrentA);
        R = CurrentA.Rmatrix(Q, CurrentA);
        // A(k) = Q(k)R(k)
        // A(k+1) = R(k)*Q(k)
        CurrentA = R * Q;
        counter++;
    }

    return CurrentA;
}

Matrix Matrix::HessenbergMatrix(const Matrix& A) {
    if (IsSpecial(A)) {
        throw std::invalid_argument("The matrix is special");
    }
    Matrix Hess(A.row, A.col);
    Hess = A;

    for (int i = 0;i < A.col - 2;i++) { // Two last columns we don't make to zero 
        Matrix a(A.row - i - 1, 1);

        for (int j = i + 1;j < A.row;j++) {
            a.matrix[j - i - 1][0] = Hess.matrix[j][i];
        }
        double sign = (a.matrix[0][0] < 0) ? -1.0 : 1.0;
        Matrix V = a + sign * (a.Norm()) * IdentityMatrix(A.row - i - 1).SplitVector(0);
        V = V.NormalizeVector(V);
        Matrix P = IdentityMatrix(A.row - i - 1) - 2 * V * (V.Transpose());

        Matrix P1(A.row, A.col); // Extended matrix

        for (int k = 0;k < A.row;k++) {
            for (int l = 0;l < A.col;l++) {
                if (k <= i || l <= i) {
                    P1.matrix[k][l] = (k == l) ? 1.0 : 0.0;
                }
                else {
                    P1.matrix[k][l] = P.matrix[k - i - 1][l - i - 1];
                }
            }
        }
        Hess = P1 * Hess * P1.Transpose();
    }
    return Hess;
}

Matrix Matrix::QRalgHess(const Matrix& A)
{
    if (IsSpecial(A)) {
        throw std::invalid_argument("The matrix is special");
    }
    Matrix I = IdentityMatrix(A.row);

    Matrix CurrentA(A.row, A.col);
    CurrentA = HessenbergMatrix(A);

    Matrix Q(A.row, A.col);
    Matrix R(A.row, A.row);

    int counter = 0;

    while (counter < 1000) {
        double shift = CurrentA.matrix[A.row - 1][A.col - 1];

        CurrentA = CurrentA - shift * I; // Ak - sk*I = Q*R

        Q = CurrentA.Qmatrix(CurrentA);
        R = CurrentA.Rmatrix(Q, CurrentA);

        CurrentA = R * Q + shift * I;// A(k+1) = R*Q +sk*I  

        counter++;
    }

    return CurrentA;
}
#pragma endregion

#pragma region SVD

Matrix Matrix::EigenValMatrix(const Matrix& A) {
    if (A.row != A.col) {
        throw std::invalid_argument("Non-square matrix !!!");
    }
    Matrix EV(A.row, A.col);

    Matrix Temp(A.row, A.col);
    Temp = QRalgHess(A);

    for (int i = 0;i < A.row;i++) {
        for (int j = 0;j < A.col;j++) {
            (i == j) ? EV.matrix[i][j] = Temp.matrix[i][j] : EV.matrix[i][j] = 0;
        }
    }
    for (int i = 0;i < A.row;i++) {
        if (abs(EV.matrix[i][i]) < 1e-9) {
            EV.matrix[i][i] = 0;
        }
    }
    return EV;
}

Matrix Matrix::EigenVectMatrix(const Matrix& A)
{
    if (A.row != A.col) {
        throw std::invalid_argument("Non-square matrix !!!");
    }
    Matrix P(A.col, A.col);

    Matrix EigVal(A.col, A.col);
    EigVal = EigenValMatrix(A);

    Matrix Ptemp(A.col, A.col);

    Matrix Temp(A.col, A.col);
    Temp = A;
    Matrix I = IdentityMatrix(A.col);
    Matrix pk(A.col, 1);

    for (int k = 0;k < A.col;k++) {
        Ptemp = Temp - EigVal.matrix[k][k] * I;
        //Ptemp.Print();
        pk = HSLE(Ptemp);
        for (int i = 0;i < A.col;i++) {
            P.matrix[i][k] = pk.matrix[i][0];
        }
    }
    return P;
    std::cout << "EigenVectMatrix" << std::endl;
}

Matrix Matrix::Smatrix(const Matrix& A) {
    Matrix E(A.col, A.col);
    E = (A.Transpose()) * A;

    Matrix D(A.col, A.col);
    D = EigenValMatrix(E);

    Matrix Sigma(A.row, A.col);

    int dimSing = std::min(A.row, A.col);
    Matrix T(dimSing, 1);
    for (int i = 0;i < dimSing;i++) {
        T.matrix[i][0] = 0;
    }
    int subi = 0;
    for (int i = 0; i < A.col; i++) {
        if (D.matrix[i][i] != 0) {
            T.matrix[subi][0] = D.matrix[i][i];
            subi++;
        }
    }
    //Bubble Sorting
    double t;
    for (int i = dimSing - 1;i > -1;i--) {
        for (int j = 0;j < i;j++) {
            if (T.matrix[j][0] < T.matrix[j + 1][0]) {
                t = T.matrix[j][0];
                T.matrix[j][0] = T.matrix[j + 1][0];
                T.matrix[j + 1][0] = t;
            }
        }
    }
    int subj = 0;
    for (int i = 0;i < A.row;i++) {
        for (int j = 0;j < A.col;j++) {
            if (i == j) {
                Sigma.matrix[i][j] = sqrt(T.matrix[subj][0]);
                subj++;
            }
            else {
                Sigma.matrix[i][j] = 0;
            }
        }
    }
    return Sigma;
}

Matrix Matrix::Vmatrix(const Matrix& A) {
    Matrix V(A.col, A.col);

    Matrix W(A.col, A.col);
    W = (A.Transpose()) * A;

    V = W.EigenVectMatrix(W);
    Matrix Vk(A.col, 1);

    for (int j = 0;j < A.col;j++) {
        Vk = V.SplitVector(j);
        Vk = Vk.NormalizeVector(Vk);

        for (int i = 0;i < A.col;i++) {
            V.matrix[i][j] = Vk.matrix[i][0];
        }
    }
    return V;

}

Matrix Matrix::Umatrix(const Matrix& A, const Matrix& V, const Matrix& S) {
    Matrix U(A.row, A.row);

    Matrix Uk(A.row, 1);

    for (int j = 0;j < A.row;j++) {
        if (S.matrix[j][j] != 0) Uk = (1 / S.matrix[j][j]) * (A * V.SplitVector(j));
        if (S.matrix[j][j] == 0) Uk = Uk.ZeroVector(A.row);
        for (int i = 0;i < A.row;i++) {
            U.matrix[i][j] = Uk.matrix[i][0];
        }
    }
    return U;
}

#pragma endregion