//#include "Gaus.cpp"
//#include "Matrices_template.cpp"
#include "QR decomposition.cpp"
#include "other.cpp"


int main()
{
    //matrix<T> A = matrix<T>(3, 3, { 1,0,4,0,1,1,1,2,1 });
//matrix<T> b = matrix<T>(3, 1, { 1, 2, 3 });
//QR_decomposition(A, b);

//std::ifstream file("..//..//labs3kurs//matrix.dat");
//matrix<T> A, b;
//try
//{
//    // Читаем первую матрицу (A)
//    A.readFromFile(file);

//    file.close(); // Закрываем файл после чтения

//}
//catch (const std::exception& e) {
//    //std::cerr << "Error: " << e.what() << std::endl;
//    std::cerr << "File opening error" << std::endl;
//    if (file.is_open())
//    {
//        file.close();
//    }
//    return 1;
//}
    
    //matrix_test<T>("..\\labs3kurs\\tests\\test1_matrix.dat", "..\\labs3kurs\\tests\\test1_vector1.dat");
	/*matrix<T> A = matrix<T>("..\\labs3kurs\\tests\\test1_matrix.dat");
	matrix<T> b = matrix<T>("..\\labs3kurs\\tests\\test1_vector.dat");
	std::cout << "Matrix A: \n" << A << '\n' << "Vector b: \n" << b << std::endl;
	auto x1 = Gauss(A, b);
	auto x2 = QR_decomposition(A, b);

	std::cout << "Gauss ans: \n" << x1 << '\n' << "QR ans: \n" << x2 << std::endl;
	auto b1 = A * x1;
	auto residual1 = mismatch(b1, b);
	std::cout << "Ans residual = " << residual1 << '\n' << std::endl;
    
	auto condV = cond(A);
    std::cout << "Matrix cond = " << condV << std::endl;

	auto b0 = matrix<T>("..\\labs3kurs\\tests\\test1_vector1.dat");
	auto x01 = Gauss(A, b0);

	auto ans_diff1 = mismatch(x01, x1);

	auto Dx = 0.01;
	auto dx1 = Dx / norm1(x1);
	auto dxinf = Dx / norminf(x1);

	auto Db = 0.01;
	auto db1 = Db / norm1(b0);
	auto dbinf = Db / norminf(x1);

	auto cond_est1 = dx1 / db1;
	auto cond_estinf = dxinf / dbinf;

	std::cout << "Matrix cond estimation with norm1 = " << cond_est1;
	std::cout << "\nMatrix cond estimation with norminf = " << cond_estinf << std::endl;

	auto A1 = reverse(A);
	auto E = A1 * A;
	std::cout << "\nE: \n" << E << std::endl;

	auto btest = b_mod(b);
	std::cout << btest;*/
	using T = double;
	matrix_test<T>("..\\labs3kurs\\tests\\test3_matrix5.dat", "..\\labs3kurs\\tests\\test3_vector5.dat");
    
	return 0;
}