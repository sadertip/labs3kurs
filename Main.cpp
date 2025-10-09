//#include "Gaus.cpp"
//#include "Matrices_template.cpp"
#include "simple_iteration_method.cpp"
#include "QR decomposition.cpp"
#include "other.cpp"


int main()
{
	using T = double;
	matrix<T> A1 = matrix<T>(3, 3, { 1,0,4,0,1,1,1,2,1 });
	matrix<T> b1 = matrix<T>(3, 1, { 1, 2, 3 });
	matrix<T> x0(3, 1, { 0.0, 0.0,0.0 });
	T tau = 0.05;

	T omega = 1.2;
	int max_iter = 100;
	T epsilon = 1e-4;

	/*matrix<T> solution = sor_method(A, b, x0, 1.0, max_iter, epsilon);

	std::cout << "\nFinal Solution:\n" << solution << "\n";

	solution = simple_iteration_method(A, b, x0, tau, max_iter, epsilon);

	std::cout << "\nFinal Solution:\n" << solution << "\n";*/

	try {

		matrix<T> A = matrix<T>(4, 4, { 15,2,-3,7,-5,11,2,-3,0,-1,7,4,12,0,-6,20 });
		matrix<T> b = matrix<T>(4, 1, { 53, -90, 107,68 });
		matrix<T> x0(4, 1, { 1.0, 1.0,1.0,1.0 });

		std::cout << "Matrix A:\n" << A << "\n";
		std::cout << "Vector b:\n" << b << "\n";

		matrix<T> solution1 = simple_iteration_method(A, b, x0, tau, max_iter, epsilon);
		matrix<T> solution2 = jacobi_method_elementwise(A, b, x0, max_iter, epsilon);
		matrix<T> solution3 = sor_method(A, b, x0, 1., max_iter, epsilon);
		matrix<T> solution4 = sor_method(A, b, x0, omega, max_iter, epsilon);

		std::cout << "\nsimple_iteration Solution:\n" << solution1 << "\n";
		std::cout << "\njacobi_method Solution:\n" << solution2 << "\n";
		std::cout << "\nSeidel Solution:\n" << solution3 << "\n";
		std::cout << "\nRelaxation Solution:\n" << solution4 << "\n";

		std::pair<matrix<T>, matrix<T>> components1 = get_simple_iteration_components(A, b, tau);

		std::cout << "\nsimple_iteration C :\n";
		std::cout << components1.first << std::endl;
		std::cout << "\nsimple_iteration y :\n";
		std::cout << components1.second << std::endl;
		std::cout << "\nnorm C :\n";
		std::cout << norminf(components1.first) << std::endl;

		std::pair<matrix<T>, matrix<T>> components2 = get_jacobi_components(A, b);

		std::cout << "\njacobi C :\n";
		std::cout << components2.first << std::endl;
		std::cout << "\njacobi y :\n";
		std::cout << components2.second << std::endl;
		std::cout << "\nnorm C :\n";
		std::cout << norminf(components2.first) << std::endl;

		SorComponents<T> components = calculate_sor_matrices(A, omega);



		std::cout << "\nCL :\n";
		std::cout << components.CL << std::endl;

		std::cout << "\nCU :\n";
		std::cout << components.CU << std::endl;


		std::cout << "\nC :\n";
		std::cout << components.C << std::endl;


		calculate_and_sum_norms(components.CL, components.CU);



	}
	catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}
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



	return 0;
}
////#include "Gaus.cpp"
////#include "Matrices_template.cpp"
//#include "QR decomposition.cpp"
//#include "other.cpp"
//
//
//int main()
//{
//    //matrix<T> A = matrix<T>(3, 3, { 1,0,4,0,1,1,1,2,1 });
////matrix<T> b = matrix<T>(3, 1, { 1, 2, 3 });
////QR_decomposition(A, b);
//
////std::ifstream file("..//..//labs3kurs//matrix.dat");
////matrix<T> A, b;
////try
////{
////    // Читаем первую матрицу (A)
////    A.readFromFile(file);
//
////    file.close(); // Закрываем файл после чтения
//
////}
////catch (const std::exception& e) {
////    //std::cerr << "Error: " << e.what() << std::endl;
////    std::cerr << "File opening error" << std::endl;
////    if (file.is_open())
////    {
////        file.close();
////    }
////    return 1;
////}
//    
//    //matrix_test<T>("..\\labs3kurs\\tests\\test1_matrix.dat", "..\\labs3kurs\\tests\\test1_vector1.dat");
//	/*matrix<T> A = matrix<T>("..\\labs3kurs\\tests\\test1_matrix.dat");
//	matrix<T> b = matrix<T>("..\\labs3kurs\\tests\\test1_vector.dat");
//	std::cout << "Matrix A: \n" << A << '\n' << "Vector b: \n" << b << std::endl;
//	auto x1 = Gauss(A, b);
//	auto x2 = QR_decomposition(A, b);
//
//	std::cout << "Gauss ans: \n" << x1 << '\n' << "QR ans: \n" << x2 << std::endl;
//	auto b1 = A * x1;
//	auto residual1 = mismatch(b1, b);
//	std::cout << "Ans residual = " << residual1 << '\n' << std::endl;
//    
//	auto condV = cond(A);
//    std::cout << "Matrix cond = " << condV << std::endl;
//
//	auto b0 = matrix<T>("..\\labs3kurs\\tests\\test1_vector1.dat");
//	auto x01 = Gauss(A, b0);
//
//	auto ans_diff1 = mismatch(x01, x1);
//
//	auto Dx = 0.01;
//	auto dx1 = Dx / norm1(x1);
//	auto dxinf = Dx / norminf(x1);
//
//	auto Db = 0.01;
//	auto db1 = Db / norm1(b0);
//	auto dbinf = Db / norminf(x1);
//
//	auto cond_est1 = dx1 / db1;
//	auto cond_estinf = dxinf / dbinf;
//
//	std::cout << "Matrix cond estimation with norm1 = " << cond_est1;
//	std::cout << "\nMatrix cond estimation with norminf = " << cond_estinf << std::endl;
//
//	auto A1 = reverse(A);
//	auto E = A1 * A;
//	std::cout << "\nE: \n" << E << std::endl;
//
//	auto btest = b_mod(b);
//	std::cout << btest;*/
//	using T = double;
//	matrix_test<T>("..\\labs3kurs\\tests\\test3_matrix5.dat", "..\\labs3kurs\\tests\\test3_vector5.dat");
//    
//	return 0;
//}