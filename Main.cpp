//#include "Gaus.cpp"
//#include "Matrices_template.cpp"
#include "QR decomposition.cpp"
//#include "other.cpp"


int main()
{
	using T = double;
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

    matrix<T> A = matrix<T>("..//..//labs3kurs//matrix.dat");
    matrix<T> b = matrix<T>("..//..//labs3kurs//vector.dat");
    std::cout << "Matrix A: \n" << A << '\n' << "Vector b: \n" << b << std::endl;
    auto x1 = Gauss(A, b);
    auto x2 = QR_decomposition(A, b);

    std::cout << "Gauss ans: \n" << x1 << "QR and: \n" << x2 << std::endl;

	return 0;
}