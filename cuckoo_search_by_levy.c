#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include<chrono> 
//#include <algorithm>
#include <ctime>
#include <map>
#include <random>
#include <string>
#include <time.h>



using namespace std;

# define M_PI 3.14159265358979323846 

#define P 150
#define PA 0.25
#define BETA 1.5
#define TMAX 300
#define N 100
#define Down -5.12
#define High 5.12




double fitness(std::vector<double> X)
{
    double f = 0;
    if (N < 3)
    {
        double x = X[0];
        double y = X[1];
        //f = (x ** 2 + y - 11) ** 2 + (x + y ** 2 - 7) ** 2 #???????????
        // f = 0.26*(x**2 + y**2) - 0.48*x*y
        f = (x + 2 * y - 7) * (x + 2 * y - 7) + (2 * x + y - 5) * (2 * x + y - 5); //???? ????????

        return f;
    }
    else
    {
        for (int i = 0; i < N - 1; i++)
        {
            f += 100 * (X[i + 1] - X[i] * X[i]) * (X[i + 1] - X[i] * X[i]) + (X[i] - 1) * (X[i] - 1); //??????????
            //f += X[i][0] ** 2 #?????
            // f += 10*N + X[i][0]**2 - 10*np.cos(2*np.pi*X[i][0]) #??????????
        }
        return f;
    }
}

std::vector<double> Comparison(std::vector<double> A, std::vector<double> B)
{
    if (fitness(A) < fitness(B))
    {
        return A;
    }
    else
    {
        return B;
    }
}
std::vector<std::vector<double>> Create_by_Levy(std::vector<std::vector<double>>& X) {
    double gv = 1;
    double gu = pow(((tgamma(1 + BETA) * sin(M_PI * BETA / 2)) / tgamma((1 + BETA) / 2) * BETA * pow(2, (BETA - 1) / 2)), (1 / BETA));
    //cout << "gu " << gu;

    unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
    default_random_engine e(seed);

    std::vector<double> u(N);
    std::vector<double> v(N);
    normal_distribution<double> distN1(0, gu * gu);
    normal_distribution<double> distN2(0, gv * gv);

    for (int i = 0; i < N; i++)
    {
        u[i] = distN1(e);
        v[i] = distN2(e);

    }




    std::vector<double> s(N);
    for (int i = 0; i < N; i++) {
        s[i] = u[i] / pow(abs(v[i]), (1 / BETA));
    }


    std::vector<double> best(N);
    best = X[0];

    for (int i = 0; i < P; i++) {
        best = Comparison(best, X[i]);
    }

    std::vector<std::vector<double>> Xnew(P, std::vector<double>(N));
    for (int i = 0; i < P; i++) {
        Xnew[i] = X[i];
        for (int j = 0; j < N; j++) {
            Xnew[i][j] += 0.01 * s[j] * (Xnew[i][j] - best[j]);
            while (!(Xnew[i][j] < 5.12 || Xnew[i][j] > -5.12)) {
                Xnew[i][j] += 0.01 * s[j] * (Xnew[i][j] - best[j]);
            }
            X[i] = Comparison(Xnew[i], X[i]);
        }
    }
    return X;
}
std::vector<std::vector<double>> Replacement(std::vector<std::vector<double>> X) {
    std::vector<std::vector<double>> Xnew = X;
    for (int i = 0; i < P; i++) {
        int d1, d2;
        d1 = rand() % (P - 1);
        d2 = rand() % (P - 1);
        for (int j = 0; j < N; j++) {
            double r = rand() / RAND_MAX;
            if (r < PA) {
                Xnew[i][j] += (rand() / RAND_MAX) * (X[d1][j] - X[d2][j]);
                // while (Xnew[i, j] > 5.12 or Xnew[i, j] < -5.12):
                //     Xnew[i, j] += np.random.rand() * (X[d1, j] - X[d2, j])
            }
            X[i] = Comparison(X[i], Xnew[i]);
        }
    }
    return X;
}

std::vector<double> Find_best(std::vector<std::vector<double>> X) {
    // best = []
    std::vector<double> best = X[0];
    // best = X[0].copy();
    for (int i = 0; i < P; i++) {
        best = Comparison(best, X[i]);
    }
    return best;
}

int main() {
    //std::vector<std::vector<double>> X;
    std::vector<std::vector<double>> best;
    std::vector<double> f_best;

   
    std::vector<std::vector<double>> X(P, std::vector<double>(N));
    for (int i = 0; i < P; i++) {
        for (int j = 0; j < N; j++) {
            X[i][j] = rand();
        }
    }
    // X = new double* [P];
     //for (int i = 0; i < P; i++) {
       //  X[i] = new std::vector<double>[N];
     //}
     // X = np.random.randn(P, N)  # ??????? ?????? ????????? ????????? ???????
     // best = Find_best(X)

     //int start_time = time.perf_counter_ns();
    std::vector<double> a;



    double time_spent = 0.0;

    clock_t begin = clock();

    srand(time(NULL));

    for (int i = 0; i < TMAX; i++) {
        // print(X)
        X = Create_by_Levy(X); 
        X = Replacement(X); 

        best.push_back(Find_best(X));
        f_best.push_back(fitness(Find_best(X)));
        
        /*
        a = Find_best(X);
        for (int i = 0; i < N; i++) {

            cout << a[i] << "; ";
        }
        cout <<endl << fitness(a) << endl;

        */
    }


    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    cout << "time: " << time_spent << endl;
    // print("BEST: ", f_best);
    //long time_elapsed = time.perf_counter_ns() - start_time;
    //cout << "time: " << time_elapsed / 1_000_000_000 << endl;
    cout << "Best result: ";
    
    std::vector<double> c = Find_best(X);
    for (int i = 0; i < N; i++) {

        cout << c[i] << "; ";
    }

}
