#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>


double approx (const double x,const double p) {
    return p * round(x / p);
}

std::vector<double> approx_vector(const std::vector<double>& x, double p) {

    std::vector<double> approximation;
    // assicuro che approximation abbia la stessa dimensione di x
    approximation.resize(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        approximation[i] = approx(x[i],p);
    }

    return approximation;
}

using potential_function = double (*)(const std::vector<double>&);

double potential_harmonic_oscillator(const std::vector<double>& x) {
    return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
}

std::vector<double> fin_diff_gradient(potential_function potential, const std::vector<double>& x, const double dx) {

    std::vector<double> gradient;

    for (size_t i = 0; i < x.size(); i++) { // ciclo sugli elementi del gradiente
        // creo i due vettori che saranno i due punti di base con cui calcolare il rapporto incrementale uguali a x
        std::vector<double> x_plus_dx = x;
        std::vector<double> x_minus_dx = x;
        // li modifico aggiungendo il dx nell'entrata corrispodente in base all'elemento del gradiente
        x_plus_dx[i] += dx;
        x_minus_dx[i] -= dx;
        // calcolo il rapporto incrementale
        double finite_difference = (potential(x_plus_dx) - potential(x_minus_dx)) / (2 * dx);
        // aggiungo al gradiente
        gradient.push_back(finite_difference);
    }

    return approx_vector(gradient, dx);
}

void print_all_i_need(const std::vector<double>& x, const std::vector<double>& v, const double t) {
    std::cout << "x: ";
    for (double element : x) {
        std::cout << element << " ";
    }
    std::cout << "v: ";
    for (double element : v) {
        std::cout << element << " ";
    }
    std::cout << "t: " << t << "\n";
}

std::vector<double> prod(const std::vector<double>& x, const double s) {

    std::vector<double> product;

    for (double element : x) {
        product.push_back(element * s);
    }

    return product;
}

std::vector<double> sum(const std::vector<double>& x, const std::vector<double>& y) {

    std::vector<double> sum;

    // controllo se x e y hanno la stessa dimensione
    if (x.size() != y.size()) {
        // errore se hanno dimensioni diverse
        std::cerr << "Error: Vectors x and y must have the same size.\n";
        // Restituisci un vettore vuoto
        return sum;
    }

    for (size_t i = 0; i < x.size(); i++) {
        sum.push_back(x[i]+y[i]);
    }

    return sum;
}

void euler(std::vector<double>& x, std::vector<double>& v, potential_function potential, double dx, double dt, double t_max) {

    int n_step = t_max / dt;
    double t = 0;
    for (int i=0; i < n_step; i++) {
        t += dt;
        x = sum(x, prod(v,dt));
        v = sum(v, prod(fin_diff_gradient(potential_harmonic_oscillator, x, dx) , -dt));
        print_all_i_need(x,v,t);
    }
}

int main (int argc, char* argv[]){ //prendo i dati iniziali da terminale inserendoli dopo il lancio del programma nella stessa riga e.g.: ./euler 1.0 0.0 0.0 0.0 1.4142 0.0

   if (argc < 3 || (argc % 2) != 1) { // controllo di avere un numero pari di parametri
        std::cerr << "Usage: " << argv[0] << " value1_1 value1_2 ... valueN-1_1 valueN-1_2\n";
        return 1;  // ritorna con errore
    }

    // ottengo il numero di elementi
    int num_elements = (argc - 1) / 2;

    // dichiaro i vettori x e v
    std::vector<double> x_0(num_elements);
    std::vector<double> v_0(num_elements);

    // popolo i vettori x e v dagli argomenti dopo il nome del programma
    for (int i = 0; i < num_elements; ++i) {
        std::istringstream(argv[1 + i]) >> x_0[i];
        std::istringstream(argv[1 + num_elements + i]) >> v_0[i];
    }

    //std::vector<double> x_0 = {1.0, 0.0, 0.0};
    //std::vector<double> v_0 = {0.0, 1.4142, 0.0};
    //std::vector<double> x;
    //std::vector<double> v;
    double dx = 0.0000000001;
    double dt = 0.00001;
    double t = 0;

    print_all_i_need(x_0, v_0, t);

    euler(x_0, v_0, potential_harmonic_oscillator, dx, dt, 10); // seguo il punto materiale per 10 unita' di tempo
}
